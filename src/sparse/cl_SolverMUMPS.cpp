/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <cstdio>

#include "cl_SolverMUMPS.hpp"

#include "mumpstools.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "commtools.hpp"
#include "fn_available_memory.hpp"

namespace belfem
{
    namespace solver
    {
        namespace
        {
            // ICNTL(14) workspace relaxation in percent. The default is what
            // init_defaults() has always forced; the cap bounds the -9 retry
            // ladder ( 30 -> 60 -> 120 -> 240 -> 480, at most four doublings )
            // and with it the persist-bounded ratchet: an escalated value
            // stays for later solves on this instance, free() restores the
            // default. The ladder runs BEHIND the ICNTL(23) cap since
            // 2026-09-01: the guide says a -9 can still occur with the cap
            // set and still asks for a larger ICNTL(14) then
            constexpr int_t gDefaultMemoryRelaxation = 30 ;
            constexpr int_t gMaxMemoryRelaxation     = 480 ;

            // share of the available memory ( per node, per rank, per live
            // instance ) handed to MUMPS as ICNTL(23). The mesh, the other
            // kernel and the page cache all live in the rest
            constexpr double gMemoryBudgetSafety = 0.5 ;

            // the MUMPS shortfall convention for INFO(2)/INFOG(2): a positive
            // value is a count of entries, a negative one is in MILLIONS.
            // Never multiplied back in int_t -- that overflows int32 at
            // 2148 million
            std::string
            missing_entries( const int_t aInfo2 )
            {
                return aInfo2 < 0 ?
                    sprint( "%li million entries", ( long ) -aInfo2 ) :
                    sprint( "%li entries", ( long ) aInfo2 );
            }

            // the same convention for a bare size ( allocation sizes, file
            // sizes in the save/restore arms )
            std::string
            decode_count( const int_t aInfo2 )
            {
                return aInfo2 < 0 ?
                    sprint( "%li million", ( long ) -aInfo2 ) :
                    sprint( "%li", ( long ) aInfo2 );
            }
        }

//------------------------------------------------------------------------------

        MUMPS::MUMPS( const SolverParameters * aParams, const proc_t aMasterRank ) :
            Wrapper( "MUMPS    " , true ),
            mParams( aParams ),
            mMasterRank( aMasterRank )
        {
#ifdef BELFEM_MUMPS
            // allocate user settings
            mIParameters.set_size( 14, 0 );
            mRParameters.set_size( 1, 0.0 );

            // allocate information vectors ( INFO rank-local,
            // INFOG rank-uniform -- see the header ). 80 entries, the
            // library's own width: the BLR memory estimate is INFOG(36)
            mInfo.set_size( 80, 0 );
            mInfoG.set_size( 80, 0 );
            mRInfoG.set_size( 20, 0 );

            this->init_defaults();

#endif
        }

//------------------------------------------------------------------------------

        MUMPS::~MUMPS()
        {
            this->free();
        }

//------------------------------------------------------------------------------

        void
        MUMPS::init_defaults()
        {
#ifdef BELFEM_MUMPS

            // job
            mIParameters( static_cast< index_t >( mumps::Parameter::Job ) ) = -1 ;

            // rank of the Host
            mIParameters( static_cast< index_t >( mumps::Parameter::MasterRank ) ) = ( int_t ) mMasterRank ;

            // Host is working
            mIParameters( static_cast< index_t >( mumps::Parameter::WorkingHost ) ) = 1 ;

            switch ( mParams->reordering_method() )
            {
                // parmetis / ptscotch select what metis / scotch already
                // mean here: the parallel library above the serial one
                case( ReorderingMethod::METIS ):
                case( ReorderingMethod::PARMETIS ):
                {
                    mIParameters( static_cast< index_t >( mumps::Parameter::SerialReordering ) )
                        = static_cast< index_t >( MumpsSerialReodrdering::METIS ) ;
                    mIParameters( static_cast< index_t >( mumps::Parameter::ParallelReordering ) )
                        = static_cast< int_t >( MumpsParallelReodrdering::PARMETIS );
                    break;
                }
                case( ReorderingMethod::SCOTCH ) :
                case( ReorderingMethod::PTSCOTCH ) :
                {
                    mIParameters( static_cast< index_t >( mumps::Parameter::SerialReordering ) )
                        = static_cast< int_t >( MumpsSerialReodrdering::SCOTCH ) ;
                    mIParameters( static_cast< index_t >( mumps::Parameter::ParallelReordering ) )
                        = static_cast< int_t >( MumpsParallelReodrdering::PTSCOTCH );
                    break;
                }
                default:
                {
                    mIParameters( static_cast< index_t >( mumps::Parameter::SerialReordering ) )
                        = static_cast< int_t >( MumpsSerialReodrdering::AUTOMATIC ) ;
                    mIParameters( static_cast< index_t >( mumps::Parameter::ParallelReordering ) )
                        = static_cast< int_t >( MumpsParallelReodrdering::AUTOMATIC ) ;
                }
            }

            //    Determinant         : 0 - no
            //                        : 1 - yes
            mIParameters( static_cast< index_t >( mumps::Parameter::ComputeDeterminant )  ) = 0 ;

            // ICNTL(14): workspace relaxation in percent. The wrapper only
            // forwards this to MUMPS when > 0 ( the guarded write in
            // mumpstools_solve ), so the old 0 left MUMPS's own default
            // ( ~20-35% ) -- too small for this parallel factorization's
            // dynamic pivoting ( INFO(1) = -9, internal work array too
            // small ). A -9 now walks the escalate_workspace() ladder up to
            // gMaxMemoryRelaxation before the soft-fail arms see it
            mIParameters( static_cast< index_t >( mumps::Parameter::MemoryRelaxation )  )
                = gDefaultMemoryRelaxation ;

            // ICNTL(23): a deck-stated budget caps every factorization from
            // the first one on; unstated, the slot stays 0 and
            // escalate_workspace() fills it from the machine on the first
            // out-of-workspace failure
            mIParameters( static_cast< index_t >( mumps::Parameter::MemoryBudget ) )
                = mParams->have_memory_budget() ?
                    ( int_t ) mParams->memory_budget() : 0 ;

            // BLR only on explicit request, mirroring the STRUMPACK
            // grouping ( 2026-08-16 ): the old inverse switch sent
            // every value except OFF — including the AUTOMATIC default —
            // into BLR, and the unconditional CNTL(7) write below turned
            // that into silent LOSSY factorization at 1e-8 absolute on
            // every MUMPS deck without a compression key. BLR is an
            // accuracy-for-memory trade and must be a stated choice.
            // default covers OFF, AUTOMATIC and any future enumerator —
            // a new value must opt IN to compression, never fall into it
            switch ( mParams->compression_method() )
            {
                case( CompressionMethod::BLR ):
                {
                    // ICNTL(35) = Automatic ( 1 ), which MUMPS resolves to
                    // FactorizationAndSolution ( 2 ); pinning 2 explicitly
                    // is a driver-jury item, not this round
                    mIParameters( static_cast< index_t >( mumps::Parameter::CompressionMode  ) )
                        = static_cast< int_t >( MumpsBlockLowRanking::Automatic ) ;
                    break;
                }
                default:
                {
                    mIParameters( static_cast< index_t >( mumps::Parameter::CompressionMode  ) )
                        = static_cast< int_t >( MumpsBlockLowRanking::Off ) ;
                }
            }



            mRParameters( 0 )  = 0.0 ;
#endif
        }

//------------------------------------------------------------------------------

        void
        MUMPS::initialize( SpMatrix & aMatrix,
                const SymmetryMode    aSymmetryMode,
                const int_t aNumRhsColumns )
        {
#ifdef BELFEM_MUMPS
            // FIRST, before anything reserves a slot. A second initialize()
            // on a live wrapper used to create a new instance, overwrite
            // mSolverID with its id -- losing the first, which no longer had
            // an owner -- and only THEN reach Wrapper::initialize(), which
            // throws on its own already-initialized check. Two leaked slots
            // per misuse, one of them unreachable. Setup path, runs once per
            // instance, so this is BELFEM_ERROR and not an assert
            BELFEM_ERROR( ! this->is_initialized(),
                "MUMPS::initialize() called on a wrapper that already holds solver id %i - "
                "free() it before re-initializing",
                ( int ) mSolverID );

            // the parent's Wrapper::initialize() is NOT called here but after
            // the instance creation below succeeds: is_initialized() returns
            // the flag it sets, and Solver::solve keys its retry-or-not on
            // that flag. Setting it before the create would turn a failed
            // create into a wrapper that reports initialized, is never
            // re-initialized, and hands solver ID 0 to the Fortran side

            // insurance only: the shim always writes the out-arg. But an id
            // must never survive a create that did not happen, because an id
            // that outlives its slot names a slot the first-free scan may
            // since have handed to somebody else
            mSolverID = 0 ;

            // warn if the MPI ranks oversubscribe this node's cores with OpenMP threads
            this->hatch_turtle() ;

            // remember symmetry mode
            // ---- BELFEM CANNOT HONOUR A SYMMETRIC MUMPS MODE. -------------
            // Under SYM != 0, MUMPS wants exactly ONE representative of each
            // symmetric coordinate. Either triangle will do -- it does NOT
            // insist on the lower one, and an entry from the opposite triangle
            // is neither ignored nor rejected. What is fatal is supplying
            // BOTH: ( i, j ) and ( j, i ) are then treated as DUPLICATES and
            // SUMMED, exactly as duplicates are in unsymmetric mode.
            //
            // BELFEM stores and hands over the FULL matrix and nothing extracts
            // a triangle, so every off-diagonal arrives twice while the
            // diagonal arrives once. MUMPS factorizes A with doubled
            // off-diagonals -- a DIFFERENT matrix.
            //
            // It does not fail. It returns a converged, plausible, wrong
            // answer: measured 2026-08-29 on a 1D Laplacian with an analytic
            // spectrum, SYM = 2 gave a first-solve residual
            // ||Ax-b||/||b|| = 7.4e16 and an eigenvalue of -2.5e-19 against a
            // true 2.46e-6, with the eigensolver reporting success throughout.
            // SYM = 0 on the same matrix: residual 1.3e-12, eigenvalue right
            // to 6e-12.
            //
            // ( MUMPS 5.5.1 user guide section 5.2.2.1. An earlier version of
            //   this comment said MUMPS "requires the lower triangle". That
            //   was wrong, and the distinction matters to whoever implements
            //   the extraction: they may pick either triangle. )
            //
            // Always-active, and an error rather than a silent downgrade to
            // SYM = 0: a caller that asked for symmetry wants the memory
            // saving, and quietly not giving it is its own kind of lie. To
            // support this properly, supply one triangle -- see
            // todo/mumps_symmetric_triangle_extraction.md
            BELFEM_ERROR( aSymmetryMode == SymmetryMode::Unsymmetric,
                "MUMPS symmetry mode %i is not supported: for SYM != 0 MUMPS wants ONE "
                "representative per symmetric coordinate and sums ( i, j ) with ( j, i ), "
                "while BELFEM supplies the full matrix - every off-diagonal would be "
                "counted twice. Use SymmetryMode::Unsymmetric, or add triangle extraction "
                "to this wrapper.",
                ( int ) aSymmetryMode );

            mSymmetryMode = aSymmetryMode ;

            // symmetry mode
            mIParameters( static_cast< index_t >( mumps::Parameter::SymmetryMode )  )
                =  static_cast< int_t >( aSymmetryMode ) ;


            // info level
            mIParameters(  static_cast< index_t >( mumps::Parameter::InfoLevel ) )
                = ( int_t ) gLog.info_level() ;


            // 5: Refinement          :<0 - fixed number of steps
            //                        : 0 - none
            //                        :>0 - maximum number of refinement steps
            // A positive value caps iterative refinement with a convergence test, so it
            // raises INFO(1)=8 if the scaled residual is not met within the cap ( common
            // on ill-conditioned HTS systems, where refinement plateaus ). Use a negative
            // value here for a fixed number of steps with no convergence test ( no warning )
            // if the +8 message becomes noise.
            mIParameters( static_cast< index_t >( mumps::Parameter::NumRefinementSteps ) )
                = aNumRhsColumns == 1 ? 20 : 0 ;


            // BLR dropping parameter — this slot is CNTL(7), NOT a
            // refinement tolerance ( the old label hid that ). Flows only
            // on explicit blr; otherwise 0.0, which MUMPS documents as
            // "full precision" ( 5.7.3 §5.19 ). With BLR off the value is
            // inert either way — the conditional is defense-in-depth and
            // honest self-documentation, not a behavior carrier
            mRParameters( 0 ) =
                mParams->compression_method() == CompressionMethod::BLR ?
                    mParams->compression_cutoff() : 0.0 ;

            int_t tInfo ;
            mumpstools_create_solver(
                    mSolverID,
                    tInfo,
                    mIParameters( static_cast< index_t >( mumps::Parameter::WorkingHost ) ),
                    mIParameters( static_cast< index_t >( mumps::Parameter::SymmetryMode ) ) );

            // BOTH halves are checked. A positive id alone is not proof of a
            // usable instance: the shim reports the JOB = -1 outcome through
            // tInfo, and a create that reserved a slot and then failed inside
            // DMUMPS would otherwise walk straight through an id-only gate.
            // The shim rolls the occupancy back and returns id -1 in exactly
            // that case, so the tInfo arm is defence in depth against the two
            // sides drifting apart -- not the only line. tInfo > 0 is a MUMPS
            // WARNING and must not abort the create
            if( mSolverID <= 0 || tInfo < 0 )
            {
                // registry exhausted, or the instance creation itself failed.
                // For a diagnostic consumer this must degrade, not abort: the
                // wrapper stays UN-initialized, so Solver::solve retries the
                // create on the next call, when a slot may have been freed.
                // The pool state is per process, and this branch is
                // rank-uniform only because every in-tree caller runs the same
                // collective create/free sequence -- NOT because anything
                // enforces it. If ranks ever disagreed about which slot is
                // free, one would return here while another entered DMUMPS on
                // MPI_COMM_WORLD; the shim says exactly that at its exhaustion
                // arm, and this side must not restate it as a guarantee
                //
                // the id is dropped either way: the shim owns the rollback of
                // the slot, and keeping a number that names a slot we do not
                // own is how a later free destroys somebody else's instance
                mSolverID = 0 ;

                if( this->soft_fail() )
                {
                    // print_soft_fail() reads mInfoG( 0 ); both arrays are
                    // stamped so neither reports whatever the last solve
                    // left behind. tInfo here is the shim's rank-local
                    // INFO(1) ( or -1000 pool exhaustion ) -- create is
                    // deliberately NOT on the INFOG export, see the audit
                    // note in mumpstools.f90
                    mInfo.fill( 0 );
                    mInfo( 0 ) = tInfo ;
                    mInfoG.fill( 0 );
                    mInfoG( 0 ) = tInfo ;

                    this->flag_failure() ;
                    this->print_soft_fail() ;
                    return ;
                }

                // -1000 is the shim's own code, not a MUMPS INFO value, and
                // error_message() has no arm for it -- only print_soft_fail()
                // does, which this branch does not reach. Decode it here or a
                // production run ( soft-fail off ) reads a bare -1000
                BELFEM_ERROR( false,
                    "couldn't initialize MUMPS solver: %s ( INFO(1) = %i )",
                    tInfo == -1000 ? "the instance registry had no free slot"
                                   : "the instance creation failed",
                    ( int ) tInfo );
            }

            // the parent's flag is set only now that an instance exists:
            // is_initialized() must never say yes for a wrapper that holds
            // no solver ID
            Wrapper::initialize();

            // ... and so is ours, AFTER the call above, so that a throw in
            // there cannot leave this flag true against a false mIsInitialized.
            //
            // mInitialized means exactly one thing: a Fortran instance exists
            // that free() owes a JOB = -2. That becomes true HERE, when the
            // create succeeds -- not at the first solve, which is where the
            // three old write sites were. An instance freed before it ever
            // solved, or used only through the matrix-RHS overload, used to
            // miss its teardown entirely and leak its slot
            mInitialized = true ;

            // measure the machine NOW, before this instance has allocated a
            // single factor: ICNTL(23) bounds the total an instance may
            // take, and a probe after a failed factorization would see the
            // remnant next to the workspace MUMPS still holds ( audit
            // finding ). Collective, like the create above; the create's
            // soft-fail return is taken uniformly, so no rank is missing
            mMeasuredBudgetMB = this->memory_budget_mb();

#endif
        }

//------------------------------------------------------------------------------

        void
        MUMPS::free()
        {
#ifdef BELFEM_MUMPS
            // JOB -2 is what actually releases the factors, so the frozen
            // scope cannot outlive this call. Dropped BEFORE the release, so
            // that a failure inside the teardown cannot leave the scope armed
            // over factors that no longer exist
            this->invalidate_factorization() ;

            if ( mInitialized )
            {
                mumpstools_free_solver( mSolverID, mInfo.data() );

                // check result
                if( mInfo( 0 ) != 0 && mMatrix != nullptr )
                {
                    // INFO(1) < 0 is an error, > 0 is a warning ( see the MUMPS manual ).
                    if( mInfo( 0 ) < 0 && this->rank() == 0 )
                    {
                        std::string tMessage = this->error_message(
                                mInfo.data(),
                                mMatrix->n_rows(),
                                mMatrix->number_of_nonzeros() );

                        BELFEM_ERROR( false,
                               "MUMPS has thrown the error: %i\n%s",
                                mInfo( 0 ),
                                tMessage.c_str() );
                    }
                    else if( mInfo( 0 ) > 0 )
                    {
                        // warnings are RANK-LOCAL: unlike errors, MUMPS does not
                        // propagate them, so reading INFO on rank 0 alone loses
                        // every warning raised only on a worker. Each rank reports
                        // its own, tagged, so the message names where it happened
                        //
                        // deliberately NOT check_warnings(): free() is teardown,
                        // reached from ~MUMPS() after JOB = -2, and a hard error
                        // here would abort during cleanup and mask whatever
                        // caused the teardown. This site stays report-only
                        Cell< string > tWarnings;
                        this->warning_message( mInfo.data(), tWarnings );

                        for( const string & tWarning : tWarnings )
                        {
                            message( InfoLevel::Minimal,
                                     "MUMPS returned a warning on proc %i ( INFO(1) = %i ): %s",
                                     ( int ) this->rank(),
                                     ( int ) mInfo( 0 ),
                                     tWarning.c_str() );
                        }
                    }
                }

                mInitialized = false ;

                // the slot is no longer ours. Keeping its number would let a
                // second free -- and there IS one, ~MUMPS() runs free() again
                // after Solver::free() -- issue JOB = -2 against whatever
                // instance the first-free scan has since put in that slot.
                // gOccupied would read 1 for the NEW tenant, so the shim would
                // not skip it: a cross-wrapper destroy, not a benign double
                // free. mInitialized alone happens to prevent it today; the id
                // is cleared so that it does not have to
                mSolverID = 0 ;
            }
#endif
            // ---- AFTER the teardown reporting above, which reads mMatrix to
            // build its error message. Clearing it earlier made the
            // `mInfo( 0 ) != 0 && mMatrix != nullptr` guard permanently false
            // and silenced every JOB -2 error and warning -- a regression
            // introduced and caught in the same session, 2026-08-29.
            //
            // It does have to be cleared, though: mMatrix is the "which matrix
            // are these factors for" record, and JOB -2 has just destroyed the
            // instance holding them. Leaving it set would make the NEXT solve
            // on the same SpMatrix select JOB 5 ( factorize + solve, reusing an
            // analysis ) against a fresh instance that has analysed nothing
            mMatrix = nullptr ;

#ifdef BELFEM_MUMPS
            // the -9 ladder's escalation is persist-bounded: it survives for
            // later solves on THIS instance, but must not leak into the next
            // one -- init_defaults() runs only in the constructor, so without
            // this write a torn-down wrapper would hand its escalated value
            // ( up to 480 ) to whatever initialize() comes next. Outside the
            // mInitialized gate ( the C++ member escalates even when the
            // Fortran instance later fails ) but inside the preprocessor
            // gate: without BELFEM_MUMPS the constructor never sizes
            // mIParameters, and the destructor reaches this line ( audit
            // finding, both auditors )
            mIParameters( static_cast< index_t >( mumps::Parameter::MemoryRelaxation )  )
                = gDefaultMemoryRelaxation ;

            // the measured ICNTL(23) cap is per instance for the same
            // reason; a deck-stated one is the deck's and stays
            mIParameters( static_cast< index_t >( mumps::Parameter::MemoryBudget ) )
                = mParams->have_memory_budget() ?
                    ( int_t ) mParams->memory_budget() : 0 ;
#endif

            // the soft-failure latch belongs to the instance being torn down,
            // not to its successor
            this->clear_failure() ;

            // call function from parent, nothing else to be done here
            Wrapper::free();
        }

//------------------------------------------------------------------------------

        bool
        MUMPS::supports_factorization_reuse() const
        {
            return true ;
        }

//------------------------------------------------------------------------------

        void
        MUMPS::freeze_factorization( const SpMatrix & aMatrix )
        {
            // arming without factors is the one mistake that would produce a
            // JOB 3 against nothing. mMatrix is set by a successful master
            // solve and cleared by every failure path, so on the master it is
            // exactly the "we hold valid factors" flag
            BELFEM_ERROR( this->is_initialized(),
                "freeze_factorization(): MUMPS is not initialized - "
                "a factorization must succeed before the frozen scope is armed" );

            mFrozen.mArmed = true ;

            // Only the rank that owns the matrix can describe it. On a worker
            // aMatrix is a local submatrix that the distributed solve never
            // reads ( the shim attaches matrix pointers on the MUMPS master
            // alone ), so recording its shape would describe the wrong object
            // and the later comparison would be meaningless
            if( this->rank() == mMasterRank )
            {
                BELFEM_ERROR( mMatrix == &aMatrix,
                    "freeze_factorization(): the matrix handed in is not the one "
                    "that was factorized - arm the scope immediately after the "
                    "solve that built the factors" );

                mFrozen.mHaveIdentity = true ;
                mFrozen.mSource       = &aMatrix ;
                mFrozen.mNumRows      = aMatrix.n_rows() ;
                mFrozen.mNumCols      = aMatrix.n_cols() ;
                mFrozen.mNumNonZeros  = aMatrix.number_of_nonzeros() ;
                mFrozen.mValues       = aMatrix.data() ;
                mFrozen.mPointers     = aMatrix.pointers() ;
            }
            else
            {
                mFrozen.mHaveIdentity = false ;
            }
        }

//------------------------------------------------------------------------------

        void
        MUMPS::unfreeze_factorization()
        {
            mFrozen = FrozenFactorization() ;
        }

//------------------------------------------------------------------------------

        bool
        MUMPS::factorization_is_frozen() const
        {
            return mFrozen.mArmed ;
        }

//------------------------------------------------------------------------------

        void
        MUMPS::invalidate_factorization()
        {
            // the factors are gone or untrustworthy, so the scope cannot
            // survive even if the caller has not left it yet. Silent: this
            // runs on failure paths that are already reporting
            mFrozen = FrozenFactorization() ;
        }

//------------------------------------------------------------------------------

        bool
        MUMPS::escalate_workspace( int_t & aJob )
        {
            // keyed on the rank-uniform INFOG: the MUMPS user guide
            // guarantees INFOG(1:2) carry the SAME code and supplementary
            // value on every rank ( unlike INFO, where only the failing rank
            // holds the true code ), so every rank takes the same branch and
            // the collective re-entry cannot strand anyone. -9: real
            // workarray S too small; -8: integer workarray IS too small;
            // -17 / -20: MPI send / reception buffer too small ( all four
            // ask for a larger ICNTL(14) per the guide ); -19: the ICNTL(23)
            // cap cannot be met
            const int_t tCode = mInfoG( 0 );

            if(    tCode != -9 && tCode != -8 && tCode != -17 && tCode != -20
                && tCode != -19 )
            {
                return false ;
            }

            // a frozen JOB 3 never retries. It also cannot legally raise a
            // -9: workspace exhaustion during the solve phase is -11, and
            // -9 belongs to the factorization -- if it appears here anyway,
            // hand it to the soft-fail arms rather than guessing
            if( aJob != 5 && aJob != 6 )
            {
                return false ;
            }

            const index_t tRelax
                = static_cast< index_t >( mumps::Parameter::MemoryRelaxation ) ;
            const index_t tBudget
                = static_cast< index_t >( mumps::Parameter::MemoryBudget ) ;

            // hang guard: a relaxation below the default cannot come from
            // this class today, but 0 doubled stays 0, the shim's .gt.-0
            // write would never fire, and the ladder would re-enter the
            // collective forever ( audit finding )
            if( mIParameters( tRelax ) < gDefaultMemoryRelaxation )
            {
                return false ;
            }

            // the machine's budget was measured once, in initialize(),
            // before MUMPS held anything; rank-uniform by construction
            const int_t tBudgetMB = mMeasuredBudgetMB ;

            // MUMPS's own estimate of the factorization's working space in
            // MB, from the analysis ( guide sec. 5.12 ): INFOG(16) for
            // full-rank factors, INFOG(36) under BLR ( ICNTL(35) = 1 or 2 ).
            // 0-based here
            const int_t tBlr = mIParameters(
                static_cast< index_t >( mumps::Parameter::CompressionMode ) );
            const int_t tBoundMB = ( tBlr == 1 || tBlr == 2 ) ?
                mInfoG( 35 ) : mInfoG( 15 );

            const mumps::WorkspaceAction tAction = mumps::next_workspace_action(
                    tCode,
                    mIParameters( tRelax ),
                    gMaxMemoryRelaxation,
                    mIParameters( tBudget ),
                    tBudgetMB,
                    tBoundMB );

            // sized for the worst-case decimal widths; the print below
            // clips to the 71-char box field, so the frame cannot be punched
            char tRow[ 128 ];

            switch ( tAction )
            {
                case mumps::WorkspaceAction::Cap :
                {
                    // EVERY rank writes the slot: the shim sets ICNTL(23)
                    // from its local parameter copy, and the value is the
                    // reduced one, so MUMPS reads the same cap everywhere
                    mIParameters( tBudget ) = tBudgetMB ;

                    std::snprintf( tRow, sizeof( tRow ),
                        " MUMPS out of workspace ( %i ): cap ICNTL(23) = %li MB, est. %li MB",
                        ( int ) tCode, ( long ) tBudgetMB, ( long ) tBoundMB );
                    break ;
                }
                case mumps::WorkspaceAction::Ladder :
                {
                    // EVERY rank escalates, for the same reason. Clamped so
                    // no writer of this slot can overshoot the cap. The
                    // escalation persists for later solves on this instance
                    // ( a hard phase clusters, so repaying the failed
                    // factorizations every step would be waste ); free()
                    // restores the default
                    mIParameters( tRelax ) *= 2 ;

                    if( mIParameters( tRelax ) > gMaxMemoryRelaxation )
                    {
                        mIParameters( tRelax ) = gMaxMemoryRelaxation ;
                    }

                    std::snprintf( tRow, sizeof( tRow ),
                        " MUMPS out of workspace ( %i ): retrying with ICNTL(14) = %i",
                        ( int ) tCode, ( int ) mIParameters( tRelax ) );
                    break ;
                }
                default :
                {
                    return false ;
                }
            }

            // repeat factorization + solve against the analysis the failed
            // attempt completed: MUMPS prescribes raising ICNTL(14) and
            // calling the factorization again for -9, and both ICNTL(14)
            // and ICNTL(23) are read at the start of the factorization
            // phase, so the change takes effect without a re-analysis. JOB
            // 2 alone would be wrong through this shim -- it copies the RHS
            // into the solution slot before the call, so a factorize-only
            // job would return the RHS as "solution". If MUMPS rejects the
            // JOB 5 after all ( -3 ), the next pass breaks on the code check
            // above and the soft-fail arms see it -- no hang
            aJob = 5 ;

            if( this->rank() == 0 )
            {
                // the retry runs INSIDE the caller's log box, so the report
                // is drawn as a box row ( inner width 71, the same field
                // print_soft_fail() writes into ). A bare message here
                // punched a hole through the frame on every escalation
                message( InfoLevel::Minimal, "   │%-71.71s│", tRow );
            }

            return true ;
        }

//------------------------------------------------------------------------------

        int_t
        MUMPS::memory_budget_mb()
        {
            int_t tLocal = 0 ;
#ifdef BELFEM_MUMPS
            const std::size_t tBytes = available_memory();

            int_t tInstances = 0 ;
            mumpstools_num_solvers( tInstances );

            const proc_t tRanksOnNode = gComm.node_size() > 0 ? gComm.node_size() : 1 ;

            if ( tBytes > 0 )
            {
                // MUMPS counts MegaBytes as 10^6 bytes ( guide sec. 5.12 )
                tLocal = ( int_t ) ( gMemoryBudgetSafety * ( double ) tBytes
                                     / ( double ) tRanksOnNode
                                     / ( double ) ( tInstances > 0 ? tInstances : 1 )
                                     / 1.0e6 );

                // a machine that was measured and has nothing left must not
                // round down into "unknown", which would send the retry up
                // the ladder to ask for more ( audit finding ): 1 MB is a
                // known budget that the estimate check turns into give-up
                if ( tLocal < 1 )
                {
                    tLocal = 1 ;
                }
            }
#endif
            // the most constrained rank decides, and a rank that could not
            // measure ( 0 ) pulls everyone to "unknown" -- the ladder --
            // rather than being outvoted into a cap it cannot honour
            int_t tGlobal = 0 ;
            allreduce_min( &tLocal, &tGlobal, 1 );

            return tGlobal ;
        }

//------------------------------------------------------------------------------

        int_t
        MUMPS::select_job( SpMatrix & aMatrix )
        {
            if( mFrozen.mArmed )
            {
                // everything we can check without an O( nnz ) scan. A values
                // CHANGE is not detectable here at any acceptable cost -- that
                // is the caller's half of the contract ( no assembly inside
                // the scope ) -- but a different matrix, a resized one, or one
                // whose storage has been reallocated underneath us all show up
                // in these pointers and counts.
                //
                // Always-active, not BELFEM_ASSERT: solving with stale factors
                // returns a plausible answer to the wrong question, and that
                // must not be a release-build behaviour
                if( mFrozen.mHaveIdentity )
                {
                    BELFEM_ERROR( mFrozen.mSource == &aMatrix,
                        "frozen factorization: solve called with a different matrix "
                        "than the one that was factorized" );

                    BELFEM_ERROR( mFrozen.mNumRows == aMatrix.n_rows()
                               && mFrozen.mNumCols == aMatrix.n_cols()
                               && mFrozen.mNumNonZeros == aMatrix.number_of_nonzeros(),
                        "frozen factorization: the matrix changed shape while frozen "
                        "( was %lu x %lu with %lu nonzeros, is %lu x %lu with %lu )",
                        ( long unsigned int ) mFrozen.mNumRows,
                        ( long unsigned int ) mFrozen.mNumCols,
                        ( long unsigned int ) mFrozen.mNumNonZeros,
                        ( long unsigned int ) aMatrix.n_rows(),
                        ( long unsigned int ) aMatrix.n_cols(),
                        ( long unsigned int ) aMatrix.number_of_nonzeros() );

                    BELFEM_ERROR( mFrozen.mValues == aMatrix.data()
                               && mFrozen.mPointers == aMatrix.pointers(),
                        "frozen factorization: the matrix storage was reallocated "
                        "while frozen" );
                }

                // JOB 3: solve only, against the factors already in the saved
                // DMUMPS_STRUC. They survive any number of JOB 3 calls and are
                // released only by JOB -2
                return 3 ;
            }

            // unfrozen: exactly the policy that was here before, unchanged
            if ( mMatrix != &aMatrix ) // matrix has not been initialized or has changed
            {
                mMatrix = &aMatrix;
                return 6 ; // initialize
            }

            return 5 ; // same matrix, just different values
        }

//------------------------------------------------------------------------------

        void
        MUMPS::solve(            SpMatrix & aMatrix,
                           Vector< real > & aLHS,
                           Vector< real > & aRHS )
        {
#ifdef BELFEM_MUMPS

            // a soft-failed create left no MUMPS instance. Bail flagged
            // before anything hands solver ID 0 to the Fortran registry --
            // the is_initialized() assert below cannot catch this in a
            // release build, and the branch is rank-uniform ( the registry
            // marches in lockstep on every rank ), so returning here cannot
            // desynchronize the ranks. The silent path exists ONLY for the
            // caller that asked for it: without soft-fail this is caller
            // misuse and must be loud
            if( mSolverID <= 0 )
            {
                if( ! this->soft_fail() )
                {
                    BELFEM_ERROR( false,
                        "MUMPS::solve() called without a solver instance ( initialize failed or was never called )" );
                }
                else
                {
                    this->flag_failure() ;
                    return ;
                }
            }

            // make sure that the wrapper has been initialited
            BELFEM_ASSERT( this->is_initialized(),
                "MUMPS has not been initialited" );

            // reset info vectors
            mInfo.fill( 0 );
            mInfoG.fill( 0 );

            if( aLHS.length() != aRHS.length() )
            {
                aLHS.set_size( aRHS.length(), 0.0 );
            }

            int_t tJob = 0 ;

            if( this->rank() == mMasterRank )
            {

                // make sure that all indices have been created
                aMatrix.create_coo_indices() ;

                // make sure that matrix is stored one-based. Restored AFTER
                // the retry ladder below -- restoring between attempts would
                // hand 0-based indices to the retry
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );

                // shared with the multiple-RHS overload below, so the two
                // cannot drift apart: JOB 3 while the frozen scope is armed,
                // otherwise the unchanged 6-for-new / 5-for-repeat policy
                tJob = this->select_job( aMatrix );

                broadcast( tJob );
            }
            else
            {
                broadcast( tJob );
            }

            // identical on every rank, so set once outside the rank branch
            mIParameters( static_cast< index_t >( mumps::Parameter::SymmetryMode )  )
                =  static_cast< int_t >( mSymmetryMode ) ;

            mIParameters( static_cast< index_t >( mumps::Parameter::SolverID ) ) = mSolverID ;

            // workspace retry. MUMPS raises -9 ( -8 ) when dynamic pivoting
            // outgrows the real ( integer ) workspace sized from the
            // analysis estimate; the recovery is to hand it the machine's
            // budget as ICNTL(23) and, behind that, to raise ICNTL(14), and
            // repeat the FACTORIZATION ( MUMPS user guide, errors -9 / -8
            // and sec. 2.11 -- both controls are read at the start of the
            // factorization, so no re-analysis ). Every decision in
            // escalate_workspace() is taken from the rank-uniform INFOG, so
            // all ranks leave or re-enter the collective call together. A
            // failure that survives the policy falls through to the
            // soft-fail arms below unchanged
            while( true )
            {
                mIParameters( static_cast< index_t >( mumps::Parameter::Job ) ) = tJob ;

                if( this->rank() == mMasterRank )
                {
                    // solve the system
                    mumpstools_solve(
                            mIParameters.data(),
                            mRParameters.data(),
                            aMatrix.n_rows(),
                            aMatrix.number_of_nonzeros(),
                            1,
                            aMatrix.rows(),
                            aMatrix.cols(),
                            aMatrix.data(),
                            aLHS.data(),
                            aRHS.data(),
                            mInfo.data(),
                            mInfoG.data(),
                            mRInfoG.data() );
                }
                else
                {
                    // solve the system as slave
                    mumpstools_solve(
                            mIParameters.data(),
                            mRParameters.data(),
                            0,
                            0,
                            1,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            mInfo.data(),
                            mInfoG.data(),
                            mRInfoG.data() );
                }

                if( ! this->escalate_workspace( tJob ) )
                {
                    break ;
                }
            }

            if( this->rank() == mMasterRank )
            {
                // restore C++ indexing for downstream consumers
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );
            }

            // check result. MUMPS convention: INFO(1)/INFOG(1) < 0 is an
            // error, > 0 a warning ( a usable solution is still returned ).
            // The error arm keys on INFOG: not just its SIGN but the CODE
            // and supplementary value are rank-uniform, so every rank takes
            // the same branch and decodes the true failure instead of the
            // propagated -1 "error on rank N". A rank-0-only throw once left
            // the other ranks marching into the next collective and the job
            // died by segfault instead of the error box ( observed
            // 2026-07-27, INFOG(1) = -10 quench )
            if( mInfoG( 0 ) < 0 )
            {
                // soft-fail contract: record and return uniformly on all
                // ranks ( the branch keys on the rank-uniform INFOG, and a
                // -9 only reaches here after the escalate_workspace() ladder
                // gave up ); the controller treats the event
                // as a failed trial and cuts the timestep. The shim instance
                // stays alive and is REUSED by the next attempt ( JOB 5/6 on
                // the same solver id ).
                //
                // mInitialized is NOT set here any more: initialize() sets it
                // when the instance is created, which is when free() starts
                // owing a JOB = -2. Setting it at the first solve was what
                // made an unsolved -- or matrix-RHS-only -- instance leak
                if ( this->soft_fail() )
                {
                    this->flag_failure() ;

                    // force a full JOB=6 ( analysis + factorization + solve )
                    // on the next attempt: a failure during analysis leaves
                    // no valid analysis for a JOB=5 reuse ( Codex round-5 )
                    mMatrix = nullptr ;

                    // and the factors this failure leaves behind, if any, are
                    // not ones a frozen JOB 3 may solve against
                    this->invalidate_factorization() ;

                    this->print_soft_fail() ;
                    return ;
                }

                std::string tMessage = this->error_message(
                        mInfoG.data(),
                        aMatrix.n_rows(),
                        aMatrix.number_of_nonzeros() );

                BELFEM_ERROR( false,
                       "MUMPS has thrown the error: %i\n%s",
                        mInfoG( 0 ),
                        tMessage.c_str() );
            }
            else if( mInfo( 0 ) > 0 )
            {
                // RANK-LOCAL, unlike the error path above: MUMPS propagates
                // errors ( non-failing ranks see INFO(1) = -1 ), but not
                // warnings, so this must run on every rank. The out-of-range
                // index bit ( +1 ) is escalated to a hard error in there
                this->check_warnings();
            }

#else
            BELFEM_ERROR( false, "We are not linked against MUMPS" );
#endif
        }

//------------------------------------------------------------------------------

        void
        MUMPS::solve(
                SpMatrix & aMatrix,
                Matrix< real > & aLHS,
                Matrix< real > & aRHS )
        {
#ifdef BELFEM_MUMPS

            // same guard as the vector overload: no instance, no Fortran
            // call -- and by the same rule, the silent return is only
            // for a caller that armed soft-fail
            if( mSolverID <= 0 )
            {
                if( ! this->soft_fail() )
                {
                    BELFEM_ERROR( false,
                        "MUMPS::solve() called without a solver instance ( initialize failed or was never called )" );
                }
                else
                {
                    this->flag_failure() ;
                    return ;
                }
            }

            // reset info vectors
            mInfo.fill( 0 );
            mInfoG.fill( 0 );

            // allocate matrix sizes of not set
            if(  aLHS.n_rows() != aRHS.n_rows() || aLHS.n_cols() != aRHS.n_cols() )
            {
                aLHS.set_size( aRHS.n_rows(), aRHS.n_cols(), 0.0 );
            }

            int_t tJob = 0 ;

#ifndef BELFEM_ARMADILLO
            // flattened copies for the Blaze path, hoisted above the retry
            // ladder so a retry does not reallocate or refill them ( the
            // shim re-copies aY into aX on every attempt itself )
            Vector< real > tX;
            Vector< real > tY;
#endif

            if( this->rank() == mMasterRank )
            {
                // make sure that all indices have been created
                aMatrix.create_coo_indices() ;

                // make sure that matrix is stored one-based. Restored AFTER
                // the retry ladder below, as in the vector overload
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );

                // same selector as the vector overload -- a reuse policy that
                // lived in only one of the two would be a silent inconsistency
                tJob = this->select_job( aMatrix );

                broadcast( tJob );

#ifndef BELFEM_ARMADILLO
                this->mat2vec( aRHS, tY );

                // allocate memory for LHS
                tX.set_size( tY.length() );
#endif
            }
            else
            {
                broadcast( tJob );
            }

            mIParameters( static_cast< index_t >( mumps::Parameter::SolverID ) ) = mSolverID ;

            // same workspace retry as the vector overload -- the policy
            // lives in escalate_workspace() so the two cannot drift apart
            while( true )
            {
                mIParameters( static_cast< index_t >( mumps::Parameter::Job ) ) = tJob ;

                if( this->rank() == mMasterRank )
                {
#ifdef BELFEM_ARMADILLO
                    // solve the system
                    mumpstools_solve(
                            mIParameters.data(),
                            mRParameters.data(),
                            aMatrix.n_rows(),
                            aMatrix.number_of_nonzeros(),
                            aRHS.n_cols(),
                            aMatrix.rows(),
                            aMatrix.cols(),
                            aMatrix.data(),
                            aLHS.data(),
                            aRHS.data(),
                            mInfo.data(),
                            mInfoG.data(),
                            mRInfoG.data() );
#else
                    // solve the system
                    mumpstools_solve(
                            mIParameters.data(),
                            mRParameters.data(),
                            aMatrix.n_rows(),
                            aMatrix.number_of_nonzeros(),
                            aRHS.n_cols(),
                            aMatrix.rows(),
                            aMatrix.cols(),
                            aMatrix.data(),
                            tX.data(),
                            tY.data(),
                            mInfo.data(),
                            mInfoG.data(),
                            mRInfoG.data() );
#endif
                }
                else
                {
                    // solve the system as slave
                    mumpstools_solve(
                            mIParameters.data(),
                            mRParameters.data(),
                            0,
                            0,
                            1,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            mInfo.data(),
                            mInfoG.data(),
                            mRInfoG.data() );
                }

                if( ! this->escalate_workspace( tJob ) )
                {
                    break ;
                }
            }

            if( this->rank() == mMasterRank )
            {
#ifndef BELFEM_ARMADILLO
                // unflatten vector to matrix
                this->vec2mat( tX, aLHS );
#endif
                // restore C++ indexing for downstream consumers
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );
            }

            // check result. MUMPS convention: INFO(1)/INFOG(1) < 0 is an
            // error, > 0 a warning ( a usable solution is still returned ).
            // The error arm keys on INFOG: not just its SIGN but the CODE
            // and supplementary value are rank-uniform, so every rank takes
            // the same branch and decodes the true failure instead of the
            // propagated -1 "error on rank N". A rank-0-only throw once left
            // the other ranks marching into the next collective and the job
            // died by segfault instead of the error box ( observed
            // 2026-07-27, INFOG(1) = -10 quench )
            if( mInfoG( 0 ) < 0 )
            {
                // soft-fail contract: record and return uniformly on all
                // ranks ( the branch keys on the rank-uniform INFOG, and a
                // -9 only reaches here after the escalate_workspace() ladder
                // gave up ); the controller treats the event
                // as a failed trial and cuts the timestep. The shim instance
                // stays alive and is REUSED by the next attempt ( JOB 5/6 on
                // the same solver id ).
                //
                // mInitialized is NOT set here any more: initialize() sets it
                // when the instance is created, which is when free() starts
                // owing a JOB = -2. Setting it at the first solve was what
                // made an unsolved -- or matrix-RHS-only -- instance leak
                if ( this->soft_fail() )
                {
                    this->flag_failure() ;

                    // force a full JOB=6 ( analysis + factorization + solve )
                    // on the next attempt: a failure during analysis leaves
                    // no valid analysis for a JOB=5 reuse ( Codex round-5 )
                    mMatrix = nullptr ;

                    // and the factors this failure leaves behind, if any, are
                    // not ones a frozen JOB 3 may solve against
                    this->invalidate_factorization() ;

                    this->print_soft_fail() ;
                    return ;
                }

                std::string tMessage = this->error_message(
                        mInfoG.data(),
                        aMatrix.n_rows(),
                        aMatrix.number_of_nonzeros() );

                BELFEM_ERROR( false,
                       "MUMPS has thrown the error: %i\n%s",
                        mInfoG( 0 ),
                        tMessage.c_str() );
            }
            else if( mInfo( 0 ) > 0 )
            {
                // RANK-LOCAL, unlike the error path above: MUMPS propagates
                // errors ( non-failing ranks see INFO(1) = -1 ), but not
                // warnings, so this must run on every rank. The out-of-range
                // index bit ( +1 ) is escalated to a hard error in there
                this->check_warnings();
            }
#else
            BELFEM_ERROR( false, "We are not linked against MUMPS" );
#endif

        }

//------------------------------------------------------------------------------

        string
        MUMPS::error_message(
                const int_t   * aInfo,
                const int_t   & aN,
                const int_t   & aNNZ )
        {
            // the string containing the error message
            string aMessage;

            switch( aInfo[ 0 ] )
            {
                case 0 :
                {
                    aMessage = "No error.";
                    break;
                }
                case -1 :
                {
                    aMessage = sprint( "An error occurred on processor %i.", aInfo[ 1 ] );
                    break;
                }
                case -2 :
                {
                    aMessage = sprint(
                            "Number of nonzeros NNZ=%i is out of range.",
                            aInfo[ 1 ] );
                    break;
                }
                case -3 :
                {
                    aMessage = "Invalid JOB index passed to MUMPS.";
                    break;
                }
                case -4 :
                {
                    aMessage = sprint(
                            "Error in user-provided permutation array PERM_IN at position %i.",
                            aInfo[ 1 ] );
                    break;
                }
                case -5 :
                {
                    aMessage = sprint(
                            "Problem of real workspace allocation of size %s during analysis.",
                            decode_count( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -6 :
                {
                    aMessage = sprint( "Matrix is singular in structure %i.", aInfo[ 1 ] );
                    break;
                }
                case -7 :
                {
                    aMessage = sprint(
                            "Problem of integer workspace allocation of size %s.",
                            decode_count( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -8 :
                {
                    aMessage = "Main internal integer workarray IS too small for factorization.";
                    break;
                }
                case -9 :
                {
                    aMessage = sprint(
                            "Main real/complex workarray S is too small ( missing %s ).",
                            missing_entries( aInfo[ 1 ] ).c_str() );

                    break;
                }
                case -10 :
                {
                    aMessage = "Matrix is singular.";
                    break;
                }
                case -11 :
                {
                    if( aInfo[ 1 ] > 0 )
                    {
                        aMessage = sprint(
                                "Internal real/complex workarray S or LWKUSER is by %i too small for solution.",
                                aInfo[ 1 ] );
                    }
                    else
                    {
                        aMessage =
                                "Internal real/complex workarray S or LWKUSER is too small for solution.";
                    }
                    break;
                }
                case -12 :
                {
                    aMessage = "Internal real/complex workarray S too small for iterative refinement.";
                    break;
                }
                case -13 :
                {
                    aMessage = sprint(
                            "Problem of workspace allocation of size %s.",
                            missing_entries( aInfo[ 1 ] ).c_str() );

                    break;
                }
                case -14 :
                {
                    aMessage = "Internal integer workarray IS too small for solution.";
                    break;
                }
                case -15 :
                {
                    aMessage = "Integer workarray IS too small for iterative refinement and/or error analysis.";
                    break;
                }
                case -16 :
                {
                    aMessage = sprint( "N=%i is out of range.", aInfo[ 1 ] );
                    break;
                }
                case -17 :
                {
                    aMessage = "The internal send buffer that was allocated dynamically by MUMPS on the processor is too small.";
                    break;
                }
                case -18 :
                {
                    aMessage = sprint( "Blocking size for multiple right-hand sides is too large; the suggested maximum is %i.", aInfo[ 1 ] );
                    break;
                }
                case -19 :
                {
                    aMessage = sprint(
                            "Maximum working memory ICNTL(23) is too small for the factorization ( missing %s ).",
                            missing_entries( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -20 :
                {
                    aMessage = sprint( "The internal reception buffer that was allocated dynamically by MUMPS is too small. Need %i.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -21 :
                {
                    aMessage = "Value of PAR=0 is not allowed because only one processor is available";
                    break;
                }
                case -22 :
                {
                    string tLabel;

                    switch( aInfo[ 1 ] )
                    {
                        case  1 : { tLabel = "IRN or ELTPTR";        break; }
                        case  2 : { tLabel = "JCN or ELTVAR";        break; }
                        case  3 : { tLabel = "PERM_IN";              break; }
                        case  4 : { tLabel = "A or AELT";            break; }
                        case  5 : { tLabel = "ROWSCA";               break; }
                        case  6 : { tLabel = "COLSCA";               break; }
                        case  7 : { tLabel = "RHS";                  break; }
                        case  8 : { tLabel = "LISTVAR_SCHUR";        break; }
                        case  9 : { tLabel = "SCHUR";                break; }
                        case 10 : { tLabel = "RHS_SPARSE";           break; }
                        case 11 : { tLabel = "IRHS_SPARSE";          break; }
                        case 12 : { tLabel = "IRHS_PTR";             break; }
                        case 13 : { tLabel = "ISOL_loc";             break; }
                        case 14 : { tLabel = "SOL_loc";              break; }
                        case 15 : { tLabel = "REDRHS";               break; }
                        case 16 : { tLabel = "IRN_loc, JCN_loc or A_loc"; break; }
                        case 17 : { tLabel = "IRHS_loc";             break; }
                        case 18 : { tLabel = "RHS_loc";              break; }
                        default : { tLabel = "unknown";             break; }
                    }

                    aMessage = sprint( "Pointer array %s is not associated, has insufficient size, or was associated but should not be.",
                                       tLabel.c_str() );

                    break;
                }
                case -23 :
                {
                    aMessage = "MPI was not initialized.";
                    break;
                }
                case -24 :
                {
                    aMessage = sprint( "NELT=%i is out of range.", aInfo[ 1 ] );
                    break;
                }
                case -25 :
                {
                    aMessage = "A problem has occurred in the initialization of the BLACS";
                    break;
                }
                case -26 :
                {
                    aMessage = sprint( "LRHS=%i is out of range.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -27 :
                {
                    aMessage = sprint( "NZ_RHS and IRHS_PTR(NRHS+1)=%i do not match.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -28 :
                {
                    aMessage = sprint( "IRHS_PTR(1)=%i is not equal to 1.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -29 :
                {
                    aMessage = sprint( "LSOL_loc=%i is smaller than %i.",
                                       aInfo[ 1 ],
                                       aInfo[ 22 ] );
                    break;
                }
                case -30 :
                {
                    aMessage = sprint( "SCHUR_LLD=%i is out of range.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -31 :
                {
                    if ( aInfo[ 1 ] > 0 )
                    {
                        aMessage = sprint( "MBLOCK=NBLOCK not fulfilled ( MBLOCK=NBLOCK+%i ).", aInfo[ 1 ] );
                    }
                    else
                    {
                        aMessage = sprint( "MBLOCK=NBLOCK not fulfilled ( MBLOCK=NBLOCK-%i ).", -aInfo[ 1 ] );
                    }
                    break;
                }
                case -32 :
                {
                    aMessage = sprint( "Value NRHS=%i not compatible with user defined setting.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -33 :
                {
                    aMessage = "ICNTL(26) was asked for during solve phase (or during the factorization – see ICNTL(32) ) but the Schur complement was not asked for at the analysis phase (ICNTL(19)).";
                    break;
                }
                case -34 :
                {
                    aMessage = sprint( "LREDRHS=%i is out of range.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -35 :
                {
                    aMessage = "The Schur expansion and reduction phases were called in an invalid order.";
                    break;
                }
                case -36 :
                {
                    aMessage = sprint( "Incompatible values of ICNTL(25)=%i and INFOG(28).",
                                       aInfo[ 1 ] );
                    break;
                }
                case -37 :
                {
                    aMessage = sprint( "Value of ICNTL(25) incompatible with some other parameter. ( INFO(2)=%i)",
                                       aInfo[ 1 ] );
                    break;
                }
                case -38 :
                {
                    aMessage = "Parallel analysis was set but neither PT-SCOTCH or ParMetis were provided.";
                    break;
                }
                case -39 :
                {
                    aMessage = "Incompatible values for ICNTL(28) and ICNTL(5) and/or ICNTL(19) and/or ICNTL(6).";
                    break;
                }
                case -40 :
                {
                    aMessage = "The matrix was indicated to be positive definite but is not";
                    break;
                }
                case -41 :
                {
                    aMessage = "Incompatible value of LWKUSER between factorization and solution phases.";
                    break;
                }
                case -42 :
                {
                    aMessage = sprint( "Forward during factorization was set but the value of NRHS=%i on the host is incorrect",
                                       aInfo[ 1 ] );
                    break;
                }
                case -43 :
                {
                    aMessage = sprint( "Incompatible values of ICNTL(32) and ICNTL(%i).",
                                       aInfo[ 1 ] );
                    break;
                }
                case -44 :
                {
                    aMessage = sprint( "The solve phase (JOB=3) cannot be performed because the factors or part of the factors are not\n"
                                       "available ( ICNTL(31)=%i)", aInfo[ 1 ] );
                    break;
                }
                case -45 :
                {
                    aMessage = sprint( "NRHS=%i must be >0.", aInfo[ 1 ] );
                    break;
                }
                case -46 :
                {
                    aMessage = sprint( "NZ_RHS=%i must be >0.", aInfo[ 1 ] );
                    break;
                }
                case -47 :
                {
                    aMessage = sprint( "Entries of A^(-1) require NRHS=N ( NRHS=%i ).",
                                       aInfo[ 1 ] );
                    break;
                }
                case -48 :
                {
                    aMessage = sprint( "Requested entries of A^(-1) are incompatible with ICNTL(30) and ICNTL(%i).",
                                       aInfo[ 1 ] );
                    break;
                }
                case -49 :
                {
                    aMessage = sprint( "Incorrect value: SIZE_SCHUR=%i", aInfo[ 1 ] );
                    break;
                }
                case -50 :
                {
                    aMessage = "An error occurred while computing the fill-reducing ordering during the analysis phase.";
                    break;
                }
                case -51 :
                {
                    aMessage = sprint( "The linked ordering library cannot handle graphs of size %s ( must be < 2^31-1).",
                                       decode_count( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -52 :
                {
                    aMessage = "The linked ordering library must have 64-bit default integers";
                    break;
                }
                case -53 :
                {
                    aMessage = "Internal error that could be due to inconsistent input data between two consecutive calls.";
                    break;
                }
                case -54 :
                {
                    aMessage = "BLR compression was requested for the factorization but not for the analysis. Rerun the analysis phase with BLR compression enabled.";
                    break;
                }
                case -55 :
                {
                    if( aInfo[ 1 ] > 0 )
                    {
                        aMessage = sprint( "Distributed right-hand side has a leading dimension LRHS_loc=%i that is too small.",
                                           aInfo[ 1 ] );
                    }
                    else
                    {
                        aMessage = sprint( "Distributed right-hand side has Nloc_RHS=%i nonzero local rows on a non-working host.",
                                           -aInfo[ 1 ] );
                    }
                    break;
                }
                case -56 :
                {
                    aMessage = sprint( "Distributed right-hand side and solution share storage, but LRHS_loc=%i is smaller than LSOL_loc.",
                                       aInfo[ 1 ] );
                    break;
                }
                case -57 :
                {
                    aMessage = sprint( "Error in the user-provided block-format interface ( detail %i ).", aInfo[ 1 ] );
                    break;
                }
                case -58 :
                {
                    aMessage = sprint( "Error in the OpenMP thread setup ( detail %i ).", aInfo[ 1 ] );
                    break;
                }
                case -60 :
                {
                    aMessage = sprint( "Error in the MPI-to-OpenMP feature ( detail %i ).", aInfo[ 1 ] );
                    break;
                }
                case -61 :
                {
                    aMessage = "Error in the MPI-to-OpenMP feature during an affinity-mask operation.";
                    break;
                }
                case -62 :
                {
                    aMessage = "Error in the MPI-to-OpenMP feature during memory allocation.";
                    break;
                }
                case -63 :
                {
                    aMessage = sprint( "The MUMPS configuration file could not be opened ( detail %i ).", aInfo[ 1 ] );
                    break;
                }
                case -64 :
                {
                    aMessage = sprint( "Read error in the MUMPS configuration file at line %i.", aInfo[ 1 ] );
                    break;
                }
                case -65 :
                {
                    aMessage = sprint( "Syntax error in the MUMPS configuration file at line %i.", aInfo[ 1 ] );
                    break;
                }
                case -66 :
                {
                    aMessage = sprint( "A file referenced in the MUMPS configuration file at line %i could not be opened or read.", aInfo[ 1 ] );
                    break;
                }
                case -67 :
                {
                    aMessage = sprint( "Memory allocation failed while reading the MUMPS configuration file at line %i.", aInfo[ 1 ] );
                    break;
                }
                case -69 :
                {
                    aMessage = sprint( "The Fortran INTEGER size does not match MUMPS_INT ( MUMPS_INT is %i bytes ). Rebuild MUMPS with a consistent integer size.", aInfo[ 1 ] );
                    break;
                }
                case -70 :
                {
                    aMessage = "The instance-save file already exists. Delete it or choose a different save path before saving.";
                    break;
                }
                case -71 :
                {
                    aMessage = "Could not create a file needed to save the MUMPS instance.";
                    break;
                }
                case -72 :
                {
                    aMessage = sprint( "Write error while saving the MUMPS instance ( write size %s ).", decode_count( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -73 :
                {
                    string tLabel;

                    switch( aInfo[ 1 ] )
                    {
                        case 1 : { tLabel = "Fortran version";                       break; }
                        case 2 : { tLabel = "integer size";                          break; }
                        case 3 : { tLabel = "MPI compatibility of the saved instance"; break; }
                        case 4 : { tLabel = "number of MPI processes";               break; }
                        case 5 : { tLabel = "arithmetic";                            break; }
                        case 6 : { tLabel = "symmetry mode SYM";                     break; }
                        case 7 : { tLabel = "host-working mode PAR";                 break; }
                        default : { tLabel = "unknown";                             break; }
                    }

                    aMessage = sprint( "Saved instance is incompatible with the current one ( parameter: %s ).",
                                       tLabel.c_str() );
                    break;
                }
                case -74 :
                {
                    aMessage = sprint( "Could not open the saved-instance file for restore on MPI process %i.", aInfo[ 1 ] );
                    break;
                }
                case -75 :
                {
                    aMessage = sprint( "Read error while restoring the MUMPS instance ( remaining read size %s ).", decode_count( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -76 :
                {
                    aMessage = sprint( "Could not delete the saved-instance files on MPI process %i.", aInfo[ 1 ] );
                    break;
                }
                case -77 :
                {
                    aMessage = sprint( "The save directory or prefix is invalid ( detail %i ).", aInfo[ 1 ] );
                    break;
                }
                case -78 :
                {
                    aMessage = sprint( "Problem of workspace allocation of size %s during the restore step.", decode_count( aInfo[ 1 ] ).c_str() );
                    break;
                }
                case -79 :
                {
                    aMessage = sprint( "MUMPS could not find a free Fortran file unit ( context %i ).", aInfo[ 1 ] );
                    break;
                }
                case -80 :
                {
                    aMessage = sprint( "The factorization phase was cancelled on user request on MPI process %i.", aInfo[ 1 ] );
                    break;
                }
                case -81 :
                {
                    aMessage = sprint( "The solution phase was cancelled on user request on MPI process %i.", aInfo[ 1 ] );
                    break;
                }
                case -88 :
                {
                    aMessage = sprint( "SCOTCH ordering failed ( SCOTCH error code %i ).", aInfo[ 1 ] );
                    break;
                }
                case -89 :
                {
                    aMessage = sprint( "SCOTCH k-way partitioning failed ( SCOTCH error code %i ). Consider making METIS available to MUMPS.", aInfo[ 1 ] );
                    break;
                }
                case -90 :
                {
                    aMessage = "Error in the out-of-core management. See the MUMPS error output for more information.";
                    break;
                }
                case -800 :
                {
                    aMessage = sprint( "Temporary error specific to this MUMPS release ( detail %i ).", aInfo[ 1 ] );
                    break;
                }
                default :
                {
                    aMessage = sprint( "Unknown MUMPS error ( INFO(1)=%i ).", aInfo[ 0 ] );
                    break;
                }

            }

            return aMessage;
        }
//------------------------------------------------------------------------------

        void
        MUMPS::warning_message(
                const int_t    * aInfo,
                Cell< string > & aWarnings )
        {
            // MUMPS encodes warnings as a sum of bit flags in INFO(1) ( e.g.
            // INFO(1)=6 combines warnings +2 and +4 ), so decode each bit
            // independently and append one line per active warning.
            const int_t tFlags = aInfo[ 0 ];

            // Only warnings +1 and +16 carry a detail in INFO(2), and INFO(2) then
            // belongs to whichever of them was raised last. If both are set we cannot
            // attribute it, so the count is suppressed for those two lines.
            const bool tDetailIsAmbiguous = ( tFlags & 1 ) && ( tFlags & 16 );

            if( tFlags & 1 )
            {
                if( tDetailIsAmbiguous )
                {
                    aWarnings.push(
                            "Indices in IRN or JCN were out of range and were ignored." );
                }
                else
                {
                    aWarnings.push( sprint(
                            "Indices in IRN or JCN were out of range and ignored ( %i faulty entries ).",
                            aInfo[ 1 ] ) );
                }
            }
            if( tFlags & 2 )
            {
                aWarnings.push(
                        "During error analysis the max-norm of the computed solution is close to zero." );
            }
            if( tFlags & 4 )
            {
                aWarnings.push(
                        "Not enough memory to compact the internal workarray at the end of the factorization." );
            }
            if( tFlags & 8 )
            {
                aWarnings.push(
                        "Iterative refinement did not converge within the allowed number of steps." );
            }
            if( tFlags & 16 )
            {
                if( tDetailIsAmbiguous )
                {
                    aWarnings.push(
                            "Rank-revealing: inertia and/or determinant may be inconsistent with the detected singularities." );
                }
                else
                {
                    aWarnings.push( sprint(
                            "Rank-revealing: inertia and/or determinant may be inconsistent with the detected singularities ( deficiency %i ).",
                            aInfo[ 1 ] ) );
                }
            }

            if( aWarnings.size() == 0 )
            {
                aWarnings.push( sprint( "Unknown MUMPS warning ( INFO(1)=%i ).", tFlags ) );
            }
        }

//------------------------------------------------------------------------------

        void
        MUMPS::print_soft_fail()
        {
            if ( this->rank() != 0 ) return ;

            // compact reason for the log box ( inner width 71, matching the
            // controller's timestep box ). Decodes the rank-uniform INFOG:
            // with the rank-local INFO, rank 0 only ever saw the propagated
            // -1 and the -9/-8 arms below were dead on a multi-rank job
            // ( observed failing ranks 7, 4, 1, 6 )
            // sized for the worst-case decimal widths; the print clips to
            // the 53-char field of the box row
            char tReason[ 96 ];
            const int tCode = ( int ) mInfoG( 0 );

            // the workspace arms: reaching them means escalate_workspace()
            // gave up, so the relaxation and the cap printed here are what
            // was in force. The numbers behind the decision go on a second
            // row below, which the 53-char reason field cannot hold
            const bool tWorkspace = tCode == -9 || tCode == -8 || tCode == -19
                                 || tCode == -17 || tCode == -20 ;
            const int_t tRelax = mIParameters(
                static_cast< index_t >( mumps::Parameter::MemoryRelaxation ) );
            const int_t tCap = mIParameters(
                static_cast< index_t >( mumps::Parameter::MemoryBudget ) );

            if ( tCode == -9 || tCode == -8 )
            {
                // three ways to arrive here, and the reason must name the
                // right one ( audit finding: the old text said "exhausted"
                // when the ladder had never been tried )
                const char * tWhat = tCode == -9 ?
                    "out of workspace" : "out of integer workspace" ;

                if ( tCap == 0 && mMeasuredBudgetMB > 0 )
                {
                    // the machine's budget was below MUMPS's estimate: no
                    // cap was set and no rung was climbed
                    std::snprintf( tReason, sizeof( tReason ),
                        "%s, budget %li MB below the estimate",
                        tWhat, ( long ) mMeasuredBudgetMB );
                }
                else if ( tRelax >= gMaxMemoryRelaxation )
                {
                    std::snprintf( tReason, sizeof( tReason ),
                        "%s ( ICNTL(14) = %i exhausted )",
                        tWhat, ( int ) tRelax );
                }
                else
                {
                    std::snprintf( tReason, sizeof( tReason ),
                        "%s ( ICNTL(14) = %i )",
                        tWhat, ( int ) tRelax );
                }
            }
            else if ( tCode == -19 )
            {
                std::snprintf( tReason, sizeof( tReason ),
                    "memory cap ICNTL(23) = %li MB cannot be met",
                    ( long ) tCap );
            }
            else if ( tCode == -17 || tCode == -20 )
            {
                std::snprintf( tReason, sizeof( tReason ),
                    "MPI %s buffer too small ( ICNTL(14) = %i%s )",
                    tCode == -17 ? "send" : "reception", ( int ) tRelax,
                    tRelax >= gMaxMemoryRelaxation ? " exhausted" : "" );
            }
            else if ( tCode == -1 )
            {
                std::snprintf( tReason, sizeof( tReason ),
                    "error on rank %i", ( int ) mInfoG( 1 ) );
            }
            else if ( tCode == -1000 )
            {
                // mumpstools shim code, not a MUMPS INFO value: the instance
                // registry had no free slot for the create
                std::snprintf( tReason, sizeof( tReason ),
                    "solver registry exhausted, no instance created" );
            }
            else
            {
                std::snprintf( tReason, sizeof( tReason ),
                    "error code %i", tCode );
            }

            message( InfoLevel::Minimal,
                "   │ MUMPS soft fail: %-53.53s│", tReason );

            if ( tWorkspace )
            {
                // estimate: INFOG(16), or INFOG(36) under BLR ( 0-based
                // here ); cap: 0 means MUMPS sized from the estimate alone
                const int_t tBlr = mIParameters(
                    static_cast< index_t >( mumps::Parameter::CompressionMode ) );
                const int_t tEstimate = ( tBlr == 1 || tBlr == 2 ) ?
                    mInfoG( 35 ) : mInfoG( 15 );

                char tRow[ 128 ];
                std::snprintf( tRow, sizeof( tRow ),
                    " estimate %li MB, cap %li MB, short by %s",
                    ( long ) tEstimate, ( long ) tCap,
                    missing_entries( mInfoG( 1 ) ).c_str() );

                message( InfoLevel::Minimal, "   │%-71.71s│", tRow );
            }

            // what happens next depends on WHO soft-failed, and this box used
            // to promise a timestep cut unconditionally. That is true for the
            // timestepping solve, whose failure the controller answers by
            // halving delta t. It is false for a failed CREATE reached from
            // the conditioning diagnostic, which arms soft-fail itself, records
            // SolverFailed and lets the step continue untouched -- telling that
            // user their timestep is about to be cut sends them after the wrong
            // subsystem, the same class of defect as the Krylov-budget mislabel
            // this round exists to fix
            message( InfoLevel::Minimal,
                "   │%-71s│",
                tCode == -1000
                    ? " no instance was created - the caller decides how to proceed"
                    : " handing the step back to the controller ( timestep will be cut )" );
        }

//------------------------------------------------------------------------------

        void
        MUMPS::check_warnings()
        {
            // INFO(1) = +1 means MUMPS found indices in IRN or JCN outside the
            // matrix, DROPPED those entries and factorized anyway. The solve then
            // succeeds on a different matrix than the one BELFEM assembled, and
            // nothing downstream can notice: the residual is measured against the
            // assembled operator, so no tolerance ever sees the missing entries.
            //
            // BELFEM never hands MUMPS an out-of-range index on purpose. Fixed
            // dofs are not discarded through a sentinel index, they live in their
            // own matrix ( mDirichletMatrix is free x fixed, mSystemMatrix is
            // free x free ), so there is no discard idiom for this to break.
            // That makes +1 an indexing defect in the assembly, not a numerical
            // condition, and a warning is the wrong tier for it.
            //
            // Test the BIT, never mInfo( 0 ) > 0: MUMPS returns the warnings as a
            // sum of flags, and +8 ( iterative refinement did not converge ) is
            // expected on ill-conditioned HTS systems. Aborting on any positive
            // INFO(1) would kill runs that are behaving normally.
            //
            // BELFEM_ERROR, not BELFEM_ASSERT: the latter compiles out under
            // NDEBUG, which would leave the release build quietly solving the
            // wrong matrix - exactly what the always-active tier exists for.
            //
            // MPI: the matrix is centralized on the host ( ICNTL(18) = 0, only
            // the master binds irn/jcn/A ), so +1 can only be raised where the
            // indices are. Aborting from that one rank is intended; unlike the
            // soft-fail path, this one does not return into another collective
            if( mInfo( 0 ) & 1 )
            {
                // INFO(2) carries the faulty-entry count, but +16 writes INFO(2)
                // as well and the value then belongs to whichever of the two was
                // raised last. If both bits are set the count cannot be
                // attributed, so it is suppressed ( same rule as warning_message )
                const string tDetail = ( mInfo( 0 ) & 16 ) ?
                        string( "" ) :
                        sprint( " ( %i faulty entries )", ( int ) mInfo( 1 ) );

                BELFEM_ERROR( false,
                        "MUMPS ignored out-of-range IRN or JCN entries%s on proc %i ( INFO(1) = %i ):\n"
                        "the factorization is not of the matrix that was assembled.\n"
                        "This is an indexing defect in the assembly, not a numerical condition.",
                        tDetail.c_str(),
                        ( int ) this->rank(),
                        ( int ) mInfo( 0 ) );
            }

            // the surviving warnings are RANK-LOCAL: unlike errors, MUMPS does
            // not propagate them, so reading INFO on rank 0 alone loses every
            // warning raised only on a worker. Each rank reports its own, tagged
            // with the proc it happened on
            Cell< string > tWarnings;
            this->warning_message( mInfo.data(), tWarnings );

            for( const string & tWarning : tWarnings )
            {
                message( InfoLevel::Minimal,
                         "MUMPS returned a warning on proc %i ( INFO(1) = %i ): %s",
                         ( int ) this->rank(),
                         ( int ) mInfo( 0 ),
                         tWarning.c_str() );
            }
        }

//------------------------------------------------------------------------------

        void
        MUMPS::set_reordering(
                const MumpsSerialReodrdering   aSerial,
                const MumpsParallelReodrdering aParallel
        )
        {
            mIParameters( static_cast< index_t > ( mumps::Parameter::SerialReordering ) )
                = static_cast< int_t >( aSerial );
            mIParameters( static_cast< index_t > ( mumps::Parameter::ParallelReordering ) )
                = static_cast< int_t >( aParallel );
        }

//------------------------------------------------------------------------------

        void
        MUMPS::set_block_low_ranking(
                const MumpsBlockLowRanking aBLK,
                const real            aEpsilon
        )
        {
            mIParameters( static_cast< index_t >( mumps::Parameter::CompressionMode) )
                = static_cast< int_t > ( aBLK );
            mRParameters( 0 ) = aEpsilon ;
        }
//------------------------------------------------------------------------------

        void
        MUMPS::set_error_analysis(
                const MumpsErrorAnalysis aSetting )
        {
            mIParameters( static_cast< index_t >( mumps::Parameter::ErrorAnalysis ) )
                = static_cast< int_t >( aSetting );
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_determinant() const
        {
            return mIParameters( static_cast< index_t >( mumps::Parameter::ComputeDeterminant )) == 0 ?
                   BELFEM_SIGNALING_NAN :
                   mRInfoG( 11 ) ;
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_cond1() const
        {
            return mIParameters( static_cast< index_t >( mumps::Parameter::ErrorAnalysis ) ) == 1 ?
                   mRInfoG( 9 )  :
                   BELFEM_SIGNALING_NAN ;
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_cond2() const
        {
            return mIParameters( static_cast< index_t >( mumps::Parameter::ErrorAnalysis ) ) == 1 ?
                   mRInfoG( 10 )  :
                   BELFEM_SIGNALING_NAN ;
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_forward_error() const
        {
            // RINFOG( 9 ), zero based. Only filled by the FULL error analysis
            // ( ICNTL( 11 ) = 1 ) -- the cheaper setting stops at the backward
            // errors and leaves this one untouched
            return mIParameters( static_cast< index_t >( mumps::Parameter::ErrorAnalysis ) ) == 1 ?
                   mRInfoG( 8 )  :
                   BELFEM_SIGNALING_NAN ;
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_backward_error() const
        {
            // RINFOG( 7 ) and RINFOG( 8 ), zero based: omega1 and omega2.
            // These are not two rival measures -- they are ONE componentwise
            // backward error split by row. Arioli, Demmel & Duff measure row i
            // as | r_i | / ( |A||x| + |b| )_i, but that denominator becomes
            // meaningless where the row nearly cancels, so rows below a
            // threshold move to omega2 with a denominator padded by
            // ||A_i||_inf * ||x||_inf instead. Rows above it stay in omega1.
            //
            // They therefore ADD rather than compete, which is also how MUMPS
            // builds the forward error it reports:
            //
            //     RINFOG( 9 )  =  omega1 * COND1  +  omega2 * COND2
            if( mIParameters( static_cast< index_t >( mumps::Parameter::ErrorAnalysis ) ) != 1 )
            {
                return BELFEM_SIGNALING_NAN ;
            }

            return mRInfoG( 6 ) + mRInfoG( 7 ) ;
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_omega2() const
        {
            // RINFOG( 8 ), zero based: omega2 on its own. It exists as its
            // own accessor because it is what says whether the COND2 beside
            // it CONTRIBUTES anything.
            //
            // MUMPS builds the forward error as
            //
            //     ERX  =  omega1 * COND1  +  omega2 * COND2
            //
            // so omega2 == 0 means the second term is empty, whatever COND2
            // holds. That is the statement this accessor supports, and the
            // only one it supports.
            //
            // It does NOT tell you whether MUMPS estimated COND2 or skipped
            // it. Those are close but not the same, and an earlier version of
            // this comment claimed the stronger thing. The membership test is
            // IW( i, 1 ), which is internal: DMUMPS_SOL_OMEGA assigns
            // IW( i, 1 ) = 2 OUTSIDE the IF ( TAU > 0 ) guard that updates
            // omega2 ( dsol_aux.F:891-900 ), so a row with TAU == 0 -- or with
            // an exactly zero residual -- joins the second category, sets
            // LCOND2, and lets the Hager estimator run, while leaving omega2
            // at its initial zero. Empty category implies omega2 == 0; the
            // converse does not hold, and MUMPS documents only the forward
            // direction ( dini_defaults.F:1062-1065 ).
            //
            // Still test THIS rather than COND2 == 1.0. The skipped path
            // leaves COND2 at the 1.0 both condition numbers are initialized
            // to ( dsol_aux.F:959-964 ), which is indistinguishable from a
            // genuine estimate landing near one, while omega2 is an exact
            // quantity MUMPS assigns as a literal
            return mIParameters( static_cast< index_t >( mumps::Parameter::ErrorAnalysis ) ) == 1 ?
                   mRInfoG( 7 )  :
                   BELFEM_SIGNALING_NAN ;
        }

//------------------------------------------------------------------------------
    }
}

//------------------------------------------------------------------------------

namespace mumps
{
    WorkspaceAction
    next_workspace_action(
            const belfem::int_t aInfoG1,
            const belfem::int_t aRelax,
            const belfem::int_t aRelaxCeiling,
            const belfem::int_t aBudgetSlot,
            const belfem::int_t aBudgetMB,
            const belfem::int_t aBoundMB )
    {
        // only the factorization's workspace codes are retryable here; -19
        // is the cap itself refusing, and no larger ICNTL(14) fixes that
        if ( aInfoG1 != -9 && aInfoG1 != -8 && aInfoG1 != -17 && aInfoG1 != -20 )
        {
            return WorkspaceAction::GiveUp ;
        }

        // the MPI buffers are sized from ICNTL(14) alone -- first, before
        // the workarrays, once a cap is in place ( guide, ICNTL(23) ) -- so
        // a cap cannot widen them and the rung is the only remedy
        if ( aInfoG1 == -17 || aInfoG1 == -20 )
        {
            return aRelax < aRelaxCeiling ?
                WorkspaceAction::Ladder : WorkspaceAction::GiveUp ;
        }

        // the cap comes first, once: with the slot empty and the machine
        // measured, hand MUMPS everything it may have -- unless MUMPS's own
        // estimate already exceeds it, in which case a ladder rung would
        // only ask for more than the machine holds
        if ( aBudgetSlot == 0 && aBudgetMB > 0 )
        {
            return aBudgetMB >= aBoundMB ?
                WorkspaceAction::Cap : WorkspaceAction::GiveUp ;
        }

        // the ladder behind the cap ( or instead of it when the budget is
        // unknown ): the guide keeps asking for a larger ICNTL(14) on a
        // residual -9 / -8 even with ICNTL(23) set
        return aRelax < aRelaxCeiling ?
            WorkspaceAction::Ladder : WorkspaceAction::GiveUp ;
    }
}
