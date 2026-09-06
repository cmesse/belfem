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
#include "commtools.hpp"
#include "arpacktools.hpp"
#ifdef BELFEM_PARPACK
#include "parpacktools.hpp"
#include "cl_SolverDistMatrix.hpp"
#endif
#include "cl_FEM_DofMgr_EigenValues.hpp"

#include <random.hpp>

#include "cl_Logger.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "cl_Timer.hpp"
#include "cl_IWG_Timestep.hpp"

// for the snprintf in report_unavailable
#include <cstdio>

namespace belfem
{

    namespace fem
    {
        namespace dofmgr
        {
//------------------------------------------------------------------------------

            EigenValues::EigenValues( DofManager * aParent ) :
                 mCommRank( comm_rank() ),
                 mCommSize( comm_size() ),
                 mParent( aParent )
            {
                mInfo.set_size( arpack::gNumInfoEntries, 0 );
            }

            EigenValues::~EigenValues()
            {
                if( mSolver != nullptr )
                {
                    delete mSolver ;
                }
                if( mShiftInvertSolver != nullptr )
                {
                    delete mShiftInvertSolver ;
                }
                if( mM != nullptr )
                {
                    delete mM ;
                }
#ifdef BELFEM_PARPACK
                if( mDistMatrix != nullptr )
                {
                    delete mDistMatrix ;
                }
#endif
            }

//------------------------------------------------------------------------------

            void
            EigenValues::link_matrix()
            {
                mK = mParent->system_matrix() ;
                BELFEM_ERROR( mK->n_rows() == mK->n_cols(),
                    "Matrix must be quadratic to compute eigenvalues. ( is %lu x %lu )",
                    ( long unsigned int ) mK->n_rows(),
                    ( long unsigned int ) mK->n_cols() );

                mN   = mK->n_rows() ;
                mNNZ = mK->number_of_nonzeros() ;

                // mM mirrors mK's pattern, so it is built BEFORE the indexing
                // base is switched -- it copies the pointer array once and
                // compute_matrices() only ever refreshes its values, so a
                // Fortran based pattern here would stick for its lifetime.
                // The rebuild tests nnz as well as the dimension: a same-size
                // pattern with more entries would otherwise be memcpy'd into
                // the old, shorter allocation
                if ( mCommRank == 0 )
                {
                    if( mM == nullptr
                        || ( int_t ) mM->n_cols() != mN
                        || ( int_t ) mM->number_of_nonzeros() != mNNZ )
                    {
                        if( mM != nullptr )
                        {
                            delete mM ;
                        }

                        mM = new SpMatrix(
                            mK->type(),
                            mN, mN,
                            mNNZ,
                            mK->type() == SpMatrixType::CSC ? mK->rows() : mK->cols(),
                            mK->pointers() );
                    }
                }

                // ARPACK is a Fortran library and reads mK's arrays directly,
                // so the base is switched here and put back by
                // restore_indexing_base() at the single exit of run().
                // run_parpack does not depend on this: DistMatrix::create_matrix
                // forces Cpp before the pattern is copied
                mOriginalBase = mK->indexing_base() ;

                if ( mOriginalBase == 0 )
                {
                    mK->set_indexing_base( SpMatrixIndexingBase::Fortran );
                }
                this->set_num_minvals( mNumMinVals );
                this->set_num_maxvals( mNumMaxVals );

            }

//------------------------------------------------------------------------------

            void
            EigenValues::set_num_minvals( const index_t aNumMinVals )
            {
                BELFEM_ERROR( aNumMinVals > 0,
                    "Number of eigenvalues must be positive." );

                mNumMinVals = aNumMinVals ;

                index_t tNumVals = std::max( mNumMinVals, mNumMaxVals ) + 1 ;
                mLambdaReal.set_size( tNumVals );
                mLambdaImag.set_size( tNumVals );

            }

//------------------------------------------------------------------------------

            void
            EigenValues::set_num_maxvals( const index_t aNumMaxVals )
            {
                BELFEM_ERROR( aNumMaxVals > 0,
                    "Number of eigenvalues must be positive." );

                mNumMaxVals = aNumMaxVals ;

                index_t tNumVals = std::max( mNumMinVals, mNumMaxVals ) + 1 ;
                mLambdaReal.set_size( tNumVals );
                mLambdaImag.set_size( tNumVals );
            }

//------------------------------------------------------------------------------

            void
            EigenValues::set_tolerance( const real aTolerance )
            {
                // ARPACK reads tol <= 0 as "use machine precision", which is
                // the one setting guaranteed to be unreachable at the small
                // end of an ill conditioned spectrum
                BELFEM_ERROR( aTolerance > 0.0,
                    "Eigenvalue tolerance must be positive ( is %e ).",
                    ( double ) aTolerance );

                mEpsilon = aTolerance ;
                mToleranceExplicit = true ;
            }

//------------------------------------------------------------------------------

            void
            EigenValues::set_subspace_size( const index_t aSubspaceSize )
            {
                // ARPACK needs ncv - nev >= 2, and nev can still grow through
                // the setters, so only the hard floor is enforced here --
                // ( p )dnaupd reports -3 if the final combination is invalid
                BELFEM_ERROR( aSubspaceSize > 2,
                    "Krylov subspace size must be larger than 2 ( is %lu ).",
                    ( long unsigned int ) aSubspaceSize );

                mSubspaceSize = aSubspaceSize ;
                mSubspaceSizeExplicit = true ;
            }

//------------------------------------------------------------------------------

            void
            EigenValues::set_symmetric( const bool aSymmetric )
            {
                mSymmetric = aSymmetric ;
            }

//------------------------------------------------------------------------------

            void
            EigenValues::set_max_iterations( const index_t aNumMaxIter )
            {
                BELFEM_ERROR( aNumMaxIter > 0,
                    "Number of restart iterations must be positive." );

                mNumMaxIter = aNumMaxIter ;
                mNumMaxIterExplicit = true ;
            }

//------------------------------------------------------------------------------

            bool
            EigenValues::configure(
                    const int_t aJob,
                    const bool  aFolded,
                    const int_t aNumEigenValues,
                    int_t     & aSubspaceSize,
                    int_t     & aNumMaxIter,
                    real      & aTolerance )
            {
                // Two different floors, and the one that matters is the larger.
                // ARPACK's hard requirement is ncv - nev >= 2, hence
                // ncv >= nev + 2 ; remark 4's ncv >= 2*nev + 1 is only a
                // recommendation. But the DRIVERS apply that recommendation
                // themselves -- both compute ncv = max( 2*nev + 1, ncvmin )
                // before allocating -- so a subspace sized against nev + 2
                // alone would be checked against the budget and then allocated
                // larger. The feasibility bound is what will actually be
                // allocated
                const int_t tSubspaceMin = std::max( aNumEigenValues + 2,
                                                     2 * aNumEigenValues + 1 );

                // what we would like: remark 4 is explicit that raising ncv at
                // fixed nev usually REDUCES the total number of OP*x products,
                // even though each restart costs more
                const int_t tSubspaceWant =
                        std::max( 2 * aNumEigenValues + 20, ( int_t ) 20 );

                // the Arnoldi basis is n x ncv doubles and is the largest
                // object in the algorithm. Both products are formed in a
                // 64-bit type first: int_t may be 32-bit, where 8*n already
                // overflows around n = 2.7e8
                const size_t tRowBytes = 8UL * ( size_t ) mN ;

                int_t tSubspaceCap = mN ;

                if( tRowBytes > 0 )
                {
                    const size_t tFits = mBasisBudget / tRowBytes ;

                    if( ( size_t ) tSubspaceCap > tFits )
                    {
                        tSubspaceCap = ( int_t ) tFits ;
                    }
                }

                // no legal subspace fits the budget, so nothing is attempted --
                // reported rather than silently allocated past
                if( tSubspaceMin > tSubspaceCap )
                {
                    return false ;
                }

                if( mSubspaceSizeExplicit )
                {
                    aSubspaceSize = mSubspaceSize ;

                    // authoritative, but not exempt: an explicit value below
                    // the floor would be raised by the driver anyway ( silently,
                    // so it was never really authoritative ), and one above the
                    // cap would allocate past the basis budget. Clamped here,
                    // where the caller can be told
                    if( aSubspaceSize < tSubspaceMin || aSubspaceSize > tSubspaceCap )
                    {
                        const int_t tClamped =
                                aSubspaceSize < tSubspaceMin ? tSubspaceMin
                                                             : tSubspaceCap ;

                        message( InfoLevel::Verbose,
                            "Krylov subspace %i is not usable at nev = %i, n = %i "
                            "( legal range %i to %i ) - using %i.\n",
                            ( int ) aSubspaceSize,
                            ( int ) aNumEigenValues,
                            ( int ) mN,
                            ( int ) tSubspaceMin,
                            ( int ) tSubspaceCap,
                            ( int ) tClamped );

                        aSubspaceSize = tClamped ;
                    }
                }
                else
                {
                    aSubspaceSize = tSubspaceWant < tSubspaceMin ? tSubspaceMin :
                                    tSubspaceWant > tSubspaceCap ? tSubspaceCap :
                                                                   tSubspaceWant ;
                }

                // maxit is a RESTART cap, NOT a product ceiling: a restart does
                // not cost ncv - nev products ( ARPACK boosts nev internally
                // between restarts, so np shrinks -- measured 10.1 and 12.2
                // products per restart at ncv - nev = 19 ). The achieved count
                // is reported from info( gInfoNumOperations ) rather than
                // modeled.
                //
                // The conditioning diagnostic asks only for exterior values
                // now, and those converge in a handful of restarts. Job 0 is
                // the unfolded INTERIOR request, which nothing in the tree
                // reaches any more but which compute_smallest_eigenvalues()
                // still exposes -- it keeps a budget wide enough to leave that
                // public entry point at least as capable as it was at the old
                // fixed maxit = 300
                // maxit is a CAP, not a target: a run that converges in 40
                // restarts pays 40 whatever the cap is, so headroom is only
                // ever paid for by a run that was going to fail anyway -- and
                // the latch bounds how often that can happen. Budget by how
                // hard the request is:
                //   job 0    unfolded INTERIOR 'SM'. Nothing in the tree reaches
                //            it now, but compute_smallest_eigenvalues() still
                //            exposes it and it is the expensive one
                //   folded   exterior, but at 1e-10 rather than 1e-4 -- six more
                //            decades of accuracy, which costs restarts
                //   plain    exterior at 1e-4, the cheap end
                const int_t tProductBudget = aJob == 0 ? 20000 :
                                             aFolded   ? 12000 :
                                                          4000 ;

                if( mNumMaxIterExplicit )
                {
                    aNumMaxIter = mNumMaxIter ;
                }
                else
                {
                    const int_t tDerived = tProductBudget
                            / std::max( ( int_t ) 1, aSubspaceSize - aNumEigenValues );

                    // mNumMaxIter is a FLOOR, not a value to be replaced. The
                    // budget may only ever RAISE the restart cap.
                    //
                    // This cost a cycle on 2026-08-28: the budget was
                    // calibrated on a 2026-08-10 measurement at n = 10428 and
                    // applied unguarded to a matrix nine times larger, which
                    // cut the exterior run from 300 restarts to 95 and broke an
                    // end that had been converging. A heuristic that can weaken
                    // a working configuration is not a heuristic, it is a
                    // regression with a formula in front of it
                    aNumMaxIter = std::max( tDerived, mNumMaxIter );
                }

                // the folded run pays for a cancellation and needs the tight
                // tolerance ; the plain large end does not. An explicit
                // set_tolerance() overrides both -- and because that rule lives
                // here alone, the chosen value is recorded rather than
                // reconstructed by whoever reports it
                aTolerance = mToleranceExplicit ? mEpsilon :
                             aFolded            ? mFoldEpsilon :
                                                  mEpsilon ;

                mEffectiveTolerance = aTolerance ;

                return true ;
            }

//------------------------------------------------------------------------------

            void
            EigenValues::compute_matrices()
            {
                this->link_matrix() ;

                IWG_Timestep * tIwg = reinterpret_cast< IWG_Timestep * > ( mParent->iwg() );

                EulerMethod tMethod = tIwg->method() ;

                tIwg->set_timestepping_method( EulerMethod::MassOnly );
                mParent->compute_jacobian_and_rhs( );

                if( mParent->parent()->is_master() )
                {
                   std::memcpy( mM->data(), mK->data(), mNNZ * sizeof( real ) );
                }

                tIwg->set_timestepping_method( EulerMethod::StiffnessOnly );
                mParent->compute_jacobian_and_rhs();

                // restore timestepping method
                tIwg->set_timestepping_method( tMethod );

            }

//------------------------------------------------------------------------------

            real
            EigenValues::run( const int_t aJob, const real aSigma )
            {
                real aResult = BELFEM_QUIET_NAN ;

                // this run's verdict, not the last one's
                mSubspaceInfeasible = false ;

                // < 2 rather than == 1: a single rank and an uninitialized
                // communicator both take the cheap path, where the matrix is
                // already whole on the calling rank and there is no row map,
                // no scatter and no allgatherv to pay for
                if ( mCommSize < 2 )
                {
                    mBackendLabel = "ARPACK" ;

                    aResult = this->run_arpack( aJob, aSigma );
                }
                else
                {
#ifdef BELFEM_PARPACK
                    // collective: every rank enters and every rank comes back
                    // with the same value
                    mBackendLabel = "PARPACK" ;

                    aResult = this->run_parpack( aJob, aSigma );
#elif defined( BELFEM_ARPACK )
                    // no distributed driver -- the master owns the whole
                    // matrix, so it computes alone and hands the answer round
                    if( mParent->parent()->is_master() )
                    {
                        aResult = this->run_arpack( aJob, aSigma );
                    }
                    comm_barrier();

                    broadcast( aResult );
#else
                    BELFEM_ERROR( false,
                        "Can not compute eigenvalues: BELFEM was built without ARPACK-ng.\n"
                        "                 Reconfigure with USE_ARPACK=ON.\n" );
#endif
                }

                // single exit, deliberately: link_matrix() re-captures
                // mOriginalBase on every call, so the restore has to happen
                // before the next one. compute_conditioning() calls run()
                // twice, and returning early from either branch would let the
                // second call capture the base the first one left behind --
                // after which the restore condition never fires again
                this->restore_indexing_base();

                return aResult ;
            }

//------------------------------------------------------------------------------

            real
            EigenValues::run_arpack( const int_t aJob, const real aSigma )
            {

                real aResult = BELFEM_QUIET_NAN ;
#ifdef BELFEM_ARPACK
                this->link_matrix();

                mInfo.fill( 0  );

                // job 0 asks for the small end of the spectrum, job 1 for
                // the large one, so the count comes from the matching
                // setter rather than from one shared value
                const int_t tNumEigenValues = aJob == 0 ?
                        mNumMinVals : mNumMaxVals ;

                // ncv, maxit and tol are sized here rather than read off a
                // hardwired member. aSigma != 0 marks the folded run, which
                // needs the tight tolerance
                int_t tSubspaceSize = 0 ;
                int_t tNumMaxIter   = 0 ;
                real  tTolerance    = 0.0 ;

                if( ! this->configure( aJob, aSigma != 0.0, tNumEigenValues,
                                       tSubspaceSize, tNumMaxIter, tTolerance ) )
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng: no legal Krylov subspace fits the basis budget at n = %i.\n"
                        "                 The conditioning estimate is not available.\n",
                        ( int ) mN );
                    mSubspaceInfeasible = true ;
                    return BELFEM_QUIET_NAN ;
                }

                // one call site, two drivers. The argument lists and the info
                // layout are identical by design, so only the symbol and the
                // decoders change
                if( mSymmetric )
                {
                    arpack::arpack_symmetric_eigen(
                    &mN,
                    &mNNZ,
                    mK->data(),
                    mK->type() == SpMatrixType::CSC ? mK->rows() : mK->cols(),
                    mK->pointers(),
                    &aJob,
                    &tNumEigenValues,
                    &tSubspaceSize,
                    &tTolerance,
                    &tNumMaxIter,
                    &aSigma,
                    mLambdaReal.data(),
                    mLambdaImag.data(),
                    mInfo.data() ) ;
                }
                else
                {
                arpack::arpack_standard_eigen(
                &mN,
                &mNNZ,
                mK->data(),
                mK->type() == SpMatrixType::CSC ? mK->rows() : mK->cols(),
                mK->pointers(),
                &aJob,
                &tNumEigenValues,
                &tSubspaceSize,
                &tTolerance,
                &tNumMaxIter,
                &aSigma,
                mLambdaReal.data(),
                mLambdaImag.data(),
                mInfo.data() ) ;
                }

                // check error messages. dsaupd and dnaupd do NOT share an
                // error table any more than dnaupd and dneupd do, so the
                // decoder follows the driver that ran
                if( mSymmetric )
                {
                    arpack::check_saupd( mInfo( arpack::gInfoNaupd ) );
                    arpack::check_seupd( mInfo( arpack::gInfoNeupd ) );
                }
                else
                {
                    arpack::check_naupd( mInfo( arpack::gInfoNaupd ) );
                    arpack::check_neupd( mInfo( arpack::gInfoNeupd ) );
                }

                // a non-fatal dnaupd exit still lands here, so the
                // converged count decides whether the value is usable
                const int_t tNumConverged = mInfo( arpack::gInfoNumConverged );


                // non-convergence is an expected algorithmic outcome, not an
                // error: this is a diagnostic, and aborting a production run
                // because a condition number failed to converge turns a
                // recoverable state into a dead job. Same contract as
                // run_parpack -- the caller sees a NaN either way
                if( tNumConverged < tNumEigenValues )
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng [ %s ] converged %i of %i in %i restarts\n"
                        "                 ( ncv = %i, %i products, tol %.1e ).\n",
                        mRunLabel,
                        ( int ) tNumConverged,
                        ( int ) tNumEigenValues,
                        ( int ) mInfo( arpack::gInfoNumIterations ),
                        ( int ) mInfo( arpack::gInfoSubspaceSize ),
                        ( int ) mInfo( arpack::gInfoNumOperations ),
                        ( double ) tTolerance );
                    return BELFEM_QUIET_NAN ;
                }

                // dneupd tops out at nev+1 converged values, which is what
                // the lambda buffers are sized for
                BELFEM_ASSERT( tNumConverged <= ( int_t ) mLambdaReal.length(),
                    "ARPACK-ng reported %i converged values, but only %lu fit into the buffer.",
                    ( int ) tNumConverged,
                    ( long unsigned int ) mLambdaReal.length() );

                // dneupd does not promise an ordering the caller can rely
                // on, so the extremum is picked from the converged set
                // rather than read off the first slot
                aResult = this->reduce_extremum( aJob, tNumConverged );
#else
                BELFEM_ERROR( false,
                    "Can not compute eigenvalues: BELFEM was built without ARPACK-ng.\n"
                    "                 Reconfigure with USE_ARPACK=ON.\n" );
#endif
                return aResult ;
            }

 //------------------------------------------------------------------------------

            real
            EigenValues::run_parpack( const int_t aJob, const real aSigma )
            {
                real aResult = BELFEM_QUIET_NAN ;
#ifdef BELFEM_PARPACK

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // COLLECTIVE from here to the end. Every rank must reach
                // every collective in the same order -- one rank returning
                // early hangs the others
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // only the master holds the full matrix ; a worker's
                // system_matrix() is its assembly submatrix, and its n_rows() is a
                // local dof count rather than N
                SpMatrix * tMatrix        = nullptr ;
                int_t      tNumGlobalRows = 0 ;
                int_t      tNumNonzeros   = 0 ;

                if( mParent->parent()->is_master() )
                {
                    this->link_matrix();

                    tMatrix        = mK ;
                    tNumGlobalRows = mN ;
                    tNumNonzeros   = mNNZ ;
                }

                // PARPACK runs one algorithm redundantly on every rank and
                // stays in lockstep only because the ranks take the same
                // branches. Everything it branches on is taken from the
                // broadcast value, never from the local member -- a divergent
                // nev or job changes the trip count inside pdnaupd and hangs
                broadcast( tNumGlobalRows );
                broadcast( tNumNonzeros );

                int_t tJob = aJob ;
                broadcast( tJob );

                int_t tNumEigenValues = tJob == 0 ? mNumMinVals : mNumMaxVals ;
                broadcast( tNumEigenValues );

                // sigma CHANGES THE OPERATOR, so a divergent value would put
                // the ranks on different problems rather than merely on
                // different schedules. It travels with job and nev for the
                // same reason they do
                real tSigma = aSigma ;
                broadcast( tSigma );

                // ncv, maxit and tol are SIZED ON THE MASTER and broadcast.
                // They must not be computed independently per rank: configure()
                // reads mN, which is a global row count on the master and a
                // local dof count on a worker, so the ranks would size
                // different subspaces and part company inside pdnaupd
                int_t tSubspaceSize = 0 ;
                int_t tNumMaxIter   = 0 ;
                real  tTolerance    = 0.0 ;
                int_t tFeasible     = 0 ;

                if( mParent->parent()->is_master() )
                {
                    tFeasible = this->configure( tJob, tSigma != 0.0, tNumEigenValues,
                                                 tSubspaceSize, tNumMaxIter,
                                                 tTolerance ) ? 1 : 0 ;
                }

                broadcast( tFeasible );
                broadcast( tNumMaxIter );
                broadcast( tSubspaceSize );
                broadcast( tTolerance );

                // collective: every rank leaves together or none does
                if( tFeasible == 0 )
                {
                    if( mCommRank == 0 )
                    {
                        message( InfoLevel::Verbose,
                            "PARPACK: no legal Krylov subspace fits the basis budget at n = %i.\n"
                            "                 The conditioning estimate is not available.\n",
                            ( int ) tNumGlobalRows );
                    }

                    // tFeasible was broadcast, so this is uniform
                    mSubspaceInfeasible = true ;
                    return BELFEM_QUIET_NAN ;
                }

                // the lambda buffers carry the same layout everywhere,
                // because pdneupd returns the same values on every rank
                index_t tNumLambda = tNumEigenValues + 1 ;
                if( mLambdaReal.length() < tNumLambda )
                {
                    mLambdaReal.set_size( tNumLambda );
                    mLambdaImag.set_size( tNumLambda );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // row distribution. The sparsity pattern is scattered only
                // when the structure actually changes -- reset() runs after
                // every Jacobian computation, so it is the wrong hook
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                int_t tRebuild = 0 ;

                if( mParent->parent()->is_master() )
                {
                    tRebuild = ( mDistMatrix    == nullptr
                              || mDistSource     != mK
                              || mDistNumRows    != tNumGlobalRows
                              || mDistNumNonzeros != tNumNonzeros ) ? 1 : 0 ;
                }

                // the decision has to be identical everywhere: the rebuild
                // path is collective
                broadcast( tRebuild );

                if( tRebuild == 1 )
                {
                    if( mDistMatrix != nullptr )
                    {
                        delete mDistMatrix ;
                        mDistMatrix = nullptr ;
                    }

                    // the params are read during construction only, so they
                    // do not outlive this block. Reordering is pinned rather
                    // than left to the default: METIS would make DistMatrix
                    // permute the matrix, which costs a second copy and buys
                    // nothing here ( the permutation is symmetric, so it
                    // would not even move the spectrum )
                    SolverParameters tDistParams( mParent->solver()->type() );
                    tDistParams.set_reordering_method( ReorderingMethod::NATURAL );

                    // NOTE for the += 1 below: the copied pattern is ALWAYS
                    // zero-based, whatever link_matrix() left mK in.
                    // DistMatrix::create_matrix forces
                    // SpMatrixIndexingBase::Cpp before the constructor reaches
                    // distribute_sparsity_pattern(), and distribute_values()
                    // forces it again on every later call. If that ever stops
                    // being true, the shift here silently becomes two-based
                    mDistMatrix = new sparse::DistMatrixCSR< int_t >(
                            &tDistParams, tMatrix );

                    const int_t   tMyRows  = mDistMatrix->n_rows();
                    const int_t * tPtrs    = mDistMatrix->pointers();
                    const int_t * tIndices = mDistMatrix->indices();
                    const int_t   tMyNnz   = tPtrs[ tMyRows ];

                    // DistMatrixCSR hands its arrays out const and builds
                    // them zero-based ; PARPACK wants one-based. Copy once
                    // per structure change, then shift with the vector
                    // operator rather than a loop
                    mParpackPointers.set_size( tMyRows + 1 );
                    std::copy( tPtrs, tPtrs + tMyRows + 1,
                               mParpackPointers.data() );
                    mParpackPointers += 1 ;

                    mParpackIndices.set_size( tMyNnz );
                    std::copy( tIndices, tIndices + tMyNnz,
                               mParpackIndices.data() );
                    mParpackIndices += 1 ;

                    mDistSource       = mK ;
                    mDistNumRows      = tNumGlobalRows ;
                    mDistNumNonzeros  = tNumNonzeros ;
                }

                // collective: rank 0 scatters the values, the others receive
                mDistMatrix->distribute_values( tMatrix );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // solve
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                const int_t tMyRows = mDistMatrix->n_rows();
                const int_t tMyNnz  = mParpackPointers( tMyRows ) - 1 ;

                mInfo.fill( 0 );

                // note on CSC: for UMFPACK / SUPERLU / MUMPS the Jacobian is
                // CSC ( fn_matrix_type.hpp ), so what is distributed here are
                // COLUMN blocks with row indices, and the operator PARPACK
                // sees is the transpose. That is harmless for eigenVALUES,
                // since the spectrum of A and A' is the same -- the serial
                // driver leans on exactly the same equivalence. It would NOT
                // be harmless for eigenvectors: flipping rvec in
                // parpacktools.f90 would silently return the LEFT
                // eigenvectors of A on those solvers
                // mSymmetric is a member and therefore already identical on
                // every rank -- it is set once at setup, never derived from
                // matrix values, so it needs no broadcast
                if( mSymmetric )
                {
                    parpack::parpack_symmetric_eigen(
                            &tMyRows,
                            &tNumGlobalRows,
                            &tMyNnz,
                            mDistMatrix->values(),
                            mParpackIndices.data(),
                            mParpackPointers.data(),
                            &tJob,
                            &tNumEigenValues,
                            &tSubspaceSize,
                            &tTolerance,
                            &tNumMaxIter,
                            &tSigma,
                            mLambdaReal.data(),
                            mLambdaImag.data(),
                            mInfo.data() ) ;
                }
                else
                {
                parpack::parpack_standard_eigen(
                        &tMyRows,
                        &tNumGlobalRows,
                        &tMyNnz,
                        mDistMatrix->values(),
                        mParpackIndices.data(),
                        mParpackPointers.data(),
                        &tJob,
                        &tNumEigenValues,
                        &tSubspaceSize,
                        &tTolerance,
                        &tNumMaxIter,
                        &tSigma,
                        mLambdaReal.data(),
                        mLambdaImag.data(),
                        mInfo.data() ) ;
                }

                // we only print the error on the first proc
                // a crash will cause the whole run to abort
                if ( mCommRank == 0 )
                {
                    if( mSymmetric )
                    {
                        arpack::check_saupd( mInfo( arpack::gInfoNaupd ) );
                        arpack::check_seupd( mInfo( arpack::gInfoNeupd ) );
                    }
                    else
                    {
                        arpack::check_naupd( mInfo( arpack::gInfoNaupd ) );
                        arpack::check_neupd( mInfo( arpack::gInfoNeupd ) );
                    }
                }

                const int_t tNumConverged = mInfo( arpack::gInfoNumConverged );

                // non-convergence is an expected algorithmic outcome, not an
                // error: this is a diagnostic, and aborting a production run
                // because a condition number failed to converge converts a
                // recoverable state into a dead job
                if( tNumConverged < tNumEigenValues )
                {
                    if( mCommRank == 0 )
                    {
                        message( InfoLevel::Verbose,
                            "PARPACK [ %s ] converged %i of %i in %i restarts\n"
                            "                 ( ncv = %i, %i products, tol %.1e ).\n",
                            mRunLabel,
                            ( int ) tNumConverged,
                            ( int ) tNumEigenValues,
                            ( int ) mInfo( arpack::gInfoNumIterations ),
                            ( int ) mInfo( arpack::gInfoSubspaceSize ),
                            ( int ) mInfo( arpack::gInfoNumOperations ),
                            ( double ) tTolerance );
                    }
                    return BELFEM_QUIET_NAN ;
                }

                BELFEM_ASSERT( tNumConverged <= ( int_t ) mLambdaReal.length(),
                    "PARPACK reported %i converged values, but only %lu fit into the buffer.",
                    ( int ) tNumConverged,
                    ( long unsigned int ) mLambdaReal.length() );

                // pdneupd returns the same values on every rank, so this
                // reduction is local and still agrees across ranks
                aResult = this->reduce_extremum( tJob, tNumConverged );
#else
                BELFEM_ERROR( false,
                    "Can not compute eigenvalues in parallel: BELFEM was built without PARPACK.\n"
                    "                 Reconfigure with USE_PARPACK=ON.\n" );
#endif
                return aResult ;
            }

            void
            EigenValues::restore_indexing_base()
            {
                if ( mParent->parent()->is_master() && mOriginalBase == 0 )
                {
                    mK->set_indexing_base( SpMatrixIndexingBase::Cpp );
                }
            }

//------------------------------------------------------------------------------

            real
            EigenValues::reduce_extremum(
                    const int_t aJob,
                    const int_t aNumConverged )
            {
                int_t tWinner = 0 ;

                real tResult = std::sqrt( mLambdaReal( 0 ) * mLambdaReal( 0 )
                    + mLambdaImag( 0 ) * mLambdaImag( 0 ) );

                for( int_t k=1; k<aNumConverged; ++k )
                {
                    real tLambda = std::sqrt( mLambdaReal( k ) * mLambdaReal( k )
                        + mLambdaImag( k ) * mLambdaImag( k ) );

                    if( aJob == 0 ? tLambda < tResult : tLambda > tResult )
                    {
                        tResult = tLambda ;
                        tWinner = k ;
                    }
                }

                // the signed parts come from the slot that WON, not from slot
                // zero and not recomputed: the returned magnitude and the two
                // members must describe one and the same Ritz value, or a
                // guard reads one value while the caller uses another
                mExtremalReal = mLambdaReal( tWinner );
                mExtremalImag = mLambdaImag( tWinner );

                return tResult ;
            }

//------------------------------------------------------------------------------

            real
            EigenValues::compute_smallest_eigenvalues()
            {
                return this->run( 0, 0.0 );
            }

//------------------------------------------------------------------------------

            real
            EigenValues::compute_largest_eigenvalues()
            {
                return this->run( 1, 0.0 );
            }

//------------------------------------------------------------------------------

            real EigenValues::compute_conditioning()
            {
                // An end of the spectrum that has been proven out of reach is
                // not tried again. Retrying it cost 4747 ms PER TIMESTEP on
                // tapestack3d ( n = 97095 ) and produced nothing on any of
                // them. The flag is only ever assigned from a broadcast
                // outcome, so every rank leaves here together
                if( mDiagnosticUnavailable )
                {
                    return BELFEM_QUIET_NAN ;
                }

                if( mMatrixFlag )
                {
                    mParent->compute_jacobian_and_rhs();
                }

                // dneupd hands back an exact zero imaginary part for a real
                // eigenvalue -- it stores a complex pair explicitly -- so this
                // threshold only has to survive roundoff, not resolve anything
                constexpr real tImagTolerance = 1e-8 ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // 1. rho = max |lambda|, the plain exterior request. This is
                //    the end Arnoldi is good at: measured 5 restarts / 61
                //    products against 122 / 1231 for the small end
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                mRunLabel = "largest |lambda|" ;

                const real tRho = this->run( 1, 0.0 );

                const int_t tRhoProducts = mInfo( arpack::gInfoNumOperations );

                // Every branch below is taken from a value the MASTER computed
                // and then broadcast. It must not be derived rank-locally: in a
                // multi-rank build without PARPACK only the master runs the
                // algorithm and only the scalar result is broadcast, so a
                // worker's Ritz buffers -- and therefore mExtremalReal /
                // mExtremalImag -- are stale by construction
                int_t tOutcome = ( int_t ) EigenOutcome::Ok ;

                if( mParent->parent()->is_master() )
                {
                    if( std::isnan( tRho ) )
                    {
                        // an infeasible subspace is not a failed iterate: no
                        // run of this size can ever be attempted, so it latches
                        // at once rather than spending its strikes
                        tOutcome = mSubspaceInfeasible ?
                                ( int_t ) EigenOutcome::Infeasible :
                                ( int_t ) EigenOutcome::NotConverged ;
                    }
                    else if( std::abs( mExtremalImag )
                             > tImagTolerance * std::abs( mExtremalReal ) )
                    {
                        // a complex extremal value means the operator is not
                        // the symmetric one this path assumes, and the ratio
                        // would not be kappa_2 even if the small end converged
                        tOutcome = ( int_t ) EigenOutcome::Complex ;
                    }
                }

                broadcast( tOutcome );

                if( tOutcome != ( int_t ) EigenOutcome::Ok )
                {
                    return this->report_unavailable( tOutcome );
                }

                // reported HERE rather than only on full success: when the fold
                // fails, the run used to say nothing at all about the end that
                // worked, so a failing log could not tell anyone what
                // lambda_max was or how cheap the first end is
                if( mCommRank == 0 )
                {
                    message( InfoLevel::Verbose,
                        "%s: lambda_max %.4e ( %i products ).\n",
                        mBackendLabel.c_str(),
                        ( double ) tRho,
                        ( int ) tRhoProducts );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // 2. the small end, by SHIFT-INVERT.
                //
                //    The spectral FOLD that used to live here is gone, and the
                //    reason is measured rather than argued. A shift is affine,
                //    so it preserves ratios of differences: folding about
                //    sigma > max|lambda| makes lambda_min the DOMINANT
                //    eigenvalue of the operator -- curing the filter problem
                //    that regular-mode 'SM' has -- but leaves its SEPARATION
                //    from lambda_2 exactly as it was. Separation is what sets
                //    the convergence rate.
                //
                //    On the tapestack3d thermal Jacobian ( n = 97095,
                //    symmetric, SPD, kappa_2 = 2.9635e7 ) the six smallest
                //    eigenvalues agree to three significant figures, giving a
                //    small-end relative gap of 9.5e-10 -- five orders below
                //    what a polynomial Krylov method can resolve. The fold
                //    spent 572 restarts and 6303 products confirming it.
                //
                //    Inversion is RATIONAL, not affine, and does change the
                //    separation: the same gap becomes 2.8e-2, an amplification
                //    of 3e7. That is the whole difference between converging
                //    in tens of products and never converging.
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                mRunLabel = "shift-invert smallest lambda" ;

                real tLambdaMin = this->run_shift_invert( tOutcome );

                // run_shift_invert broadcasts its own outcome, so every rank
                // reaches this branch with the same value
                if( tOutcome == ( int_t ) EigenOutcome::Ok )
                {
                    if( std::isnan( tLambdaMin ) )
                    {
                        tOutcome = ( int_t ) EigenOutcome::NotConverged ;
                    }
                    else if( tLambdaMin <= 0.0 )
                    {
                        // the matrix is not positive definite. Reported rather
                        // than returned as a ratio: lambda_max / lambda_min is
                        // only kappa_2 for a symmetric POSITIVE DEFINITE
                        // matrix, and a negative lambda_min would make the
                        // printed number meaningless rather than merely
                        // imprecise
                        tOutcome = ( int_t ) EigenOutcome::NotPositiveDefinite ;
                    }
                    else if( tLambdaMin > tRho )
                    {
                        // lambda_min <= lambda_max <= rho always
                        tOutcome = ( int_t ) EigenOutcome::Inconsistent ;
                    }
                }

                if( tOutcome != ( int_t ) EigenOutcome::Ok )
                {
                    return this->report_unavailable( tOutcome );
                }

                mOutcome      = EigenOutcome::Ok ;
                mFailureCount = 0 ;

                // lambda_min > 0 has now been PROVEN, so every eigenvalue lies
                // in ( 0, rho ] and lambda_max is rho itself
                const real tKappa = tRho / tLambdaMin ;

                // reporting only, rank 0. message() does not filter by rank --
                // every other message in this file guards it explicitly -- and
                // under PARPACK configure() runs on the master alone, so
                // mEffectiveTolerance is only meaningful here
                if( mCommRank == 0 )
                {
                    // the achieved product count, not a modeled budget: this is
                    // the number that tells anyone re-tuning the heuristic what
                    // the folded end actually costs on THIS matrix
                    message( InfoLevel::Verbose,
                        "%s: lambda_min %.4e ( %i products ), kappa %.3e.\n",
                        mBackendLabel.c_str(),
                        ( double ) tLambdaMin,
                        ( int ) mInfo( arpack::gInfoNumOperations ),
                        ( double ) ( tRho / tLambdaMin ) );

                    // lambda_min falls out of a cancellation, so its relative
                    // error is about max( tol, eps_mach ) * kappa. Past that the
                    // digits are gone and the number is a resolution-limited
                    // estimate -- NOT a proved bound, and said plainly rather
                    // than dressed up. The tolerance used is the one configure()
                    // actually applied, which an explicit set_tolerance()
                    // changes
                    if( tKappa * mEffectiveTolerance > 0.1 )
                    {
                        message( InfoLevel::Verbose,
                            "%s: kappa ~ %.2e is at the resolution limit of the folded\n"
                            "                 estimate ( tol %.1e x kappa ). Read it as an order of magnitude.\n",
                            mBackendLabel.c_str(),
                            ( double ) tKappa,
                            ( double ) mEffectiveTolerance );
                    }
                }

                return tKappa ;
            }

//------------------------------------------------------------------------------

            namespace
            {
                /**
                 * drops a frozen factorization on the way out of a scope,
                 * whatever the way out is. The ordinary failure break already
                 * unfreezes explicitly ; this exists for the paths that do not
                 * run to the end of the function -- a hard error inside a
                 * frozen solve, or freeze_factorization() itself throwing on
                 * its identity check. Leaving the scope armed there would make
                 * an unrelated later solve on the same wrapper silently reuse
                 * dead factors
                 */
                class FrozenScopeGuard
                {
                    solver::Wrapper * mWrapper = nullptr ;

                public:

                    FrozenScopeGuard() = default ;

                    ~FrozenScopeGuard()
                    {
                        if( mWrapper != nullptr )
                        {
                            mWrapper->unfreeze_factorization();
                        }
                    }

                    void
                    arm( solver::Wrapper * aWrapper )
                    {
                        mWrapper = aWrapper ;
                    }

                    void
                    release()
                    {
                        mWrapper = nullptr ;
                    }
                };
            }

//------------------------------------------------------------------------------

            real
            EigenValues::run_shift_invert( int_t & aOutcome )
            {
                aOutcome = ( int_t ) EigenOutcome::Ok ;

#ifdef BELFEM_MUMPS
                const bool tIsMaster = mParent->parent()->is_master() ;

                // ---- sizing. Master decides, everyone is told: mN is a
                // global row count on the master and a LOCAL dof count on a
                // worker, so a rank that sized this itself would size a
                // different problem
                int_t tSubspaceSize = 0 ;
                int_t tNumMaxIter   = 0 ;
                real  tTolerance    = 0.0 ;
                int_t tNumRows      = 0 ;
                int_t tFeasible     = 0 ;

                if( tIsMaster )
                {
                    tNumRows  = mN ;
                    tFeasible = this->configure( 1, true, mNumMinVals,
                                                 tSubspaceSize, tNumMaxIter,
                                                 tTolerance ) ? 1 : 0 ;
                }

                broadcast( tFeasible );
                broadcast( tNumRows );
                broadcast( tSubspaceSize );
                broadcast( tNumMaxIter );
                broadcast( tTolerance );

                if( tFeasible == 0 )
                {
                    aOutcome = ( int_t ) EigenOutcome::Infeasible ;
                    return BELFEM_QUIET_NAN ;
                }

                const int_t tNumEigenValues = mNumMinVals ;
                const int_t tLdv            = tNumRows ;
                const int_t tLworkl         = tSubspaceSize * ( tSubspaceSize + 8 ) ;

                // ---- the dedicated solver. Built here rather than reusing
                // mSolver: this one must be MUMPS whatever the field's
                // production solver is, because MUMPS is the wrapper that can
                // freeze a factorization. Every rank builds it -- Solver's
                // constructor and first solve are collective
                if( mShiftInvertSolver == nullptr )
                {
                    mShiftInvertSolver = new Solver( SolverType::MUMPS );

                    // UNSYMMETRIC ( SYM = 0 ), even though the matrix is
                    // symmetric and MUMPS could halve the factor memory.
                    //
                    // BELFEM stores the FULL matrix and hands it to MUMPS
                    // whole: nothing anywhere extracts a triangle
                    // ( cl_SolverMUMPS.cpp passes SymmetryMode straight into
                    // SYM ). MUMPS with SYM = 1 or 2 expects only the lower
                    // triangle, so a full matrix double-counts every
                    // off-diagonal and the factorization is of a DIFFERENT
                    // matrix. It does not fail -- it returns a plausible
                    // wrong answer.
                    //
                    // Measured 2026-08-29 on a 1D Laplacian with an analytic
                    // spectrum: SYM = 2 gave a first-solve residual
                    // ||Ax-b||/||b|| = 7.4e16 and a converged lambda_min of
                    // -2.5e-19 against a true 2.46e-6, with ARPACK reporting
                    // info 0 and nconv 1 throughout. SYM = 0 on the same
                    // matrix gives residual 1.3e-12 and lambda_min correct to
                    // 6e-12.
                    //
                    // This is the whole production convention too -- every
                    // BELFEM solve runs SYM = 0. Do not "optimize" it without
                    // first making the caller supply a triangle
                    mShiftInvertSolver->set_symmetry_mode(
                            SymmetryMode::Unsymmetric );
                }

                // A DIAGNOSTIC MUST NOT KILL THE RUN. Without this the wrapper
                // raises an always-active error on a factorization failure and
                // the whole job dies to report a conditioning number ; with it
                // the failure becomes a flag, this function reports
                // SolverFailed, and the solve the user actually cares about
                // carries on. Set on every rank, before any solve
                mShiftInvertSolver->wrapper()->set_soft_fail( true );

                // the latch is sticky by design, so a failure recorded by an
                // earlier call would make every later one look failed
                mShiftInvertSolver->wrapper()->clear_failure();

                // ---- ARPACK state. Master-only storage: a worker never
                // enters the Fortran side
                if( tIsMaster )
                {
                    mSiResid.set_size( tNumRows, 0.0 );
                    mSiBasis.set_size( tLdv * tSubspaceSize, 0.0 );
                    mSiWorkD.set_size( 3 * tNumRows, 0.0 );
                    mSiWorkL.set_size( tLworkl, 0.0 );
                    mSiIparam.set_size( 11, 0 );
                    mSiIpntr.set_size( 11, 0 );
                    mSiLambda.set_size( tNumEigenValues, 0.0 );

                    mSiIparam( 0 ) = 1 ;            // exact shifts, so ido = 3 cannot occur
                    mSiIparam( 2 ) = tNumMaxIter ;  // restart cap
                    mSiIparam( 6 ) = 3 ;            // MODE 3: shift-invert
                }

                // the solve vectors are needed on EVERY rank: the solve is
                // collective and each rank passes its own pair, even though
                // only the master's carry data
                mSiRhs.set_size( tNumRows, 0.0 );
                mSiLhs.set_size( tNumRows, 0.0 );

                // ... and so is a dereferenceable matrix. MUMPS ignores
                // worker-side storage, but the collective solve below still
                // does *mK on every rank. The capture is UNCONDITIONAL: a
                // once-captured mK can be stale after SolverData::reset()
                // deletes the system matrix -- on any rank, since
                // link_matrix() runs unguarded in compute_matrices() and
                // run_arpack() too, not only on the master. The worker's own
                // assembly submatrix is exactly what the production solve
                // passes on that rank, so it is used here too
                mK = mParent->system_matrix() ;

                BELFEM_ERROR( mK != nullptr,
                    "EigenValues::run_shift_invert() : system_matrix() returned null - "
                    "the system matrix has been reset and not rebuilt" );

                int_t tIdo   = 0 ;
                int_t tInfo  = 0 ;
                real  tSigma = 0.0 ;   // shift zero: we want the eigenvalues nearest the origin

                // ---- the reverse-communication loop.
                // Master decides; the command word is broadcast so that every
                // rank runs the SAME sequence of collective solves. No rank
                // may take this branch on its own evidence
                // THREE-valued, not two. The workers never see tIdo -- it is
                // master-only state -- so every branch that decides whether a
                // COLLECTIVE happens must be decided from the broadcast word.
                // An earlier version folded ido = 2 in with the solve requests
                // and deadlocked: the master took its copy branch while every
                // worker sat inside the collective solve waiting for a master
                // that was never coming
                constexpr int_t tCmdStop = 0 ;
                constexpr int_t tCmdSolve = 1 ;
                constexpr int_t tCmdCopy = 2 ;

                FrozenScopeGuard tGuard ;
                bool  tFrozen  = false ;
                int_t tCommand = tCmdStop ;

                while( true )
                {
                    if( tIsMaster )
                    {
                        arpack::arpack_si_step(
                                &tIdo,
                                &tNumRows,
                                &tNumEigenValues,
                                &tSubspaceSize,
                                &tTolerance,
                                mSiResid.data(),
                                mSiBasis.data(),
                                &tLdv,
                                mSiIparam.data(),
                                mSiIpntr.data(),
                                mSiWorkD.data(),
                                mSiWorkL.data(),
                                &tLworkl,
                                &tInfo );

                        tCommand = ( tIdo == -1 || tIdo == 1 ) ? tCmdSolve :
                                     tIdo == 2                   ? tCmdCopy  :
                                                                   tCmdStop ;
                    }

                    broadcast( tCommand );

                    if( tCommand == tCmdStop )
                    {
                        break ;
                    }

                    if( tCommand == tCmdCopy )
                    {
                        // ---- ido = 2 is y = B*x, and B is the identity here,
                        // so it reduces to a copy of master-local data. NO
                        // collective: every rank must agree to skip the solve,
                        // which is why this is reached through the broadcast
                        // word and not through tIdo
                        if( tIsMaster )
                        {
                            const int_t tX = mSiIpntr( 0 ) - 1 ;  // ipntr is 1-based
                            const int_t tY = mSiIpntr( 1 ) - 1 ;

                            for( int_t k=0; k<tNumRows; ++k )
                            {
                                mSiWorkD( tY + k ) = mSiWorkD( tX + k );
                            }
                        }

                        continue ;
                    }

                    {
                        // ---- the inverse application, y = A^-1 * x.
                        //
                        // THE OPERAND SLOT DEPENDS ON ido, and getting it
                        // wrong corrupts the basis silently rather than
                        // failing: on ido = 1 in mode 3 the operand is B*x at
                        // ipntr( 3 ), NOT x at ipntr( 1 ). arpack-ng 3.9.1
                        // dsaupd.f says so in as many words -- "In mode 3,4
                        // and 5, the vector B * X is already available in
                        // WORKD(ipntr(3))"
                        if( tIsMaster )
                        {
                            const int_t tSrc = ( tIdo == 1 ? mSiIpntr( 2 )
                                                           : mSiIpntr( 0 ) ) - 1 ;

                            for( int_t k=0; k<tNumRows; ++k )
                            {
                                mSiRhs( k ) = mSiWorkD( tSrc + k );
                            }
                        }

                        // COLLECTIVE. The first one runs JOB 6 and builds the
                        // factorization ; every later one runs JOB 3 against
                        // it. There is no dummy solve -- the factorizing solve
                        // IS the first operator application
                        mShiftInvertSolver->solve( *mK, mSiLhs, mSiRhs );

                        int_t tSolveFailed =
                                mShiftInvertSolver->wrapper()->failed() ? 1 : 0 ;
                        broadcast( tSolveFailed );

                        if( tSolveFailed == 1 )
                        {
                            aOutcome = ( int_t ) EigenOutcome::SolverFailed ;
                            break ;
                        }

                        // freeze on the way out of the FIRST successful solve,
                        // so every later application reuses these factors
                        if( ! tFrozen )
                        {
                            // guard armed BEFORE the freeze: if
                            // freeze_factorization() throws on its identity
                            // check, the scope must still be dropped on the
                            // way out
                            tGuard.arm( mShiftInvertSolver->wrapper() );

                            mShiftInvertSolver->wrapper()
                                    ->freeze_factorization( *mK );
                            tFrozen = true ;
                        }

                        if( tIsMaster )
                        {
                            const int_t tY = mSiIpntr( 1 ) - 1 ;

                            for( int_t k=0; k<tNumRows; ++k )
                            {
                                mSiWorkD( tY + k ) = mSiLhs( k );
                            }
                        }
                    }
                }

                // ---- the scope ends here, on EVERY path including the
                // failure break above. Nothing below may run a frozen solve
                if( tFrozen )
                {
                    // RELEASE FIRST. If the explicit unfreeze threw while the
                    // guard was still armed, unwinding would call the guard's
                    // destructor, which would call unfreeze again -- a second
                    // exception during unwinding terminates the process. The
                    // MUMPS override cannot throw today, so this is latent
                    // rather than live, but the ordering costs nothing and the
                    // failure mode is unrecoverable
                    tGuard.release();
                    mShiftInvertSolver->wrapper()->unfreeze_factorization();
                }

                // ---- extraction and outcome, decided on the master and
                // broadcast before anyone branches on it
                real tLambdaMin = BELFEM_QUIET_NAN ;

                if( aOutcome == ( int_t ) EigenOutcome::Ok )
                {
                    if( tIsMaster )
                    {
                        if( tInfo < 0 )
                        {
                            aOutcome = ( int_t ) EigenOutcome::NotConverged ;
                        }
                        else
                        {
                            int_t tExtractInfo = 0 ;

                            arpack::arpack_si_extract(
                                    &tSigma,
                                    &tNumRows,
                                    &tNumEigenValues,
                                    &tSubspaceSize,
                                    &tTolerance,
                                    mSiResid.data(),
                                    mSiBasis.data(),
                                    &tLdv,
                                    mSiIparam.data(),
                                    mSiIpntr.data(),
                                    mSiWorkD.data(),
                                    mSiWorkL.data(),
                                    &tLworkl,
                                    mSiLambda.data(),
                                    &tExtractInfo );

                            mInfo.fill( 0 );
                            mInfo( arpack::gInfoNaupd )         = tInfo ;
                            mInfo( arpack::gInfoNeupd )         = tExtractInfo ;
                            mInfo( arpack::gInfoNumConverged )  = mSiIparam( 4 ) ;
                            mInfo( arpack::gInfoNumIterations ) = mSiIparam( 2 ) ;
                            mInfo( arpack::gInfoNumOperations ) = mSiIparam( 8 ) ;
                            mInfo( arpack::gInfoSubspaceSize )  = tSubspaceSize ;

                            if( tExtractInfo != 0
                                || mSiIparam( 4 ) < tNumEigenValues )
                            {
                                aOutcome = ( int_t ) EigenOutcome::NotConverged ;
                            }
                            else
                            {
                                // ALREADY eigenvalues of A: dseupd applied
                                // lambda = 1/theta + sigma itself. A
                                // reciprocal here would invert them twice
                                tLambdaMin = mSiLambda( 0 );

                                for( int_t k=1; k<mSiIparam( 4 )
                                             && k<( int_t ) mSiLambda.length(); ++k )
                                {
                                    if( mSiLambda( k ) < tLambdaMin )
                                    {
                                        tLambdaMin = mSiLambda( k );
                                    }
                                }
                            }
                        }
                    }
                }

                broadcast( aOutcome );
                broadcast( tLambdaMin );

                // the MUMPS instance is deliberately NOT freed here. This runs
                // once per timestep, and a per-call free/re-init burns a fresh
                // slot in the mumpstools registry each time ( the production
                // solver keeps the pool from ever draining, so slot IDs used
                // to climb until the pool ran dry -- observed as the
                // tapestack3d step-8 abort ). Keeping the instance also lets
                // select_job() reuse the analysis ( JOB 5 ) on the next step.
                // The factors of the last step stay resident until
                // ~EigenValues; that memory cost is the documented price of
                // the opt-in diagnostic

                return tLambdaMin ;
#else
                // no MUMPS in this build, so there is no wrapper that can hold
                // a factorization for repeated solves. Reported, not faked
                aOutcome = ( int_t ) EigenOutcome::Infeasible ;
                return BELFEM_QUIET_NAN ;
#endif
            }

//------------------------------------------------------------------------------

            real
            EigenValues::report_unavailable( const int_t aOutcome )
            {
                mOutcome = ( EigenOutcome ) aOutcome ;

                // Only Infeasible latches on its first occurrence, because it
                // is the one outcome that depends on the problem SIZE and the
                // basis budget rather than on the numbers in the matrix -- it
                // cannot come right at the next timestep. That rests on the dof
                // count being fixed for the lifetime of the DofManager, which
                // it is today; if dynamic remeshing or reinitialization ever
                // changes it mid-run, this latch has to be cleared or keyed by
                // size, because link_matrix() would refresh mN and the latch
                // would stop the refresh from ever being acted on.
                //
                // Everything else gets the strike count, INCLUDING
                // NotPositiveDefinite. An earlier version latched that one at
                // once, on the claim that a non-positive spectrum is "a
                // property of the problem, not of this timestep". That claim is
                // false: the Jacobian is rebuilt every timestep, the deck is
                // nonlinear and delta t varies, so the spectrum moves with it.
                // Nothing says an indefinite iterate at step k implies one at
                // step k+1, and a single transient sign flip would have
                // switched the diagnostic off for the rest of the run.
                // A spectrum that really is indefinite by construction still
                // latches -- it just takes mMaxFailures steps to prove it
                // instead of being assumed
                ++mFailureCount ;

                const bool tLatch =
                        mOutcome == EigenOutcome::Infeasible
                     || mFailureCount >= mMaxFailures ;

                if( mCommRank == 0 )
                {
                    // EVERY enumerator gets its own arm, and the default names
                    // the unmatched value instead of reading as a diagnosis.
                    // Until 2026-08-29 Infeasible had no arm: it WAS the
                    // default, so SolverFailed -- added to the enum later --
                    // inherited the basis-budget sentence, and every failed or
                    // exhausted solver was reported to the user as "no legal
                    // Krylov subspace fits the basis budget". Wrong subsystem,
                    // in the one message a user reads when the diagnostic goes
                    // quiet. A default that cannot be mistaken for a reason is
                    // what stops the next enumerator repeating it
                    char tUnhandled[ 96 ];
                    std::snprintf( tUnhandled, sizeof( tUnhandled ),
                        "unhandled eigen outcome %i - a missing arm, not a diagnosis",
                        ( int ) mOutcome );

                    const char * tReason =
                        mOutcome == EigenOutcome::NotConverged        ?
                            "ARPACK did not converge the requested value" :
                        mOutcome == EigenOutcome::Complex             ?
                            "the extremal eigenvalue is complex, so the spectral ratio is undefined" :
                        mOutcome == EigenOutcome::NotPositiveDefinite ?
                            // NOT "no spectral fold can reach it": the fold was
                            // retired in favour of shift-invert, which CAN
                            // return a negative lambda_min. What the outcome
                            // now means is that the value came back and was not
                            // positive, so the ratio is not a condition number
                            "the smallest eigenvalue came back non-positive, so the spectrum is not\n"
                            "                 positive definite and the ratio is not a condition number" :
                        mOutcome == EigenOutcome::Inconsistent        ?
                            "the fold returned lambda_min > lambda_max, so the spectrum is not real" :
                        mOutcome == EigenOutcome::Infeasible          ?
                            "no legal Krylov subspace fits the basis budget" :
                        mOutcome == EigenOutcome::SolverFailed        ?
                            // deliberately NOT "the factorization failed":
                            // wrapper()->failed() is raised by a failed or
                            // registry-exhausted CREATE as well as by a failed
                            // factorization or frozen solve, so naming factors
                            // would point at the wrong subsystem a second time
                            "the linear solver behind the shift-invert could not supply an\n"
                            "                 instance or complete a solve" :
                            tUnhandled ;

                    // Detailed, NOT Verbose. InfoLevel::Detailed is 3 and
                    // Verbose is 4, message() prints only when the level is
                    // <= the logger's, and both debug executables construct
                    // gLog at Detailed -- so a Verbose reason line does not
                    // print in the configuration where somebody is actually
                    // debugging a quiet diagnostic. The per-step ARPACK
                    // telemetry around this function is correctly Verbose;
                    // this is the conclusion, not telemetry
                    if( tLatch )
                    {
                        message( InfoLevel::Detailed,
                            "%s: %s.\n"
                            "                 The conditioning diagnostic is switched OFF for the rest of\n"
                            "                 this run rather than retried every timestep.\n",
                            mBackendLabel.c_str(),
                            tReason );
                    }
                    else
                    {
                        message( InfoLevel::Detailed,
                            "%s: %s.\n"
                            "                 No conditioning estimate for this step ( strike %i of %i ).\n",
                            mBackendLabel.c_str(),
                            tReason,
                            ( int ) mFailureCount,
                            ( int ) mMaxFailures );
                    }
                }

                mDiagnosticUnavailable = tLatch ;

                return BELFEM_QUIET_NAN ;
            }

//------------------------------------------------------------------------------

            real
            EigenValues::compute_lambda_max()
            {
                this->compute_matrices();
                mMatrixFlag = true ;

                if( mParent->parent()->is_master() )
                {
                    if( mFirstRun )
                    {
                        mX.set_size( mN );
                        mY.set_size( mN );
                        mZ.set_size( mN );
                        mXreal.set_size( mN );
                        mYreal.set_size( mN );
                        mZreal.set_size( mN );


                        random_seed();
                        for( int_t k=0; k<mN; ++k )
                        {
                            mXreal( k ) = rand();
                        }

                        mK->multiply( mXreal, mYreal );
                        mXreal /= std::sqrt( std::abs( dot( mXreal, mYreal ) ) );

                        mXimag.set_size( mN, 0.0 );
                        mYimag.set_size( mN, 0.0 );
                        mZimag.set_size( mN );

                        for( int_t i=0; i<mN; ++i )
                        {
                            mX( i ) = cplx( mXreal( i ), mXimag( i ) );
                        }

                        mSolver = new Solver( mParent->solver()->type() );
                        mFirstRun = false ;

                    }

                    mLambdaMax = 0 ;
                }
                else
                {
                    if( mFirstRun )
                    {
                        mSolver = new Solver( mParent->solver()->type() );
                        mFirstRun = false ;
                    }
                }


                real tErr = 1.0 ;

                for( int_t k=0; k < mNumMaxIter && tErr > mEpsilon ; ++k )
                {
                    real tNormZi = 0.0 ;
                    if( mParent->parent()->is_master() )
                    {
                        mM->multiply( mXreal, mZreal );



                        if( norm( mXimag ) > BELFEM_EPSILON )
                        {
                            mM->multiply( mXimag, mZimag );
                            tNormZi = norm( mZimag );
                        }
                    }

                    broadcast( tNormZi );

                    mSolver->solve( *mM, mYreal, mZreal );

                    if( tNormZi > BELFEM_EPSILON )
                    {
                        mSolver->solve( *mM, mYimag, mZimag );
                    }
                    else
                    {
                        mYimag.fill( 0.0 );
                    }

                    if( mParent->parent()->is_master() )
                    {
                        mK->multiply( mXreal, mZreal );
                        mK->multiply( mXimag, mZimag ) ;

                        for( int_t i=0; i<mN; ++i )
                        {
                            mY( i ) = cplx( mYreal( i ), mYimag( i ) );
                        }

                        for( int_t i=0; i<mN; ++i )
                        {
                            mZ( i ) = cplx( mZreal( i ), mZimag( i ) );
                        }

                        for( int_t i=0; i<mN; ++i )
                        {
                            mXreal( i ) = std::real( mX( i ) ) ;
                        }
                        for( int_t i=0; i<mN; ++i )
                        {
                            mXimag( i ) = std::imag( mX( i ) ) ;
                        }


                        mK->multiply( mXreal, mYreal );
                        mK->multiply( mXimag, mYimag );

                        mM->multiply( mXreal, mZreal );
                        mM->multiply( mXimag, mZimag );

                        for( int_t i=0; i<mN; ++i )
                        {
                            mY( i ) = std::complex< real >( mYreal( i ), mYimag( i ) );
                        }

                        for( int_t i=0; i<mN; ++i )
                        {
                            mZ( i ) = std::complex< real >( mZreal( i ), mZimag( i ) );
                        }

                        mX = mY / std::sqrt( dot( mX, mZ ) );

                        cplx tLambda = dot( mX, mY ) / dot( mX, mZ );

#ifdef BELFEM_ARMADILLO
                        tErr = norm(  mY.vector_data() - tLambda * mZ.vector_data() );
#else
                        tErr = 0.0 ;
                        for( int_t i=0; i<mN; ++i )
                        {
                            tErr += std::abs( mY( i ) - tLambda * mZ( i ) );
                        }
                        tErr = std::sqrt( tErr );
#endif
                        mLambdaMax = std::abs( tLambda );
                    }

                    comm_barrier();
                    broadcast( tErr );
                    broadcast( mLambdaMax );

                    if( tErr < mEpsilon )
                    {
                        break ;
                    }
                }
                return mLambdaMax ;
            }

//------------------------------------------------------------------------------

            void
            EigenValues::reset()
            {
                mMatrixFlag = false ;
            }

//------------------------------------------------------------------------------
        }
    }
}
