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

#ifdef BELFEM_PETSC
#include <algorithm>   // for create_indices
#endif
#include "commtools.hpp"
#include "cl_Logger.hpp"

#include "petsctools.hpp"
#include "cl_SolverPETSC.hpp"
#include "st_SolverPetscData.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        PETSC::PETSC( const SolverParameters * aParams ) :
                Wrapper( "PETSC    ", true ),
                mParams( aParams )
        {
#ifdef BELFEM_PETSC

            mPreconditioner = mParams->preconditioner() ;
            mKrylovMethod   = mParams->krylov_method() ;
            mEpsilon        = mParams->relative_tolerance() ;

            // use world as default communicator
            mData.mComm = PETSC_COMM_WORLD ;
#endif
        }

//------------------------------------------------------------------------------

        PETSC::~PETSC()
        {
            this->free();
        }

//------------------------------------------------------------------------------

        void
        PETSC::set(
                const Preconditioner aPreconditioner,
                const KrylovMethod   aKrylovMethod,
                const real           aEpsilon )
        {
            mPreconditioner = aPreconditioner ;
            mKrylovMethod   = aKrylovMethod ;
            mEpsilon        = aEpsilon ;
        }

//------------------------------------------------------------------------------

        void
        PETSC::solve(
                SpMatrix      & aMatrix,
                Vector <real> & aLHS,
                Vector <real> & aRHS )
        {
#ifdef BELFEM_PETSC

            Vector< index_t > tLengths( 2 );
            if( gComm.rank() == 0 )
            {
                // enforce zero based indexing on matrix
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

                tLengths( 0 ) = aRHS.length() ;
                tLengths( 1 ) = aLHS.length() ;
            }

            broadcast( tLengths );

            index_t tRhsLength = tLengths( 0 );

            // bool set flag if a nonzero LHS is given
            bool tHaveLHS = tLengths( 0 ) == tLengths( 1 );

            // no LHS of matching length: allocate it ( it then cannot serve as initial guess )
            if( ! tHaveLHS )
            {
                // allocate LHS
                aLHS.set_size( tRhsLength, 0.0 );
            }

            // on non-master, also allocate RHS
            if( tRhsLength != aRHS.length() )
            {
                aRHS.set_size( tRhsLength );
            }

            // check if we need to initialize the solver
            if( ! this->is_initialized() )
            {
                this->initialize( aMatrix );
            }

            // set solver options (must be done AFTER initialization)
            this->set_preconditioner( mPreconditioner );
            this->set_krylovmethod( mKrylovMethod );
            this->set_matrix_ordering( mParams->reordering_method() );

            // set initial guess flag: respect SolverParameters setting if use_initial_guess is true,
            // otherwise determine based on whether LHS was provided
            if( mParams->use_initial_guess() && tHaveLHS )
            {
                this->set_initial_guess_flag( true );
            }
            else
            {
                this->set_initial_guess_flag( false );
            }

            // update data containers
            if( this->comm_size() > 1 ) // parallel mode
            {
                // distribute matrix values
                mDistMatrix->distribute_values( &aMatrix );

                // update matrix in PETSc
                petsctools_update_matrix( mDistMatrix, mData.mMat );

                if( tHaveLHS )
                {
                    mDistMatrix->distribute_lhs( aLHS );

                    // set local portion of LHS vector
                    petsctools_set_local_vector(
                        mDistMatrix->lhs_vector(),
                        mData.mLHS );
                }
                else
                {
                    // Zero the LHS vector explicitly when no initial guess
                    VecZeroEntries( mData.mLHS );
                }
                mDistMatrix->distribute_rhs( aRHS );

                // set local portion of RHS vector
                petsctools_set_local_vector(
                    mDistMatrix->rhs_vector(),
                    mData.mRHS );
            }
            else  // sequential mode
            {
                // use LHS as initial guess
                if( tHaveLHS )
                {
                    petsctools_set_vector( aLHS, mData.mVectorIndices, mData.mLHS );
                }

                // update RHS data
                petsctools_set_vector( aRHS, mData.mVectorIndices, mData.mRHS );
            }

             // solve the system
            PetscErrorCode tStatus = KSPSolve( mData.mKSP, mData.mRHS, mData.mLHS );

            proc_t tBadProc = 0 ;

            if ( this->rank() > 0 )
            {
                send( tStatus );
                broadcast( tBadProc );
            }
            else
            {
                Vector< PetscErrorCode > tErrors( this->comm_size() );
                collect( tErrors, tStatus );

                for ( proc_t p=0 ; p<this->comm_size() ; ++p )
                {
                    if ( tErrors( p ) != 0 )
                    {
                        tStatus = tErrors( p ) ;
                        tBadProc = p ;
                        break ;
                    }
                }
                broadcast( tBadProc );
            }

            if ( this->rank() == tBadProc )
            {
                // check error
                BELFEM_ERROR( tStatus == 0,
                           "PETSC has thrown the error %i in function KSPSolve()\n%s",
                            tStatus,
                            petsctools_error_message( tStatus ).c_str() );
            }

            // The error code above only says the CALL did not fail. Whether the
            // Krylov iteration converged is reported through the converged
            // reason, which nothing read until now — a solve that diverged or
            // stopped at maxits returned code 0 and was accepted silently, the
            // same defect class as the MUMPS INFO(1)=+1 escalation.
            // An unconverged or zero-iteration solve fed into the nonlinear
            // loop is the leading hypothesis for the thermal Picard freeze, so
            // the iteration count is reported for the record: a stalled
            // nonlinear residual with "0 iterations" printed here means the
            // initial guess already met the KSP relative tolerance and the
            // solve returned its own input.
            // initialized in case KSPGetConvergedReason itself fails — an
            // uninitialized reason feeding the branch below would be UB.
            // The call's own failure is folded into the HARD verdict below
            // rather than raised here: a per-rank error immediately before
            // a collective is exactly the split this block exists to close
            KSPConvergedReason tReason = KSP_CONVERGED_ITERATING ;
            PetscErrorCode tReasonStatus =
                KSPGetConvergedReason( mData.mKSP, & tReason );

            PetscInt tNumIterations = 0 ;
            KSPGetIterationNumber( mData.mKSP, & tNumIterations );

            // A conservative subset of the divergence reasons is treated as
            // an EXPECTED algorithmic failure ( budget exhausted, residual
            // grew, Krylov breakdown — NULL is GMRES's Hessenberg
            // breakdown, same family ): the controller above has a retry
            // policy ( cut Δt, reassemble, solve again ), and the KSP
            // object is safely reusable — the next attempt re-updates the
            // matrix, which triggers a preconditioner rebuild. Everything
            // else — NaN/Inf, preconditioner failures, and any reason not
            // listed — stays HARD: reuse of the KSP after those is not
            // established ( it would need KSPReset plus re-setting the
            // operators ), so they abort until someone wires that path
            const bool tSoftReason =
                   ( tReason == KSP_DIVERGED_ITS )
                || ( tReason == KSP_DIVERGED_DTOL )
                || ( tReason == KSP_DIVERGED_NULL )
                || ( tReason == KSP_DIVERGED_BREAKDOWN )
                || ( tReason == KSP_DIVERGED_BREAKDOWN_BICG ) ;

            // PETSc documents the reason as a property of the collective KSP
            // object, hence uniform across ranks — but neither the soft
            // return below nor the hard abort may bet the run on a
            // documentation claim: one rank returning while its peers enter
            // the collectives further down ( KSPView, collect_lhs, the
            // barrier ) is a deadlock, and under the test/debug THROW
            // reaction even the hard path does not kill the peers. So BOTH
            // verdicts are reduced — any rank failing makes every rank
            // fail, and any rank seeing a hard reason makes the failure
            // hard everywhere. STRUMPACK return codes made the
            // uniformity bet once and lost; two ints per solve is cheap
            int tVerdict[ 2 ] ;
            tVerdict[ 0 ] = ( tReason < 0 || tReasonStatus != 0 ) ? 1 : 0 ;
            tVerdict[ 1 ] = ( tReasonStatus != 0
                || ( tReason < 0 && ! tSoftReason ) ) ? 1 : 0 ;
            if ( this->comm_size() > 1 )
            {
                int tVerdictGlobal[ 2 ] = { 0, 0 } ;
                allreduce( tVerdict, tVerdictGlobal, 2 );
                tVerdict[ 0 ] = tVerdictGlobal[ 0 ] ;
                tVerdict[ 1 ] = tVerdictGlobal[ 1 ] ;
            }

            if ( tVerdict[ 0 ] != 0 )
            {
                if ( this->soft_fail() && tVerdict[ 1 ] == 0 )
                {
                    this->flag_failure() ;

                    if ( this->rank() == 0 )
                    {
                        // This lands INSIDE the controller's iteration box,
                        // so it is drawn as a box line: "   │ " + 69 columns
                        // + "│", matching Controller::print_line(). Written
                        // straight to cout for the same reason the box is —
                        // message() would prepend its own decoration.
                        // rank 0 may be here because a PEER diverged while
                        // its own reason is a converged one — printing that
                        // would be a diagnostic lie
                        string tText = sprint(
                            "PETSc solve failed softly ( %s after %i iterations )",
                            tReason < 0 ?
                                KSPConvergedReasons[ tReason ] :
                                "divergence on a peer rank",
                            ( int ) tNumIterations );

                        // truncate rather than break the frame if a future
                        // reason string is longer than the box
                        if ( tText.size() > 69 )
                        {
                            tText = tText.substr( 0, 69 );
                        }

                        std::cout << "   │ " << tText
                                  << std::string( 69 - tText.size(), ' ' )
                                  << " │" << std::endl ;
                    }

                    // uniform return: the verdict is Allreduced, so every rank
                    // takes this branch and the skipped collectives below
                    // are skipped everywhere. The caller's field distribute
                    // still runs on every rank — the SolverData master
                    // checks failed() before consuming the solution
                    return ;
                }

                BELFEM_ERROR( false,
                    "PETSc KSPSolve did not converge ( %s after %i iterations%s )",
                    tReason < 0 ?
                        KSPConvergedReasons[ tReason ] :
                        "hard divergence on a peer rank",
                    ( int ) tNumIterations,
                    tReasonStatus != 0 ?
                        ", and KSPGetConvergedReason itself failed" : "" );
            }

            if ( this->rank() == 0 )
            {
                // Verbose, not Detailed: this fires once per linear solve, so
                // several times per nonlinear iterate, and drowns the timestep
                // boxes at the level a normal debug run uses. It sits between
                // Detailed ( ordinary diagnostics ) and the KSPView dump gated
                // at Everything below. The hard error on a negative reason
                // above is NOT gated and stays active at every level —
                // silencing the count must not silence the failure.
                //
                // message() prints on every rank, so gate on the master
                message( InfoLevel::Verbose,
                    "    PETSc: %s after %i iterations",
                    KSPConvergedReasons[ tReason ],
                    ( int ) tNumIterations );
            }

            // check if we should print information
            if( gLog.info_level() >= 5 )
            {
                KSPView( mData.mKSP, PETSC_VIEWER_STDOUT_WORLD );
            }

            // recover the solution vector
            if( this->comm_size() > 1 ) // parallel mode
            {
                // extract local portion of the solution vector
                petsctools_get_local_vector( mData.mLHS, mDistMatrix->lhs_vector() );

                mDistMatrix->collect_lhs( aLHS );

#if !defined( NDEBUG ) || defined( DEBUG )
                // Diagnostic: Check for NaN in collected solution
                if( gComm.rank() == 0 )
                {
                    for( index_t i = 0; i < aLHS.length(); ++i )
                    {
                        if( std::isnan( aLHS(i) ) || std::isinf( aLHS(i) ) )
                        {
                            fprintf( stderr, "ERROR: NaN or Inf detected in collected LHS at index %lu (value = %e)\n",
                                   (unsigned long)i, aLHS(i) );
                            fprintf( stderr, "       This occurred after PETSc solve in parallel mode.\n" );
                            BELFEM_ERROR( false, "NaN detected in PETSc solution after collection" );
                        }
                    }
                }
#endif
            }
            else
            {
                petsctools_get_vector( mData.mLHS, mData.mVectorIndices, aLHS );
            }
            comm_barrier() ;
#endif
        }

//------------------------------------------------------------------------------

        void
        PETSC::initialize(
                SpMatrix           & aMatrix,
                const SymmetryMode   aSymmetryMode,
                const int_t          aNumRhsColumns )
        {
#ifdef BELFEM_PETSC
            // call initialize function from parent
            Wrapper::initialize();

            BELFEM_ERROR( mParams->distributed_matrix_type() == DistributedMatrixType::AIJ || gComm.size() < 2,
                          "Need AIJ as matrix type for PETSc %s",
                          to_string( mParams->distributed_matrix_type() ).c_str()
                          );

            // Only check matrix type on rank 0 where aMatrix is valid
            if ( gComm.rank() == 0 )
            {
                BELFEM_ERROR( aMatrix.type() == SpMatrixType::CSR ||
                         aSymmetryMode != SymmetryMode::Unsymmetric,
                         "CSC Matrix must be symmetric if solved with PETSc. Use CSR instead." );
            }

            BELFEM_ERROR(aNumRhsColumns == 1,
                       "Multiple RHS columns not implemented in BELFEM while using PETSc" );


            mNumCols = aNumRhsColumns ;

            if ( this->comm_size() > 1 )
            {
                mDistMatrix = new sparse::PETScAIJ( mParams, &aMatrix );

                // CRITICAL: distribute values BEFORE creating PETSc matrix!
                // The DistMatrixAIJ constructor only sets up structure with zero values.
                // We must fill in actual values before PETSc copies them.
                mDistMatrix->distribute_values( &aMatrix );

                petsctools_create_matrix( mDistMatrix, mData.mMat );

                // get local vector size and global offset
                PetscInt tLocalLength = mDistMatrix->n_rows();
                PetscInt tGlobalOffset = mDistMatrix->dist()[ gComm.rank() ];

                // get global size from distributed matrix (aMatrix is only valid on rank 0!)
                PetscInt tGlobalSize = mDistMatrix->size();

                // create vector indices for parallel mode (global indices)
                this->create_indices( tLocalLength, tGlobalOffset );

                // allocate vectors for parallel mode
                petsctools_allocate_vector(
                    mData.mComm,
                    mData.mLHS,
                    tGlobalSize,
                    tLocalLength );

                petsctools_allocate_vector(
                    mData.mComm,
                    mData.mRHS,
                    tGlobalSize,
                    tLocalLength );
            }
            else
            {
                // get size of matrix
                PetscInt tN = aMatrix.n_rows() ;

                // create the indices
                this->create_indices( tN );

                // allocate vector for right hand side
                petsctools_allocate_vector( mData.mComm, mData.mRHS, tN );

                // allocate vector for left hand side
                petsctools_allocate_vector( mData.mComm, mData.mLHS, tN );

                // link with input matrix
                this->link_matrix( aMatrix );
            }

            // create the preconditioner context and the Krylov subspace solver
            this->create_pc_and_ksp() ;

#endif
        }

//------------------------------------------------------------------------------

        void
        PETSC::free()
        {
#ifdef BELFEM_PETSC
            if( this->is_initialized() )
            {
                if ( mDistMatrix != nullptr )
                {
                    delete mDistMatrix;
                }

                // tidy up
                KSPDestroy( & mData.mKSP );
                MatDestroy( & mData.mMat );

                if( mNumCols == 1 )
                {
                    VecDestroy( & mData.mLHS );
                    VecDestroy( & mData.mRHS );
                }
            }

#endif
            // call function from parent
            Wrapper::free();

        }

//------------------------------------------------------------------------------

        void
        PETSC::create_indices( const PetscInt & aLength, const PetscInt aOffset )
        {
#ifdef BELFEM_PETSC
            mData.mVectorIndices.set_size( aLength );

            std::generate(
                    mData.mVectorIndices.begin(),
                    mData.mVectorIndices.end(),
                    [ n=aOffset-1 ]() mutable { n++; return n; }
                    );
#endif
        }

//------------------------------------------------------------------------------

        void
        PETSC::create_pc_and_ksp()
        {
#ifdef BELFEM_PETSC
            // create the solver
            KSPCreate( mData.mComm, & mData.mKSP );

            // set matrix and preconditioning matrix
            KSPSetOperators( mData.mKSP, mData.mMat, mData.mMat );

            // grab PC pointer
            KSPGetPC( mData.mKSP, & mData.mPC );

            // set epsilon and tolerances. A deck-stated "max iterations"
            // replaces PETSc's default budget of 10000; zero must never
            // reach PETSc — max_it == 0 means ZERO iterations there, not
            // "default". KSPSetFromOptions below still lets a command-line
            // -ksp_max_it override the deck value
            // atol only when the deck states one — PETSc's own default
            // ( 1e-50 ) is already far below any floor, unlike STRUMPACK's
            // library-default absolute tolerance; an unrequested override would
            // change long-verified
            // behavior for no benefit
            PetscErrorCode tStatus = KSPSetTolerances(
                    mData.mKSP,
                    mEpsilon,
                    mParams->have_absolute_tolerance() ?
                        ( PetscReal ) mParams->absolute_tolerance() :
                        PETSC_DEFAULT,
                    PETSC_DEFAULT,
                    mParams->max_iterations() > 0 ?
                        ( PetscInt ) mParams->max_iterations() :
                        PETSC_DEFAULT );
            BELFEM_ERROR( tStatus == 0,
                "KSPSetTolerances failed with error %i", ( int ) tStatus );

            // set the runtime options
            KSPSetFromOptions( mData.mKSP );
#endif
        }

//------------------------------------------------------------------------------

        void
        PETSC::link_matrix( SpMatrix & aMatrix )
        {
#ifdef BELFEM_PETSC
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            BELFEM_ERROR( this->comm_size() < 2, "PETSC::link_matrix() must not be called in parallel mode");

            MatCreateSeqAIJWithArrays(
                    mData.mComm,
                    aMatrix.n_rows(),
                    aMatrix.n_cols(),
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    & mData.mMat );
#endif
        }

//------------------------------------------------------------------------------

        PetscErrorCode
        PETSC::set_preconditioner( const Preconditioner aPreconditioner )
        {
#ifdef BELFEM_PETSC
            PetscErrorCode aStatus = PCSetType( mData.mPC,
                    to_string( aPreconditioner ).c_str() ) ;

            BELFEM_ASSERT( aStatus==0,
                    "PETSc has thrown error %i during PCSetType(): %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

            // Configure GASM (restricted additive Schwarz) if requested
            if( aPreconditioner == Preconditioner::GASM )
            {
                aStatus = PCASMSetType( mData.mPC, PC_ASM_RESTRICT );
                BELFEM_ASSERT( aStatus==0,
                        "PETSc has thrown error %i during PCASMSetType(): %s",
                        ( int ) aStatus,
                        petsctools_error_message( aStatus ).c_str() );
            }

            //  for direct solver, there is never an initial guess
            if( aPreconditioner == Preconditioner::LU )
            {
                this->set_initial_guess_flag( false );
            }
            else
            {
                this->set_initial_guess_flag( true );
            }

            return aStatus ;
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        PetscErrorCode
        PETSC::set_krylovmethod( const KrylovMethod aKrylovMethod )
        {
#ifdef BELFEM_PETSC
            // AUTO is a BELFEM sentinel ( "let the wrapper choose" ), not a PETSc
            // KSP type -- handing it to KSPSetType throws error 86. Resolve it
            // here at the PETSc boundary ( STRUMPACK resolves the same enum on
            // its own side ): a direct LU preconditioner wants KSPPREONLY, any
            // iterative preconditioner wants GMRES ( see KrylovMethod in
            // en_SolverEnums.hpp ).
            KrylovMethod tKrylovMethod = aKrylovMethod ;
            if( tKrylovMethod == KrylovMethod::AUTO )
            {
                tKrylovMethod = mPreconditioner == Preconditioner::LU
                    ? KrylovMethod::PREONLY
                    : KrylovMethod::GMRES ;
            }

            PetscErrorCode aStatus = KSPSetType( mData.mKSP,
                    to_string( tKrylovMethod ).c_str() );

            BELFEM_ASSERT( aStatus==0,
                    "PETSc has thrown error %i during KSPSetType(): %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

            return aStatus ;
#else
            return 0 ;
#endif
        }
//------------------------------------------------------------------------------

        PetscErrorCode
        PETSC::set_initial_guess_flag( const bool aSwitch )
        {
#ifdef BELFEM_PETSC
            PetscErrorCode aStatus ;

            if( aSwitch )
            {
                 aStatus = KSPSetInitialGuessNonzero( mData.mKSP, PETSC_TRUE );
            }
            else
            {
                aStatus = KSPSetInitialGuessNonzero( mData.mKSP, PETSC_FALSE );
            }

            BELFEM_ASSERT( aStatus==0,
                    "PETSc has thrown error %i during KSPSetInitialGuessNonzero(): %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

            return aStatus ;
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        PetscErrorCode
        PETSC::set_matrix_ordering( const ReorderingMethod aReorderingMethod )
        {
#ifdef BELFEM_PETSC
            // Get the factorization preconditioner
            //PC tFactorPC; <- not used
            PetscErrorCode aStatus;

            // Check if we're using a factorization-based preconditioner
            // (LU, ILU, ICC, Cholesky)
            PCType tPCType;
            aStatus = PCGetType( mData.mPC, &tPCType );

            if( aStatus != 0 )
            {
                return aStatus;
            }


            /* Available ordering schemes
             * see https://petsc.org/release/manualpages/MatGraphOperations/MatGetOrdering/
             *  MATORDERINGNATURAL_OR_ND - Nested dissection unless matrix is SBAIJ then it is natural
               MATORDERINGNATURAL - Natural
               MATORDERINGND - Nested Dissection
               MATORDERING1WD - One-way Dissection
               MATORDERINGRCM - Reverse Cuthill-McKee
               MATORDERINGQMD - Quotient Minimum Degree
               MATORDERINGEXTERNAL - Use an ordering internal to the factorzation package and do not compute or use PETSc's
             */
            // Map BELFEM ordering to PETSc ordering
            MatOrderingType tPetscOrdering;
            switch( aReorderingMethod )
            {
                // NATURAL for every handled method is BY DESIGN, not a
                // placeholder: the DofManager has already ordered the free
                // dofs before this matrix exists — symrcm for PETSc and
                // every other non-STRUMPACK solver
                // ( cl_FEM_DofMgr_DofData.cpp, reorder step ) — so
                // "natural" here means KEEP that upstream RCM. Re-running
                // RCM inside PETSc would permute an already-banded matrix
                // to no effect, and a fill-reducing ordering ( ND/METIS )
                // would actively hurt an incomplete factorization: it
                // scatters fill away from the diagonal, which is exactly
                // what ILU(0) then discards. The one thing a PETSc-side
                // ordering could add is a sub-RCM on the OVERLAPPED ASM
                // subdomain graph in parallel ( the overlap rows sit
                // outside the local band ) — expected gain is small and
                // measurable without code via
                // -sub_pc_factor_mat_ordering_type rcm
                case ReorderingMethod::AUTOMATIC:
                case ReorderingMethod::NATURAL:
                case ReorderingMethod::METIS:
                case ReorderingMethod::SCOTCH:
                case ReorderingMethod::PARMETIS:
                case ReorderingMethod::PTSCOTCH:
                {
                    tPetscOrdering = MATORDERINGNATURAL;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid reordering method for PETSc: %s ",
                        to_string( aReorderingMethod ).c_str() );
                    tPetscOrdering = MATORDERINGEXTERNAL;
                }
            }

            // Set the ordering for factorization
            aStatus = PCFactorSetMatOrderingType( mData.mPC, tPetscOrdering );

            BELFEM_ASSERT( aStatus == 0,
                    "PETSc has thrown error %i during PCFactorSetMatOrderingType(): %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

            return aStatus;
#else
            return 0;
#endif
        }

//------------------------------------------------------------------------------
    }
}
