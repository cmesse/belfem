//
// Created by Christian Messe on 11.07.20.
//
#ifdef BELFEM_PETSC
#include <algorithm>   // for create_indices
#endif
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_SolverPETSC.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        PETSC::PETSC() :
                Wrapper( "PETSC   " )
        {
#ifdef BELFEM_PETSC
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

            // enforce zero based indexing on matrix
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            Vector< index_t > tLengths( 2 );

            if( gComm.rank() == 0 )
            {
                tLengths( 0 ) = aRHS.length() ;
                tLengths( 1 ) = aLHS.length() ;

                if( gComm.size() > 1 )
                {
                    send_same( mDistributor->comm_list(), tLengths );
                }
            }
            else
            {
                receive( 0, tLengths );
            }
            index_t tRhsLength = tLengths( 0 );
            index_t tLhsLength = tLengths( 1 );

            // bool set flag if a nonzero LHS is given
            bool tHaveLHS = tLengths( 0 ) == tLengths( 1 );

            this->set_preconditioner( mPreconditioner );
            this->set_krylovmethod( mKrylovMethod );

            // check if LHS exists, use direct solver if not
            if( ! tHaveLHS )
            {
                // allocate LHS
                aLHS.set_size( tRhsLength, 0.0 );

                /*// use direct solver
                if( this->comm_size() == 0 )
                {
                    this->set_preconditioner( Preconditioner::LU );
                }
                else
                {
                    // FIXME: add something here
                }
                this->set_krylovmethod( KrylovMethod::PREONLY ); */

                // do not have initial guess
                this->set_initial_guess_flag( false );
            }
            else
            {
                // use input as initial guess
                this->set_initial_guess_flag( true );
            }

            // on non-master, also allocate RHS
            if( tRhsLength != aRHS.length() )
            {
                aRHS.set_size( tRhsLength );
            }

            // update data containers
            if( this->comm_size() > 1 ) // parallel mode
            {
                // synchronize data and set vectors
                mDistributor->synchronize( aLHS, aRHS, tHaveLHS );
            }
            else  // sequential mode
            {
                // use LHS as initial guess
                if( tLhsLength == tRhsLength )
                {
                    petsctools_set_vector( aLHS, mData.mVectorIndices, mData.mLHS );
                }

                // update RHS data
                petsctools_set_vector( aRHS, mData.mVectorIndices, mData.mRHS );
            }

            // MatView( mData.mMat, PETSC_VIEWER_STDOUT_( mData.mComm) );
            //VecView( mData.mLHS,  PETSC_VIEWER_STDOUT_(  mData.mComm ) );

             // solve the system
            PetscErrorCode tStatus = KSPSolve( mData.mKSP, mData.mRHS, mData.mLHS );

            // VecView( mData.mRHS,  PETSC_VIEWER_STDOUT_(  mData.mComm ) );

            // check error
            if( tStatus != 0 )
            {
                string tMessage = petsctools_error_message( tStatus );

                BELFEM_ERROR( false,
                       "PETSC has thrown the error %i in function KSPSolve()\n%s",
                        tStatus,
                        tMessage.c_str() );
            }

            // check if we should print information
            if( gLog.info_level() >= 5 )
            {
                KSPView( mData.mKSP, PETSC_VIEWER_STDOUT_WORLD );
            }

            // recover the solution vector
            if( this->comm_size() > 1 ) // parallel mode
            {
                mDistributor->recover_lhs( aLHS );
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
                const int            aNumRhsColumns )
        {
#ifdef BELFEM_PETSC
            // call initialize function from parent
            Wrapper::initialize();

            // sanity check
            if ( this->rank() == 0 )
            {
                // make sure that matrix has the right format
                BELFEM_ERROR( aMatrix.type() == SpMatrixType::CSR ||
                    aSymmetryMode != SymmetryMode::Unsymmetric,
                    "CSC Matrix must be symmetric if solved with PETSc. Use CSR instead." );


                 //
                BELFEM_ERROR(aNumRhsColumns == 1,
                        "Multiple RHS columns not implemented in BELFEM while using PETSC" );
            }

            // copy number of columns
            mNumCols = aNumRhsColumns ;

            if( this->comm_size() > 1 ) // parallel mode
            {
                // create the distributor
                mDistributor = new PetscDistributor( mData, aMatrix );
            }
            else // sequential mode
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
                if( this->comm_size() > 1 )
                {
                    delete mDistributor ;
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
        PETSC::create_indices( const PetscInt & aLength  )
        {
#ifdef BELFEM_PETSC
            mData.mVectorIndices.set_size( aLength );

            std::generate(
                    mData.mVectorIndices.begin(),
                    mData.mVectorIndices.end(),
                    [ n=-1 ]() mutable { n++; return n; }
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

            // set epsilon and default tolerances
            KSPSetTolerances(
                    mData.mKSP,
                    mEpsilon,
                    PETSC_DEFAULT,
                    PETSC_DEFAULT,
                    PETSC_DEFAULT);

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

            if ( this->comm_size() > 1 ) // parallel mode
            {
                BELFEM_ERROR( false, "Parellel mode is not implemented yet");
            }
            else
            {
                MatCreateSeqAIJWithArrays(
                        mData.mComm,
                        aMatrix.n_rows(),
                        aMatrix.n_cols(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        & mData.mMat );
            }
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

            // only for direct solver, there is no initial guess
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
            PetscErrorCode aStatus = KSPSetType( mData.mKSP,
                    to_string( aKrylovMethod ).c_str() );

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
    }
}