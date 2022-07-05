//
// Created by Christian Messe on 10.07.20.
//


#ifdef BELFEM_SUITESPARSE
#include <umfpack.h>
#endif
#include "cl_SolverUMFPACK.hpp"
namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        UMFPACK::UMFPACK() :
            Wrapper( "UMFPACK " )
        {
        }

//------------------------------------------------------------------------------

        UMFPACK::~UMFPACK()
        {
            this->free();
        }

//------------------------------------------------------------------------------

        void
        UMFPACK::initialize( SpMatrix & aMatrix,
            const SymmetryMode aSymmetryMode,
            const int aNumRhsColumns )
        {
#ifdef BELFEM_SUITESPARSE
            // call initialize function from parent
            Wrapper::initialize();

            // check matrix type and select setting
            mTransposedFlag = 0;

            switch ( aMatrix.type() )
            {
                case( SpMatrixType::CSC ) :
                {
                    mTransposedFlag = UMFPACK_A;
                    break;
                }
                case( SpMatrixType::CSR ) :
                {
                    // if the system is not symmetric, we must tell
                    // UMFPACK to transpose the matrix.
                    // Otherwise, it does not matter.
                    if( aSymmetryMode == SymmetryMode::Unsymmetric )
                    {
                        mTransposedFlag = UMFPACK_At;
                    }
                    else
                    {
                        mTransposedFlag = UMFPACK_A;
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid matrix type for UMFPACK.");
                    break;
                }
            }

            // make sure that matrix is stored zero-based
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            // create a null pointer
            double *null = ( double * ) nullptr;

            // create symbolic factorization
            int tStatus = umfpack_di_symbolic (
                    aMatrix.n_rows(),
                    aMatrix.n_cols(),
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    null,
                    &mSymbolic,
                    null,
                    null );

            // check for error
            if( tStatus != 0 )
            {
                // create error message
                string tMessage = this->error_message( tStatus );

                // throw error
                BELFEM_ERROR( false,
                    "UMFPACK has thrown the error: %i at umfpack_di_symbolic():\n%s",
                    tStatus,
                    tMessage.c_str() );
            }

#endif
        }

//------------------------------------------------------------------------------

        void
        UMFPACK::free()
        {
#ifdef BELFEM_SUITESPARSE
            if( this->is_initialized() )
            {
                umfpack_di_free_symbolic ( &mSymbolic );
                mSymbolic = nullptr ;
            }
#endif
            // call fre function from parent
            Wrapper::free();
        }

//------------------------------------------------------------------------------

        void
        UMFPACK::solve(
                SpMatrix       & aMatrix,
                Vector< real > & aLHS,
                Vector< real > & aRHS )
        {
#ifdef BELFEM_SUITESPARSE
            // make sure that the wrapper has been initialized
            BELFEM_ASSERT( this->is_initialized(),
                "UMFPACK has not been initialited" );

            // allocate space for LHS
            if( aLHS.length() != aRHS.length() )
            {
                aLHS.set_size( aRHS.length() );
            }

            // make sure that matrix is stored zero-based
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            // create a null pointer
            double *null = ( double * ) nullptr;

            // numeric factorization
            void * tNumeric  = nullptr ;

            // From the symbolic factorization information, carry out the numeric factorization.
            int tStatus = umfpack_di_numeric (
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    mSymbolic,
                    &tNumeric,
                    null,
                    null );

            // check for error
            if( tStatus != 0 )
            {
                // create error message
                string tMessage = this->error_message( tStatus );

                // throw error
                BELFEM_ERROR( false,
                             "UMFPACK has thrown the error: %i at  umfpack_di_numeric():\n%s",
                             tStatus,
                             tMessage.c_str() );
            }

            // Using the numeric factorization, solve the linear system.
            tStatus = umfpack_di_solve (
                    mTransposedFlag,
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    aLHS.data(),
                    aRHS.data(),
                    tNumeric,
                    null,
                    null );

            // Free the numeric factorization.
            umfpack_di_free_numeric ( &tNumeric );

            // check for error
            if( tStatus != 0 )
            {
                // create error message
                string tMessage = this->error_message( tStatus );

                // throw error
                BELFEM_ERROR( tStatus == 0,
                             "UMFPACK has thrown the error: %i at  umfpack_di_solve():\n%s",
                             tStatus,
                             tMessage.c_str() );
            }
#else
            BELFEM_ERROR( false, "We are not linked against UMFPACK." );
#endif
        }

//------------------------------------------------------------------------------

        void
        UMFPACK::solve(
                SpMatrix & aMatrix,
                Matrix< real > & aLHS,
                Matrix< real > & aRHS )
        {
#ifdef BELFEM_SUITESPARSE
            // make sure that the wrapper has been initialized
            BELFEM_ASSERT( this->is_initialized(),
                          "UMFPACK has not been initialited" );

            // allocate space for LHS
            if( aLHS.n_rows() != aRHS.n_rows() ||
                aLHS.n_cols() != aRHS.n_cols() )
            {
                aLHS.set_size( aRHS.n_rows(), aRHS.n_cols() );
            }

            // make sure that matrix is stored zero-based
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            // create a null pointer
            double *null = ( double * ) nullptr;

            // numeric factorization
            void * tNumeric  = nullptr ;

            // From the symbolic factorization information, carry out the numeric factorization.
            int tStatus = umfpack_di_numeric (
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    mSymbolic,
                    &tNumeric,
                    null,
                    null );

            // check for error
            if( tStatus != 0 )
            {
                // create error message
                string tMessage = this->error_message( tStatus );

                // throw error
                BELFEM_ERROR( false,
                             "UMFPACK has thrown the error: %i at  umfpack_di_numeric():\n%s",
                             tStatus,
                             tMessage.c_str() );
            }

            // temporary vector for LHS
            Vector< real > tX( aRHS.n_rows() );

            // temporary vector for RHS
            Vector< real > tY( aRHS.n_rows() );

            // get number of cols
            uint tNcols = aRHS.n_cols() ;

            // loop over all columns
            for( uint k=0; k<tNcols; ++k )
            {
                // poulate RHS
                tY = aRHS.col( k );

                // Using the numeric factorization, solve the linear system.
                tStatus = umfpack_di_solve(
                        mTransposedFlag,
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        tX.data(),
                        tY.data(),
                        tNumeric,
                        null,
                        null );

                // check for error
                if ( tStatus != 0 )
                {
                    std::string tMessage = this->error_message( tStatus );

                    BELFEM_ERROR( false,
                                 "UMFPACK has thrown the error: %i at umfpack_di_solve(): %i\n%s",
                                 tStatus,
                                 tMessage.c_str() );
                }

                // now we copy the right hand side into the output
                aLHS.set_col( k, tX );
            }

            // Free the numeric factorization.
            umfpack_di_free_numeric ( &tNumeric );
#else
            BELFEM_ERROR( false, "We are not linked against UMFPACK." );
#endif
        }

//------------------------------------------------------------------------------

        string
        UMFPACK::error_message( const int aStatus ) const
        {
            switch( aStatus )
            {
                case ( 0 ) :
                {
                    return "No error";
                }
                case ( 1 ) :
                {
                    return "Matrix is singular.";
                }
                case ( 2 ) :
                {
                    return "The determinant is nonzero, but smaller in magnitude than the smallest positive floating-point number.";
                }
                case ( 3 ) :
                {
                    return "The determinant is larger in magnitude than the largest positive floating-point number.";
                }
                case ( -3 ) :
                {
                    return "Memory leak at numeric object?";
                }
                case ( -4 ) :
                {
                    return "Memory leak at symbolic object?";
                }
                case ( -5 ) :
                {
                    return "Passed nullptr but need argument.";
                }
                case ( -6 ) :
                {
                    return "Matrix size must be > 0.";
                }
                case ( -8 ) :
                {
                    return "The matrix is invalid.";
                }
                case ( -11 ) :
                {
                    return "Wrong pattern.";
                }
                case ( -13 ) :
                {
                    return "The sys argument provided to one of the solve routines is invalid.";
                }
                case ( -15 ) :
                {
                    return "The permutation vector provided as input is invalid.";
                }
                case ( -17 ) :
                {
                    return "Error in reading and writing objects from file";
                }
                case ( -18 ) :
                {
                    return "The ordering method failed.";
                }
                case ( -911 ) :
                {
                    return "An internal error in UMFPACK has occurred, of unknown cause.";
                }
                default:
                {
                    return "Unknown Error";
                }
            }
        }

//------------------------------------------------------------------------------
    }
}