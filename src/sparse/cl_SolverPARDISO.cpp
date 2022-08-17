//
// Created by Christian Messe on 11.07.20.
//

#ifdef BELFEM_PARDISO
#include "pardisotools.hpp"
#endif
#include "cl_Logger.hpp"
#include "cl_SolverPARDISO.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        PARDISO::PARDISO() :
                Wrapper( "PARDISO  " )
        {
            mParameters.set_size( 8, 0 );
            mInfo.set_size( 8, 0 );

        }

//------------------------------------------------------------------------------

        PARDISO::~PARDISO()
        {
            this->free();
        }

//------------------------------------------------------------------------------

        void
        PARDISO::solve(
                SpMatrix & aMatrix,
                Vector <real> & aLHS,
                Vector <real> & aRHS )
        {
#ifdef BELFEM_PARDISO
            BELFEM_ASSERT( this->is_initialized(),
                          "PARDISO has not been initialited" );

            // allocate lhs
            if( aLHS.length() != aRHS.length() )
            {
                aLHS.set_size( aRHS.length(), 0.0 );
            }

            aMatrix.set_indexing_base( mIndexingBase );

            int tStatus = pardisotools_symbolic_factorization(
                aMatrix.n_rows(),
                aMatrix.number_of_nonzeros(),
                1,
                aMatrix.pointers(),
                aMatrix.indices(),
                aMatrix.data() );

            this->check_status( tStatus );

            tStatus = pardisotools_solve(
                    aMatrix.n_rows(),
                    aMatrix.number_of_nonzeros(),
                    1,
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    aLHS.data(),
                    aRHS.data(),
                    mInfo.data() );

            // check for error
            this->check_status( tStatus );

#else
            BELFEM_ERROR( false, "We are not linked against PARDISO" );
#endif
        }

//------------------------------------------------------------------------------

        void
        PARDISO::solve(
                SpMatrix & aMatrix,
                Matrix <real> & aLHS,
                Matrix <real> & aRHS )
        {
#ifdef BELFEM_PARDISO
            BELFEM_ASSERT( this->is_initialized(),
                          "PARDISO has not been initialized" );

            aMatrix.set_indexing_base( mIndexingBase );

            // allocate lhs
            if( aLHS.n_rows() != aRHS.n_rows() ||
                aLHS.n_cols() != aRHS.n_cols()  )
            {
                aLHS.set_size( aRHS.n_rows(), aRHS.n_cols(), 0.0 );
            }

#ifdef BELFEM_ARMADILLO

        int tStatus = pardisotools_solve(
                    aMatrix.n_rows(),
                    aMatrix.number_of_nonzeros(),
                    aRHS.n_cols(),
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    aLHS.data(),
                    aRHS.data(),
                    mInfo.data() );
#else
           // flatten matrix to vector
           Vector< real > tX;
           Vector< real > tY;
           this->mat2vec( aRHS, tY );

           // allocate memory for LHS
           tX.set_size( tY.length() );

           // call pardiso
           int tStatus = pardisotools_solve(
                        aMatrix.n_rows(),
                        aMatrix.number_of_nonzeros(),
                        aRHS.n_cols(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        tX.data(),
                        tY.data(),
                        mInfo.data() );

           // unflatten vector to matrix
           this->vec2mat( tX, aLHS );
#endif
           // check for error
           this->check_status( tStatus );
#else
            BELFEM_ERROR( false, "We are not linked against PARDISO" );
#endif
        }

//------------------------------------------------------------------------------

        real
        PARDISO::get_determinant() const
        {
#ifdef BELFEM_PARDISO
            return mParameters( 7 ) == 0 ?
                    BELFEM_SIGNALING_NAN :
                    pardisotools_get_determinant() ;
#else
            BELFEM_ERROR( false, "We are not linked against PARDISO" );
            return BELFEM_SIGNALING_NAN ;
#endif
        }

//------------------------------------------------------------------------------

        void
        PARDISO::free()
        {
#ifdef BELFEM_PARDISO
            if( this->is_initialized() )
            {
                // tidy up memory
                int tStatus = pardisotools_free();

                if( tStatus != 0 )
                {
                    // create the error message
                    string tMessage = this->error_message( tStatus );

                    // throw error
                    BELFEM_ERROR( tStatus == 0,
                                 "PARDISO has thrown the error: %i during memory freeing:\n%s",
                                 tStatus,
                                 tMessage.c_str() );
                }
            }
#endif
            // call function from parent
            Wrapper::free();
        }

//------------------------------------------------------------------------------

        void
        PARDISO::initialize( SpMatrix & aMatrix,
                             const SymmetryMode aSymmetryMode,
                             const int aNumRhsColumns )
        {
#ifdef BELFEM_PARDISO
            // call initialize function from parent
            Wrapper::initialize();

            aMatrix.set_indexing_base( mIndexingBase );

            if( aSymmetryMode == SymmetryMode::Unsymmetric )
            {
                // set the matrix type
                switch ( aMatrix.type() )
                {
                    case ( SpMatrixType::CSC ):
                    {
                        mParameters( 0 ) = 1;
                        break;
                    }
                    case ( SpMatrixType::CSR ) :
                    {
                        mParameters( 0 ) = 0;
                        break;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Matrix Type not supported by PARDISO." );
                    }
                }
            }
            else
            {
                mParameters( 0 ) = 0;
            }

            // set the indexing base
            mParameters( 1 ) = ( int ) aMatrix.indexing_base() ;

            // symmetry must be general symmetric, others did not work so far
            mParameters( 2 ) = 11;

            // get infolevel from logger and set verbosity
            if( gLog.info_level() >= 5 )
            {
                mParameters( 3 ) = 1;
            }

            // preconditionong
            mParameters( 4 ) = 6 ;

            // max number of refinement steps, 0 is auto
            mParameters( 5 ) = 9 ;

            // reduced odering scheme
            mParameters( 6 ) = 0 ;

            // do not compute determinant
            mParameters( 7 ) = 0 ;

            // create the parameter list
            int tStatus = pardisotools_initialize_parameters( mParameters.data() );

            // check for error
            if( tStatus != 0 )
            {
                // create the error message
                string tMessage = this->error_message( tStatus );

                // throw error
                BELFEM_ERROR( tStatus == 0,
                             "PARDISO has thrown the error: %i during parameter initialization:\n%s",
                             tStatus,
                             tMessage.c_str() );
            }

            // perform the symbolic factorization
            aMatrix.set_indexing_base( mIndexingBase );

            tStatus = pardisotools_symbolic_factorization(
                aMatrix.n_rows(),
                aMatrix.number_of_nonzeros(),
                aNumRhsColumns,
                aMatrix.pointers(),
                aMatrix.indices(),
                aMatrix.data() );

            if( tStatus != 0 )
            {
                // create the error message
                string tMessage = this->error_message( tStatus );

                // throw error
                BELFEM_ERROR( tStatus == 0,
                             "PARDISO has thrown the error: %i during symbolic factorization:\n%s",
                             tStatus,
                             tMessage.c_str() );
            }
#endif
        }


//-----------------------------------------------------------------------------

        void
        PARDISO::check_status( const int aStatus )
        {
            if( aStatus != 0 )
            {
                // create the error message
                string tMessage = this->error_message( aStatus );

                string tPhaseString = "" ;

                // grab phase
                if ( mInfo( 0 ) == 22 )
                {
                    tPhaseString = "during numeric factorization, phase " ;
                }
                else if ( mInfo( 0 ) == 33 )
                {
                    tPhaseString = "during numeric solving, phase " ;
                }
                else
                {
                    tPhaseString = "during phase " ;
                }

                // throw error
                BELFEM_ERROR( aStatus == 0,
                             "PARDISO has thrown the error: %i %s %i:\n%s",
                             aStatus,
                             tPhaseString.c_str(),
                             mInfo( 0 ),
                             tMessage.c_str() );
            }
        }

//-----------------------------------------------------------------------------

        string
        PARDISO::error_message( const int aStatus ) const
        {
            switch( aStatus )
            {
                case(  0 ) :
                {
                    return "No error";
                }
                case( -1 ) :
                {
                    return "Input inconsistent.";
                }
                case( -2 ) :
                {
                    return "Not enough memory.";
                }
                case( -3 ) :
                {
                    return "Reordering problem.";
                }
                case( -4 ) :
                {
                    return "Zero pivot, numerical fact. or iterative refinement problem.";
                }
                case( -5 ) :
                {
                    return "Unclassified (internal) error.";
                }
                case( -6 ) :
                {
                    return "Preordering failed.";
                }
                case( -7 ) :
                {
                    return "Diagonal matrix problem.";
                }
                case( -8 ) :
                {
                    return "32-bit integer overflow problem.";
                }
                case( -10 ) :
                {
                    return "No license file found.";
                }
                case( -11 ) :
                {
                    return "License has expired";
                }
                case( -12 ) :
                {
                    return "Wrong username or hostname";
                }
                case( -100 ) :
                {
                    return "Reached maximum number of Krylov-subspace iteration in iterative solver.";
                }
                case( -101 ) :
                {
                    return "No sufficient convergence in Krylov-subspace iteration within 25 iterations.";
                }
                case( -102 ) :
                {
                    return "Error in Krylov-subspace iteration";
                }
                case( -103 ) :
                {
                    return "Break-Down in Krylov-subspace iteration.";
                }
                default:
                {
                    return "Unknown Error";
                }
            }
        }


    }
}