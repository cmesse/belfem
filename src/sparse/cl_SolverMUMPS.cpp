//
// Created by Christian Messe on 10.07.20.
//

#include "cl_SolverMUMPS.hpp"

#ifdef BELFEM_MUMPS
#include "mumpstools.hpp"
#endif

#include "cl_Logger.hpp"
#include "commtools.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        MUMPS::MUMPS( const proc_t aMasterRank ) :
#ifdef BELFEM_MUMPS
                Wrapper( "MUMPS    " ),
                mMasterRank( aMasterRank )
#else
                Wrapper( "MUMPS   " )
#endif
        {
            // allocate user settings
            mIParameters.set_size( 10, 0 );
            mRParameters.set_size( 1, 0.0 );

            // allocate information vector
            mInfo.set_size( 40, 0 );

            this->init_defaults();
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
            // rank of the Host
            mIParameters( 0 ) = ( int ) mMasterRank ;

            // Host is working
            mIParameters( 1 ) = 1 ;

            // Serial Reordering      : 0 - AMD
            //                        : 1 - user pivot
            //                        : 2 - AMF
            //                        : 3 - SCOTCH
            //                        : 4 - PORD
            //                        : 5 - METIS
            //                        : 6 - QAMD
            //                        : 7 - AUTOMATIC
            mIParameters( 4 ) = 7 ;

            //    Parallel Reordering : 0 - AUTOMATIC
            //                        : 1 - PT-SCOTCH
            //                        : 2 - PARMETIS
            mIParameters( 5 ) = 0 ;

            //    Determinant         : 0 - no
            //                        : 1 - yes
            mIParameters( 7 ) = 0 ;

            // memory relaxation, 0: default
            mIParameters( 8 ) = 0 ;

            // blr setting
            mIParameters( 9 )  = 0 ;
            mRParameters( 0 )  = 0.0 ;
#endif
        }

//------------------------------------------------------------------------------

        void
        MUMPS::initialize( SpMatrix & aMatrix,
                const SymmetryMode    aSymmetryMode,
                const int             aNumRhsColumns )
        {
#ifdef BELFEM_MUMPS
            // call initialize function from parent
            Wrapper::initialize();


            // symmetry mode
            mIParameters( 2 ) = 0; // sould be ( int ) aSymmetryMode,
                                  // but didn't work, todo: find out why

            // info level
            mIParameters( 3 ) = ( int ) gLog.info_level() ;


            // 5: Refinement          :<0 - fixed number of steps
            //                        : 0 - none
            //                        :>0 - maximum number of refinement steps
            mIParameters( 6 ) = aNumRhsColumns == 1 ? 9 : 0 ;
#endif
        }

//------------------------------------------------------------------------------

        real
        MUMPS::get_determinant() const
        {
#ifdef BELFEM_MUMPS
            return mIParameters( 7 ) == 0 ?
                    BELFEM_SIGNALING_NAN :
                    mumpstools_get_determinant() ;
#else
            BELFEM_ERROR( false, "We are not linked against MUMPS" );
            return BELFEM_SIGNALING_NAN ;
#endif
        }

//------------------------------------------------------------------------------

        void
        MUMPS::free()
        {
            // call function from parent, nothing else to be done here
            Wrapper::free();
        }

//------------------------------------------------------------------------------

        void MUMPS::solve( SpMatrix & aMatrix,
                           Vector< real > & aLHS,
                           Vector< real > & aRHS )
        {
#ifdef BELFEM_MUMPS

            // make sure that the wrapper has been initialited
            BELFEM_ASSERT( this->is_initialized(),
                "MUMPS has not been initialited" );

            // reset info vector
            mInfo.fill( 0 );

            if( aLHS.length() != aRHS.length() )
            {
                aLHS.set_size( aRHS.length(), 0.0 );
            }

            if( this->rank() == mMasterRank )
            {
                // make sure that all indices have been created
                aMatrix.create_coo_indices() ;

                // make sure that matrix is stored one-based
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );

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
                        mInfo.data() );

                /*call_mumps(
                        aMatrix.n_rows(),
                        aMatrix.number_of_nonzeros(),
                        1,
                        aMatrix.rows(),
                        aMatrix.cols(),
                        aMatrix.data(),
                        aLHS.data(),
                        aRHS.data(),
                        5,
                        0,
                        mInfo.data() ); */

            }
            else
            {
                // solve the system as servant
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
                        mInfo.data() );
            }

            // check result
            if( mInfo( 0 ) != 0 )
            {
                std::string tMessage = this->error_message(
                        mInfo.data(),
                        aMatrix.n_rows(),
                        aMatrix.number_of_nonzeros() );

                BELFEM_ERROR( false,
                       "MUMPS has thrown the error: %i\n%s",
                        mInfo( 0 ),
                        tMessage.c_str() );
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

            // reset info vector
            mInfo.fill( 0 );

            // allocate matrix sizes of not set
            if(  aLHS.n_rows() != aRHS.n_rows() || aLHS.n_cols() != aRHS.n_cols() )
            {
                aLHS.set_size( aRHS.n_rows(), aRHS.n_cols(), 0.0 );
            }

            // make sure that matrix is stored one-based
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );

            // make sure that all indices have been created
            aMatrix.create_coo_indices() ;

            if( this->rank() == mMasterRank )
            {
                // make sure that all indices have been created
                aMatrix.create_coo_indices() ;

                // make sure that matrix is stored one-based
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );
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
                        mInfo.data() );
#else
               Vector< real > tX;
               Vector< real > tY;
               this->mat2vec( aRHS, tY );

               // allocate memory for LHS
               tX.set_size( tY.length() );

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
                        mInfo.data() );

               // unflatten vector to matrix
               this->vec2mat( tX, aLHS );
#endif
            }
            else
            {
                // solve the system as servant
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
                        mInfo.data() );
            }

            // check result
            if( mInfo( 0 ) != 0 )
            {
                std::string tMessage = this->error_message(
                        mInfo.data(),
                        aMatrix.n_rows(),
                        aMatrix.number_of_nonzeros() );

                BELFEM_ERROR( false,
                       "MUMPS has thrown the error: %i\n%s",
                        mInfo( 0 ),
                        tMessage.c_str() );
            }
#else
            BELFEM_ERROR( false, "We are not linked against MUMPS" );
#endif

        }

//------------------------------------------------------------------------------

        string
        MUMPS::error_message(
                const int   * aInfo,
                const int   & aN,
                const int   & aNNZ )
        {
            // the string containing the error message
            string aMessage;

            switch( aInfo[ 0 ] )
            {
                case( 0 ) :
                {
                    aMessage = "No error";
                    break;
                }
                case( -1 ) :
                {
                    aMessage = sprint( "An error occurred on proc %i.", aInfo[ 1 ] );
                    break;
                }
                case( -2 ) :
                {
                    if( aInfo[ 3 ] == aN )
                    {
                        aMessage = "Matrix Size out of range.";
                    }
                    else if( aInfo[ 3 ] == aNNZ )
                    {
                        aMessage = "Number of nonzeros out of range.";
                    }
                    else
                    {
                        aMessage = "unknown";
                    }
                    break;
                }
                case( -3 ) :
                {
                    aMessage = "Invalid JOB index passed to MUMPS.";
                    break;
                }
                case( -4 ) :
                {
                    aMessage =  sprint(
                            "Error in user-provided permutation array PERM IN at position %i.",
                            aInfo[ 1 ] );
                    break;
                }
                case( -5 ) :
                {
                    aMessage =  sprint(
                            "Problem of real workspace allocation of size %i during analysis.",
                            aInfo[ 1 ] );
                    break;
                }
                case( -6 ) :
                {
                    aMessage =  sprint(
                            "Matrix is singular in structure %i.",
                            aInfo[ 1 ] );
                    break;
                }
                case( -7 ) :
                {
                    aMessage =  sprint(
                            "Problem of integer workspace allocation of size %i.",
                            aInfo[ 1 ] );
                    break;
                }
                case( -8 ) :
                {
                    aMessage = "Main internal integer workarray IS too small for factorization.";
                    break;
                }
                case( -9 ) :
                {
                    int tVal;
                    if( aInfo[ 1 ] < 0 )
                    {
                        tVal = -aInfo[ 1 ]*1e6;
                    }
                    else
                    {
                        tVal = aInfo[ 2 ];
                    }

                    aMessage =  sprint(
                            "Main workarray of size %i is to small.",
                            tVal );

                    break;
                }
                case( -10 ) :
                {
                    aMessage = "Matrix is singular.";
                    break;
                }
                case( -11 ) :
                {
                    if( aInfo[ 1 ] > 0 )
                    {
                        aMessage = sprint(
                                "Internal real/complex workarray S or LWK_USER is by %i too small for solution.",
                                aInfo[ 1 ] );
                    }
                    else
                    {
                        aMessage =
                                "Internal real/complex workarray S or LWK_USER is by %i too small for solution.";
                    }
                    break;
                }
                case( -12 ) :
                {
                    aMessage = "Internal real/complex workarray S too small for iterative refinement.";
                    break;
                }
                case( -13 ) :
                {
                    int tVal;
                    if( aInfo[ 1 ] < 0 )
                    {
                        tVal = -aInfo[ 1 ]*1e6;
                    }
                    else
                    {
                        tVal = aInfo[ 1 ];
                    }

                    aMessage =  sprint(
                            "Problem of workspace allocation of size %i.",
                            tVal );

                    break;
                }
                case( -14 ) :
                {
                    aMessage = "Internal integer workarray IS too small for solution.";
                    break;
                }
                case( -15 ) :
                {
                    aMessage = "Integer workarray IS too small for iterative refinement and/or error analysis.";
                    break;
                }
                case( -16 ) :
                {
                    aMessage = sprint(
                            "N=%i is out of range",
                            aInfo[ 1 ] );
                    break;
                }
                case( -17 ) :
                {
                    aMessage = "The internal send buffer that was allocated dynamically by MUMPS on the processor is too small.";
                    break;
                }
                case( -20 ) :
                {
                    aMessage = sprint( "The internal reception buffer that was allocated dynamically by MUMPS is too small. Need %i.",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -21 ) :
                {
                    aMessage = "Value of PAR=0 is not allowed because only one processor is available";
                    break;
                }
                case( -22 ) :
                {
                    std::string tLabel;

                    switch( aInfo[ 1 ] )
                    {
                        case(  1 ) :
                        {
                            tLabel = "IRN";
                            break;
                        }
                        case(  2 ) :
                        {
                            tLabel = "JCN";
                            break;
                        }
                        case(  3 ) :
                        {
                            tLabel = "R\n"
                                     "PERM_IN";
                            break;
                        }
                        case(  4 ) :
                        {
                            tLabel = "A";
                            break;
                        }
                        case(  5 ) :
                        {
                            tLabel = "ROWSCA";
                            break;
                        }
                        case(  6 ) :
                        {
                            tLabel = "COLSCA";
                            break;
                        }
                        case(  7 ) :
                        {
                            tLabel = "RHS";
                            break;
                        }
                        case(  8 ) :
                        {
                            tLabel = "LISTVAR_SHUR";
                            break;
                        }
                        case(  9 ) :
                        {
                            tLabel = "SHUR";
                            break;
                        }
                        case( 10 ) :
                        {
                            tLabel = "RHS_SPARSE";
                            break;
                        }
                        case( 11 ) :
                        {
                            tLabel = "IRHS_SPARSE";
                            break;
                        }
                        case( 12 ) :
                        {
                            tLabel = "IRHS_PTR";
                            break;
                        }
                        case( 13 ) :
                        {
                            tLabel = "ISOL_loc";
                            break;
                        }
                        case( 14 ) :
                        {
                            tLabel = "SOL_loc";
                            break;
                        }
                        case( 15 ) :
                        {
                            tLabel = "REDRHS";
                            break;
                        }
                        case( 16 ) :
                        {
                            tLabel = "IRN_loc/JCN_loc";
                            break;
                        }
                        default:
                        {
                            tLabel = "unknown";
                            break;
                        }
                    }

                    aMessage = sprint( "Pointer array %s not associated, has unsufficent size or was associated but should not be.",
                                       tLabel.c_str() );

                    break;
                }
                case( -23 ) :
                {
                    aMessage = "MPI was not initialized.";
                    break;
                }
                case( -24 ) :
                {
                    aMessage = sprint( "NELT=%i is out of range.", aInfo[ 1 ] );
                    break;
                }
                case( -25 ) :
                {
                    aMessage = "A problem has occurred in the initialization of the BLACS";
                    break;
                }
                case( -26 ) :
                {
                    aMessage = sprint( "LRHS=%i is out of range-",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -27 ) :
                {
                    aMessage = sprint( "NZ_RHS and IRHS_PTR(NZ_RHS+1)=%i do not match.",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -28 ) :
                {
                    aMessage = sprint( "IRHS_PTR=%i is not equal to 1.",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -29 ) :
                {
                    aMessage = sprint( "LSOL_loc=%i is smaller than %i.",
                                       aInfo[ 1 ],
                                       aInfo[ 22 ] );
                    break;
                }
                case( -30 ) :
                {
                    aMessage = sprint( "SHUR_LLD=%i is out of range.",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -31 ) :
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
                case( -32 ) :
                {
                    aMessage = sprint( "Value NHRS=%i not compatible with user defined setting.",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -33 ) :
                {
                    aMessage = "ICNTL(26) was asked for during solve phase (or during the factorization – see ICNTL(32) ) but the Schur complement was not asked for at the analysis phase (ICNTL(19)).";
                    break;
                }
                case( -34 ) :
                {
                    aMessage = sprint( "LREDRHS=%i is out of range.",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -36 ) :
                {
                    aMessage = sprint( "Incompativle values of ICNTL(26)=%i and INFOG(28).",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -37 ) :
                {
                    aMessage = sprint( "Value of ICNTL(25) incompatible with some other parameter. ( INFO(2)=%i)",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -38 ) :
                {
                    aMessage = "Parallel analysis was set but neither PT-SCOTCH or ParMetis were provided.";
                    break;
                }
                case( -39 ) :
                {
                    aMessage = "Incompatible values for ICNTL(28) and ICNTL(5) and/or ICNTL(19) and/or ICNTL(6).";
                    break;
                }
                case( -40 ) :
                {
                    aMessage = "The matrix was indicated to be positive definite but is not";
                    break;
                }
                case( -41 ) :
                {
                    aMessage = "Incompatible value of LWK USER between factorization and solution phases.";
                    break;
                }
                case( -42 ) :
                {
                    aMessage = sprint( "Forward during factorization was set but but the value of NRHS=%i on the host is incorrect",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -43 ) :
                {
                    aMessage = sprint( " Incompatible values of ICNTL(32) and ICNTL(%i).",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -44 ) :
                {
                    aMessage = sprint( "The solve phase (JOB=3) cannot be performed because the factors or part of the factors are not\n"
                                       "available ( ICNTL(32)=%i)", aInfo[ 1 ] );
                    break;
                }
                case( -45 ) :
                {
                    aMessage = sprint( "RHS=%i must be >0.", aInfo[ 1 ] );
                    break;
                }
                case( -46 ) :
                {
                    aMessage = sprint( "NZ_RHS=%i must be >0.", aInfo[ 1 ] );
                    break;
                }
                case( -47 ) :
                {
                    aMessage = sprint( "Entries of A^(−1) were requested during the solve phase but the constraint NRS=%i is not resprected",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -48 ) :
                {
                    aMessage = sprint( "A^(−1) Incompatible values if ICNTL(30) and ICNTL(%i).",
                                       aInfo[ 1 ] );
                    break;
                }
                case( -49 ) :
                {
                    aMessage = sprint( "Incorrect value: SIZE_SCHUR=%i", aInfo[ 1 ] );
                    break;
                }
                case( -50 ) :
                {
                    aMessage = " An error occurred while computing the fill-reducing ordering during the analysis phase.";
                    break;
                }
                case( -51 ) :
                {
                    int tVal;
                    if( aInfo[ 1 ] < 0 )
                    {
                        tVal = -aInfo[ 1 ]*1e6;
                    }
                    else
                    {
                        tVal = aInfo[ 1 ];
                    }

                    aMessage = sprint( "The linked ordering library cannot hadle graphs of size %i ( must be < 2^31-1).",
                                       tVal );
                    break;
                }
                case( -52 ) :
                {
                    aMessage = "The linked ordering library must have 64-bit default integers";
                    break;
                }
                case( -53 ) :
                {
                    aMessage = "Internal error that could be due to inconsistent input data between two consecutive calls.";
                    break;
                }
                default:
                {
                    aMessage = "unknown error";
                    break;
                }

            }

            return aMessage;
        }

//------------------------------------------------------------------------------

        void
        MUMPS::set_reordering(
                const SerialReodrdering   aSerial,
                const ParallelReodrdering aParallel
        )
        {
            mIParameters( 4 ) = static_cast< int >( aSerial );
            mIParameters( 5 ) = static_cast< int >( aParallel );
        }

//------------------------------------------------------------------------------

        void
        MUMPS::set_block_low_ranking(
                const BlockLowRanking aBLK,
                const real            aEpsilon
        )
        {
            mIParameters( 9 ) = static_cast< int > ( aBLK );
            mRParameters( 0 ) = aEpsilon ;
        }

//------------------------------------------------------------------------------
    }
}
