//
// Created by christian on 8/16/22.
//
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "cl_SolverSTRUMPACK.hpp"
#include "assert.hpp"
// Externally Defined Global Communicator
extern belfem::Logger gLog;
extern belfem::Communicator gComm;
namespace belfem
{
    namespace solver
    {
//----------------------------------------------------------------

        STRUMPACK::STRUMPACK() :
                Wrapper( "STRUMPACK" )
        {

        }
//----------------------------------------------------------------

        STRUMPACK::~STRUMPACK()
        {
#ifdef BELFEM_STRUMPACK
            if( mSolver != nullptr )
            {
                delete mSolver;
            }
#endif
            if( mArgV != nullptr )
            {
                delete mArgV ;
            }

        }

//------------------------------------------------------------------------------

        void
        STRUMPACK::initialize( SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode,
                    const int aNumRhsColumns )
        {
            BELFEM_ASSERT(  aNumRhsColumns == 1,
                            "number of RHS cols must be 1 for STRUMPACK" );
#ifdef BELFEM_STRUMPACK
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            mArgC = gComm.arguments().size() ;

            mArgV = new StringList( mArgC );

            for( string tString : gComm.arguments( ) )
            {
                mArgV->push( tString );
            }

            // see https://portal.nersc.gov/project/sparse/strumpack/v6.3.1/sparse_example_usage.html#autotoc_md9
            BELFEM_ERROR( gComm.size() == 1 , "Parallel STRUMPACK not supportet yet");

            // create the solver
            mSolver = new strumpack::StrumpackSparseSolver< double, int >(
                    mArgC,
                    mArgV->data(),
                    gLog.info_level() >= 5 );

            // set the options
            mSolver->options().set_from_command_line( mArgC, mArgV->data() );

            mSolver->set_csr_matrix(
                    aMatrix.n_rows(),
                    aMatrix.pointers(),
                    aMatrix.indices(),
                    aMatrix.data(),
                    false );

            // call initialize function from parent
            Wrapper::initialize();
#endif
        }

//------------------------------------------------------------------------------

        void
        STRUMPACK::solve(
                SpMatrix & aMatrix,
                Vector< real > & aLHS,
                Vector< real > & aRHS )
        {
            BELFEM_ASSERT(  aMatrix.type() == SpMatrixType::CSR,
                            "Matrix must be of type CSR" );

#ifdef BELFEM_STRUMPACK


            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            if( aRHS.length() != aLHS.length() )
            {
                 aLHS.set_size( aRHS.length() );
            }

            if( ! this->is_initialized() )
            {
                this->initialize( aMatrix );
            }
            else
            {
                mSolver->update_matrix_values(
                        aMatrix.n_rows(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        false );
            }

            mSolver->factor();
            mSolver->solve( aRHS.data(), aLHS.data() );
#else
            BELFEM_ERROR( false, "We are not linked against STRUMPACK");
#endif
        }

//----------------------------------------------------------------
    }
}