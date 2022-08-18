//
// Created by christian on 8/16/22.
//
#include "commtools.hpp"
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
            if( mDistributor != nullptr )
            {
                delete mDistributor ;
            }

#ifdef BELFEM_STRUMPACK
            if( mSolver != nullptr )
            {
                delete mSolver;
            }
            if( mDistSolver != nullptr )
            {
                delete mDistSolver ;
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

            // call initialize function from parent
            Wrapper::initialize();

#ifdef BELFEM_STRUMPACK
            mArgC = gComm.arguments().size() ;

            mArgV = new StringList( mArgC );

            for( string tString : gComm.arguments( ) )
            {
                mArgV->push( tString );
            }

            if( gComm.size() == 1 )
            {
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

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
            }
            else if ( gComm.rank() == 0 )
            {
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );
                mDistributor = new StrumpackDistributor( aMatrix );
                mDistSolver = new strumpack::StrumpackSparseSolverMPIDist<real, int>( gComm.world(),
                        mArgC,
                        mArgV->data(),
                        gLog.info_level() >= 5 );

                mDistributor->send_values ( aMatrix );
                mDistSolver->set_distributed_csr_matrix(
                        mDistributor->n_rows(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        mDistributor->dist() );


            }
            else
            {
                mDistributor = new StrumpackDistributor();
                mDistSolver = new strumpack::StrumpackSparseSolverMPIDist<real, int>( gComm.world(),
                                                                                       mArgC,
                                                                                       mArgV->data(),
                                                                                       gLog.info_level() >= 5 );
                mDistributor->receive_values();
                mDistSolver->set_distributed_csr_matrix(
                        mDistributor->n_rows(),
                        mDistributor->pointers(),
                        mDistributor->indices(),
                        mDistributor->values(),
                        mDistributor->dist() );
            }
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

            if( aRHS.length() != aLHS.length() && gComm.rank() == 0 )
            {
                 aLHS.set_size( aRHS.length() );
            }

            if( gComm.size() == 1 )
            {
                if ( ! this->is_initialized() )
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
                mSolver->solve( aRHS.data(), aLHS.data());
            }
            else
            {
                if ( ! this->is_initialized() )
                {
                    this->initialize( aMatrix );
                }
                else if( gComm.rank() == 0 )
                {
                    mDistributor->send_values ( aMatrix );
                    mDistSolver->update_matrix_values(
                            mDistributor->n_rows(),
                            aMatrix.pointers(),
                            aMatrix.indices(),
                            aMatrix.data(),
                            mDistributor->dist() );
                }
                else
                {
                    mDistributor->receive_values();
                    mDistSolver->update_matrix_values(
                            mDistributor->n_rows(),
                            mDistributor->pointers(),
                            mDistributor->indices(),
                            mDistributor->values(),
                            mDistributor->dist() );
                };
                mDistributor->distribute_rhs( aRHS );
                comm_barrier();
                if( gComm.rank() == 0 )
                {
                    mDistSolver->solve( aRHS.data(), aLHS.data() );
                }
                else
                {
                    mDistSolver->solve( mDistributor->rhs(), mDistributor->lhs());
                }
                mDistributor->collect_lhs( aLHS );
            }
#else
            BELFEM_ERROR( false, "We are not linked against STRUMPACK");
#endif
        }

//----------------------------------------------------------------
    }
}