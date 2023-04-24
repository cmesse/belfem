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
#ifdef BELFEM_STRUMPACK
            if( mDistributor != nullptr )
            {
                delete mDistributor ;
            }

            if( mSolver != nullptr )
            {
                delete mSolver;
            }
            if( mDistSolver != nullptr )
            {
                delete mDistSolver ;
            }

            if( mArgV != nullptr )
            {
                delete mArgV ;
            }
#endif
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

            // get the arguments from the communicator
            // we need them to set custom settings for STRUMPACK
            mArgC = gComm.arguments().size() ;
            mArgV = new StringList( mArgC );
            for( string tString : gComm.arguments( ) )
            {
                mArgV->push( tString );
            }

            // make sure that the matrix is set to zero-indexing
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            // check if we need to build a parallel solver
            if( gComm.size() == 1 )
            {
                // create the solver
                mSolver = new strumpack::StrumpackSparseSolver< double, int >(
                        mArgC,
                        mArgV->data(),
                        gLog.info_level() >= 5 );

                // set the options from the command line
                mSolver->options().set_from_command_line( mArgC, mArgV->data() );

                // link the matrix to the solver
                mSolver->set_csr_matrix(
                        aMatrix.n_rows(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        false );
            }
            else if ( gComm.rank() == 0 )
            {
                // build the distribution object
                mDistributor = new StrumpackDistributor( aMatrix );

                // split up the matrix amongt the other procs
                mDistributor->send_values ( aMatrix );

                // create a parallel solver
                mDistSolver = new strumpack::StrumpackSparseSolverMPIDist<real, int>( gComm.world(),
                                                                                      mArgC,
                                                                                      mArgV->data(),
                                                                                      gLog.info_level() >= 5 );

                // link this distributed matrix to the solver
                mDistSolver->set_distributed_csr_matrix(
                        mDistributor->n_rows(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        mDistributor->dist() );
            }
            else
            {
                // build the distribution object
                mDistributor = new StrumpackDistributor();

                // get the matrix components form the main proc
                mDistributor->receive_values();

                // create a parallel solver
                mDistSolver = new strumpack::StrumpackSparseSolverMPIDist<real, int>( gComm.world(),
                                                                                       mArgC,
                                                                                       mArgV->data(),
                                                                                       gLog.info_level() >= 5 );
                // link this distributed matrix to the solver
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

            // make sure that matrix is zero-indexed
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            // on the main proc, make sure that the LHS vector is properly assigned
            if( aRHS.length() != aLHS.length() && gComm.rank() == 0 )
            {
                 aLHS.set_size( aRHS.length() );
            }

            // check if we are in serial or parallel mode
            if( gComm.size() == 1 )
            {
                // check if we need to initialize the solver
                if ( ! this->is_initialized() )
                {
                    this->initialize( aMatrix );
                }
                else
                {
                    // update the values from the matrix
                    mSolver->update_matrix_values(
                            aMatrix.n_rows(),
                            aMatrix.pointers(),
                            aMatrix.indices(),
                            aMatrix.data(),
                            false );
                }

                // solve the system
                mSolver->solve( aRHS.data(), aLHS.data() );
            }
            else // run in parallel mode
            {
                // check if we need to initialize the solver
                if ( ! this->is_initialized() )
                {
                    this->initialize( aMatrix );
                }
                else if( gComm.rank() == 0 )
                {
                    // update the matrix values
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
                    // update the matrix values
                    mDistributor->receive_values();
                    mDistSolver->update_matrix_values(
                            mDistributor->n_rows(),
                            mDistributor->pointers(),
                            mDistributor->indices(),
                            mDistributor->values(),
                            mDistributor->dist() );
                };

                //  distribute the right hand side
                mDistributor->distribute_rhs( aRHS );

                // wait for all procs
                comm_barrier();

                // solve the system
                if( gComm.rank() == 0 )
                {
                    mDistSolver->solve( aRHS.data(), aLHS.data() );
                }
                else
                {
                    mDistSolver->solve( mDistributor->rhs(), mDistributor->lhs() );
                }

                // unite the left hand side on the master prock
                mDistributor->collect_lhs( aLHS );
            }
#else
            BELFEM_ERROR( false, "We are not linked against STRUMPACK");
#endif
        }

//----------------------------------------------------------------
    }
}
