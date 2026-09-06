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

#include <string>

#include "commtools.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "cl_SolverSTRUMPACK.hpp"
#include "assert.hpp"
// Externally Defined Global Communicator
extern belfem::Logger gLog;
extern belfem::Communicator gComm;

#ifdef BELFEM_STRUMPACK
#include "strumpacktools.hpp"
#endif

namespace belfem
{
    namespace solver
    {
//----------------------------------------------------------------

        STRUMPACK::STRUMPACK( const SolverParameters * aParams ) :
                Wrapper( "STRUMPACK" , true ),
                mParams( aParams )
        {

        }

//----------------------------------------------------------------

        STRUMPACK::~STRUMPACK()
        {

           STRUMPACK::free();

        }

//----------------------------------------------------------------

        void
        STRUMPACK::free()
        {
#ifdef BELFEM_STRUMPACK
            if( mDistMatrix != nullptr )
            {
                delete mDistMatrix ;
                mDistMatrix = nullptr ;
            }

            if( mSolver != nullptr )
            {
                delete mSolver;
                mSolver = nullptr ;
            }
            if( mDistSolver != nullptr )
            {
                delete mDistSolver ;
                mDistSolver = nullptr ;
            }

            if( mArgV != nullptr )
            {
                delete mArgV ;
                mArgV = nullptr ;
            }
#endif
            Wrapper::free();
        }

//------------------------------------------------------------------------------

        void
        STRUMPACK::initialize( SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode,
                    const int_t aNumRhsColumns )
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

            // check if we need to build a parallel solver
            if( gComm.size() == 1 )
            {
                // make sure that the matrix is set to zero-indexing
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

                // create the solver
                mSolver = new strumpack::StrumpackSparseSolver< real, int >(
                        mArgC,
                        mArgV->data(),
                        gLog.info_level() >= 5 );

                // convert options to strumpack format
                sparse::set_strumpack_options( *mParams, mSolver->options(), aMatrix.n_rows() );

                // allow user to override setings using the command line
                mSolver->options().set_from_command_line( mArgC, mArgV->data() );

                // link the matrix to the solver
                mSolver->set_csr_matrix(
                        aMatrix.n_rows(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        false );

            }
            else
            {
                // Note: aMatrix is only populated on rank 0 (master rank assembles the full matrix)
                // Non-root ranks never access it (DistMatrix handles distribution internally)
                if ( this->rank() == 0 )
                {
                    aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );
                }

                switch ( mParams->distributed_matrix_type() )
                {
                    case DistributedMatrixType::CSR :
                    {
                        mDistMatrix = new sparse::StrumpackCSR( mParams, &aMatrix );
                        break ;
                    }
                    case DistributedMatrixType::AIJ :
                    {
                        mDistMatrix = new sparse::StrumpackAIJ( mParams, &aMatrix );
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false,
                            "Invalid distributed matrix type for STRUMPACK %s",
                            to_string( mParams->distributed_matrix_type() ).c_str()
                            );
                    }
                }

                // distribute actual matrix values
                // (DistMatrix constructor only distributes sparsity pattern with zero values)
                mDistMatrix->distribute_values( &aMatrix );

                // create a parallel solver
                mDistSolver = new strumpack::StrumpackSparseSolverMPIDist<real, int>( gComm.world(),
                                                                                      mArgC,
                                                                                      mArgV->data(),
                                                                                      gLog.info_level() >= 5 );
                index_t N = mDistMatrix->size();

                // warn if the MPI ranks oversubscribe this node's cores with OpenMP threads
                this->hatch_turtle();

                // convert options to strumpack format
                sparse::set_strumpack_options( *mParams, mDistSolver->options(), N );

                // allow user to override setings using the command line
                mDistSolver->options().set_from_command_line( mArgC, mArgV->data() );

                // Link this distributed matrix to the solver
                // Use compile-time switch to choose format
                if ( mParams->distributed_matrix_type() == DistributedMatrixType::CSR )
                {
                    auto * tCSR = static_cast< sparse::StrumpackCSR * >( mDistMatrix );

                    mDistSolver->set_distributed_csr_matrix(
                        tCSR->n_rows(),
                        tCSR->pointers(),
                        tCSR->indices(),
                        tCSR->values(),
                        tCSR->dist(),
                        false );  // not symmetric
                }
                else
                {
                    auto * tAIJ = static_cast< sparse::StrumpackAIJ * >( mDistMatrix );
                    // MPIAIJ format (now with proper garray!)
                    mDistSolver->set_MPIAIJ_matrix(
                        tAIJ->n_rows(),
                        tAIJ->diagonal_pointers(),
                        tAIJ->diagonal_indices(),
                        tAIJ->diagonal_values(),
                        tAIJ->offdiagonal_pointers(),
                        tAIJ->offdiagonal_indices(),
                        tAIJ->offdiagonal_values(),
                        tAIJ->garray() );
                }

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


#ifdef BELFEM_STRUMPACK

            // on the main proc, make sure that the LHS vector is properly assigned

            if ( gComm.rank() == 0 )
            {
                BELFEM_ASSERT(  aMatrix.type() == SpMatrixType::CSR,
                      "Matrix must be of type CSR" );

                // make sure that matrix is zero-indexed
                aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

                if( aRHS.length() != aLHS.length() )
                {
                    aLHS.set_size( aRHS.length() );
                }
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
                    mSolver->update_matrix_values(
                        aMatrix.n_rows(),
                        aMatrix.pointers(),
                        aMatrix.indices(),
                        aMatrix.data(),
                        false );
                }

                strumpack::ReturnCode tReturnCode
                    = mSolver->solve( aRHS.data(), aLHS.data(), mParams->use_initial_guess() );

                if ( tReturnCode != strumpack::ReturnCode::SUCCESS )
                {
                    // soft-fail contract: record and return; the controller
                    // treats the event as a failed trial ( cf. MUMPS )
                    if ( this->soft_fail() )
                    {
                        this->flag_failure() ;
                        message( InfoLevel::Minimal,
                            "STRUMPACK solve failed softly ( code %d ) - handing back to the controller",
                            static_cast<int>( tReturnCode ) );
                        return ;
                    }

                    BELFEM_ERROR( false,
                                 "STRUMPACK serial solve failed with code %d : %s",
                                  static_cast<int>(tReturnCode),
                                  sparse::strumpack_message( tReturnCode ).c_str()
                                   );
                }

            }
            else // run in parallel mode
            {
                // check if we need to initialize the solver
                if ( ! this->is_initialized() )
                {
                    this->initialize( aMatrix );
                }
                else
                {

                    mDistMatrix->distribute_values( &aMatrix );

                    if ( mParams->distributed_matrix_type() == DistributedMatrixType::CSR )
                    {
                        sparse::StrumpackCSR * tCSR = static_cast< sparse::StrumpackCSR * >( mDistMatrix );
                        mDistSolver->update_matrix_values(
                            tCSR->n_rows(),
                            tCSR->pointers(),
                            tCSR->indices(),
                            tCSR->values(),
                            tCSR->dist(),
                            false );  // not symmetric
                    }
                    else
                    {
                        sparse::StrumpackAIJ * tAIJ = static_cast< sparse::StrumpackAIJ * >( mDistMatrix );
                        mDistSolver->update_MPIAIJ_matrix_values(
                            tAIJ->n_rows(),
                            tAIJ->diagonal_pointers(),
                            tAIJ->diagonal_indices(),
                            tAIJ->diagonal_values(),
                            tAIJ->offdiagonal_pointers(),
                            tAIJ->offdiagonal_indices(),
                            tAIJ->offdiagonal_values(),
                            tAIJ->garray() );
                    }

                }

                mDistMatrix->distribute_rhs( aRHS );

                // pass the caller's guess to the local vectors
                // ( STRUMPACK only uses it with an iterative krylov method )
                if ( mParams->use_initial_guess() )
                {
                    mDistMatrix->distribute_lhs( aLHS );
                }

                comm_barrier();

                strumpack::ReturnCode tReturnCode
                    = mDistSolver->solve( mDistMatrix->rhs(), mDistMatrix->lhs(), mParams->use_initial_guess() );

                // synchronize the verdict: STRUMPACK propagates front errors
                // locally and does NOT reduce the return code across ranks
                // ( verified in the STRUMPACK sources, round-5 audit ), so
                // codes can diverge; a split between the soft return and
                // collect_lhs would deadlock
                int tFailLocal  = ( tReturnCode == strumpack::ReturnCode::SUCCESS ) ? 0 : 1 ;
                int tFailGlobal = tFailLocal ;
                allreduce( &tFailLocal, &tFailGlobal, 1 );

                if ( tFailGlobal != 0 )
                {
                    // soft-fail contract: all ranks record and return
                    // uniformly; the controller treats the event as a failed
                    // trial ( cf. MUMPS )
                    if ( this->soft_fail() )
                    {
                        this->flag_failure() ;
                        // the local code can read SUCCESS on ranks where the
                        // failure happened remotely ( codes are not reduced
                        // by STRUMPACK — hence the Allreduce above )
                        if ( this->rank() == 0 )
                        {
                            message( InfoLevel::Minimal,
                                "STRUMPACK solve failed softly on at least one rank ( local code %d ) - handing back to the controller",
                                static_cast<int>( tReturnCode ) );
                        }
                        comm_barrier();
                        return ;
                    }

                    BELFEM_ERROR( false,
                                 "STRUMPACK parallel solve failed ( local code on rank %d is %d ) : %s",
                                 this->rank(), static_cast<int>(tReturnCode),
                                 sparse::strumpack_message( tReturnCode ).c_str()
                                  );
                }

                mDistMatrix->collect_lhs( aLHS );

                // wait for all procs
                comm_barrier();
            }
#else
            BELFEM_ERROR( false, "We are not linked against STRUMPACK");
#endif
        }

//----------------------------------------------------------------
    }
}
