//
// Created by christian on 8/17/22.
//
#include "commtools.hpp"
#include "cl_SolverStrumpackDistributor.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        StrumpackDistributor::StrumpackDistributor( SpMatrix & aMatrix ) :
            mMyRank( comm_rank() ),
            mCommSize( comm_size() )
        {
            BELFEM_ASSERT( mMyRank == 0, "Rank must be 0 for this constructor" );

            this->send_sparsity_pattern( aMatrix ) ;

            mMyRhs.set_size( mMyNumRows, 0.0 );
            mMyLhs.set_size( mMyNumRows, 0.0 );
        }

//------------------------------------------------------------------------------

        StrumpackDistributor::StrumpackDistributor() :
                mMyRank( comm_rank() ),
                mCommSize( comm_size() )
        {
            BELFEM_ASSERT( mMyRank != 0, "Rank must not be 0 for this constructor" );

            this->receive_sparsity_pattern();

            mMyRhs.set_size( mMyNumRows, 0.0 );
            mMyLhs.set_size( mMyNumRows, 0.0 );
        }

//------------------------------------------------------------------------------

        void
        StrumpackDistributor::send_sparsity_pattern( SpMatrix & aMatrix )
        {
            BELFEM_ASSERT( mMyRank == 0, "Rank must be 0 in order to receive data" );

            // make sure that matrix is in C++ mode
            aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            // get number of rows per matrix
            mNumRowsOther = std::ceil( ( real ) aMatrix.n_rows() / ( real ) mCommSize );

            // get number of rows for first proc
            mMyNumRows = aMatrix.n_rows() - mNumRowsOther * ( mCommSize - 1 );
            //mMyPointers.set_size( mMyNumRows + 1 );
            //std::copy( aMatrix.pointers(), aMatrix.pointers() + mMyNumRows + 1, mMyPointers.data() );

            uint tOffset = mMyNumRows;
            for( proc_t p=1; p<mCommSize; ++p )
            {
                send( p, mNumRowsOther + 1, aMatrix.pointers() + tOffset );
                tOffset += mNumRowsOther ;
            }

            tOffset = 0 ;
            mNumNnzPerProc.set_size( comm_size() );
            uint tNumRows = mMyNumRows ;

            // get indices from matrix
            const int * tPointers = aMatrix.pointers() ;

            mDist.set_size( mCommSize+1 );
            mDist( 0 ) = 0 ;

            // loop over all procs
            for( proc_t p=0; p<mCommSize; ++p )
            {
                mNumNnzPerProc( p ) = 0 ;

                // loop over all rows of this proc
                for( int r=0; r<tNumRows; ++r )
                {
                    // count number of nnz
                    mNumNnzPerProc( p ) += tPointers[ r + tOffset + 1 ] - tPointers[ r + tOffset ];
                }

                // increment offset
                tOffset += tNumRows ;

                mDist( p+1 ) = tOffset ;

                // reset number of rows
                tNumRows = mNumRowsOther ;
            }




            // remember my nnz
            mMyNNZ = mNumNnzPerProc( 0 );

            //mMyIndices.set_size( mMyNNZ );
            //std::copy( aMatrix.indices(), aMatrix.indices() + mMyNNZ, mMyIndices.data() );

            // send indices
            tOffset = mNumNnzPerProc( 0 ) ;
            for( proc_t p=1; p<mCommSize; ++p )
            {
                send( p, mDist );
                send( p, mNumNnzPerProc( p ), aMatrix.indices() + tOffset );
                tOffset += mNumNnzPerProc( p );
            }
        }

//------------------------------------------------------------------------------

        void
        StrumpackDistributor::receive_sparsity_pattern()
        {
            BELFEM_ASSERT( mMyRank != 0, "Rank must not be 0 in order to receive data" );
            receive( 0 , mMyPointers );
            mMyNumRows = mMyPointers.length() - 1 ;
            mMyPointers -= mMyPointers( 0 );
            receive( 0, mDist );
            receive( 0, mMyIndices );
            mMyNNZ = mMyIndices.length() ;
        }

//------------------------------------------------------------------------------

        void
        StrumpackDistributor::send_values( const SpMatrix & aMatrix )
        {
            BELFEM_ASSERT( mMyRank == 0, "Rank must be 0 in order to send data" );

            // send indices
            int tOffset = mNumNnzPerProc( 0 ) ;
            for( proc_t p=1; p<mCommSize; ++p )
            {
                send( p, mNumNnzPerProc( p ), aMatrix.data() + tOffset );
                tOffset += mNumNnzPerProc( p );
            }
        }

//----------------------------------------------------------------------------

        void
        StrumpackDistributor::receive_values()
        {
            BELFEM_ASSERT( mMyRank != 0, "Rank must not be 0 in order to receive data" );
            receive( 0, mMyValues );
        }


//------------------------------------------------------------------------------

        void
        StrumpackDistributor::distribute_rhs( Vector< real> & aRHS )
        {
            if ( mMyRank == 0 )
            {
                int tOff = mMyNumRows ;
                for( uint p=1; p<mCommSize; ++p )
                {
                    send( p, mNumRowsOther, aRHS.data() + tOff );
                    tOff += mNumRowsOther ;
                }
            }
            else
            {
                receive( 0, mMyRhs );
            }
        }

//------------------------------------------------------------------------------

        void
        StrumpackDistributor::collect_lhs( Vector< real> & aLHS )
        {
            if ( mMyRank == 0 )
            {
                int tOff = mMyNumRows ;
                for( uint p=1; p<mCommSize; ++p )
                {
                    receive( p, aLHS.data() + tOff );
                    tOff += mNumRowsOther ;
                }
            }
            else
            {
                send( 0, mMyLhs );
            }
        }

//------------------------------------------------------------------------------
    }
}