/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_SolverDistMatrix.hpp"
#include "fn_create_graph_from_matrix.hpp"
#include "graphtools.hpp"
#include "fn_Graph_METIS.hpp"
#include "fn_Graph_ParMETIS.hpp"
#include "fn_Graph_PTSCOTCH.hpp"
namespace belfem
{
    namespace sparse
    {

        DistMatrix::DistMatrix( const SolverParameters * aParams, SpMatrix * aMatrix ) :
                   mCommRank( comm_rank() ),
                   mCommSize( comm_size() ),
                   mPermutationSwitch( this->set_permutation_switch( aParams ) ),
                   mReorderingMethod( this->set_reordering_method( aParams ) )
        {
            // aMatrix is only populated on rank 0 (master rank assembles the full matrix)
            // Non-root ranks can pass nullptr or an empty SpMatrix object
            BELFEM_ASSERT( aMatrix != nullptr || mCommRank != 0 , "Matrix must not be null on rank 0" );

            this->create_matrix( aMatrix );
            this->determine_sizes_and_offsets();
        }

        DistMatrix::~DistMatrix()
        {
            if ( mPermutationSwitch && mMatrix != 0 ) delete mMatrix;
        }

        bool
        DistMatrix::set_permutation_switch( const SolverParameters * aParams )
        {
            int tSwitch = 0 ;
            if ( gComm.rank() == 0 )
            {
                // STRUMPACK runs its own nested dissection and redistributes
                // internally; a BELFEM-side pre-permutation would only
                // reorder the matrix twice
                tSwitch =    aParams->type() != SolverType::STRUMPACK
                          && (    aParams->reordering_method() == ReorderingMethod::METIS
                               || aParams->reordering_method() == ReorderingMethod::PARMETIS
                               || aParams->reordering_method() == ReorderingMethod::PTSCOTCH ) ;
            }
            broadcast( tSwitch );
            return tSwitch == 1 ;
        }

        ReorderingMethod
        DistMatrix::set_reordering_method( const SolverParameters * aParams )
        {
            int tMethod = static_cast< int >( ReorderingMethod::METIS );
            if ( gComm.rank() == 0 )
            {
                tMethod = static_cast< int >( aParams->reordering_method() );
            }
            broadcast( tMethod );
            return static_cast< ReorderingMethod >( tMethod );
        }

        void
        DistMatrix::create_matrix( SpMatrix * aMatrix )
        {
            // Every rank stays in this function up to and including
            // order_graph(): the parallel orderings are MPI-collective, and
            // the old root-only early return here would hang them. aMatrix
            // is still touched on root only
            if ( mCommRank == 0 )
            {
                mNumRows = aMatrix->n_rows();
            }
            comm_barrier() ;
            broadcast( mNumRows );

            if ( mCommRank == 0 )
            {
                // enforce zero based indexing
                aMatrix->set_indexing_base( SpMatrixIndexingBase::Cpp );
            }

            if ( ! mPermutationSwitch )
            {
                if ( mCommRank == 0 )
                {
                    mMatrix = aMatrix;
                }
                return;
            }

            // root holds the complete graph, non-root ranks an empty one:
            // the input contract of parmetis_nd() / ptscotch_nd()
            Graph tGraph ;

            if ( mCommRank == 0 )
            {
                // Step 1: Create graph from original matrix
                sparse::create_graph_from_matrix( *aMatrix, tGraph );
            }

            // Step 2: Apply reordering to the graph ( collective )
            this->order_graph( tGraph );

            if ( mCommRank != 0 )
            {
                return;
            }

            mForwardPermutation.set_size( tGraph.size() );
            mBackwardPermutation.set_size( tGraph.size() );

            for( graph::Vertex * tVertex : tGraph )
            {
                mForwardPermutation( tVertex->index() )     = tVertex->id();
                mBackwardPermutation( tVertex->id() ) = tVertex->index();
            }

            // Step 4: Create the permuted SpMatrix
            mMatrix = new SpMatrix( tGraph,
                aMatrix->type(),
                aMatrix->n_rows(),
                aMatrix->n_cols() );

            // Step 5: Build index permutation AFTER SpMatrix creation
            // so it matches the actual matrix structure (post-sort)
            this->compute_index_permutation( aMatrix );

            // Clean up graph now that we're done with it
            for ( graph::Vertex * tVertex : tGraph )
            {
                delete tVertex;
            }
            tGraph.clear();
        }

        void
        DistMatrix::order_graph( Graph & aGraph )
        {
            switch ( mReorderingMethod )
            {
                case ReorderingMethod::PARMETIS :
                {
                    if ( mCommRank == 0 )
                    {
                        graph::block_distribution( aGraph, mCommSize );
                    }
                    graph::parmetis_nd( aGraph );
                    break;
                }
                case ReorderingMethod::PTSCOTCH :
                {
                    if ( mCommRank == 0 )
                    {
                        graph::block_distribution( aGraph, mCommSize );
                    }
                    graph::ptscotch_nd( aGraph );
                    break;
                }
                default :
                {
                    // serial nested dissection on root; a no-op on the
                    // empty non-root graph
                    graph::metis_ndp( aGraph, mCommSize );
                }
            }
        }

        void
        DistMatrix::compute_index_permutation( const SpMatrix * aMatrix )
        {
            mIndexPermutation.set_size( aMatrix->number_of_nonzeros() );

            bool tIsCSR = aMatrix->type() == SpMatrixType::CSR;
            const int_t * tOldPointers     = aMatrix->pointers();
            const int_t * tOldIndices      = aMatrix->indices();

            // Now iterate using mMatrix's structure (which reflects the sorted graph order)
            const int_t * tNewPointers = mMatrix->pointers();
            const int_t * tNewIndices  = mMatrix->indices();
            int_t tN = mMatrix->n_rows();

            index_t tNewPos = 0;
            for( int_t iNew = 0; iNew < tN; ++iNew )
            {
                // Get the original row index for this new row
                int_t iOld = mForwardPermutation( iNew );

                // Iterate through columns in the NEW matrix order
                for( int_t k = tNewPointers[iNew]; k < tNewPointers[iNew + 1]; ++k )
                {
                    int_t jNew = tNewIndices[k];  // New column index
                    int_t jOld = mForwardPermutation( jNew );  // Original column index

                    // Find position in original matrix
                    int_t iSearch = tIsCSR ? iOld : jOld;
                    int_t jSearch = tIsCSR ? jOld : iOld;

                    int_t start = tOldPointers[iSearch];
                    int_t end = tOldPointers[iSearch + 1];

                    bool found = false;
                    
                    for( int_t pos = start; pos < end; ++pos )
                    {
                        if( tOldIndices[pos] == jSearch )
                        {
                            mIndexPermutation( tNewPos ) = pos;
                            found = true;
                            break;
                        }
                    }

                    BELFEM_ERROR( found,
                        "Could not find old position for entry (%d, %d) -> (%d, %d)",
                        (int)iNew, (int)jNew, (int)iOld, (int)jOld );
                    ++tNewPos;
                }
            }
        }

        void
        DistMatrix::determine_sizes_and_offsets()
        {
            if ( mCommRank != 0 )
            {
                // receive as int_t to match the type sent by distribute(mSizes)
                int_t tNumRows;
                receive( tNumRows );
                mMyNumRows = tNumRows;
                mMyLhs.set_size( mMyNumRows );
                mMyRhs.set_size( mMyNumRows );
                return;
            }

            index_t tN = mMatrix->n_rows();
            index_t tM = mCommSize ;
            index_t tMod = tN % tM ;
            index_t tDiv = tN / tM ;
            index_t tSplit = tM - tMod ;

            mSizes.set_size( mCommSize, 0 );
            for ( index_t p=0; p<tSplit; ++p )
            {
                mSizes( p ) = tDiv ;
            }
            for ( index_t p=tSplit; p<tM; ++p )
            {
                mSizes( p ) = tDiv + 1 ;
            }

            distribute( mSizes );
            mMyNumRows = mSizes( 0 );
            mMyLhs.set_size( mMyNumRows );
            mMyRhs.set_size( mMyNumRows );

            mDist.set_size( tM + 1, 0 );

            for ( index_t p=0; p<tM; ++p )
            {
                mDist( p+1 ) = mDist( p ) + mSizes( p );
            }

            mOffsets.set_size( tM + 1, 0 );

            const int_t * tRawPtr = mMatrix->pointers();

            for ( index_t p=0; p<=tM; ++p )
            {
                mOffsets( p ) = tRawPtr[ mDist( p )] ;
            }

            // set the number of rows for rank 0



        }

        void
        DistMatrix::distribute_rhs( const Vector< real > & aRhs )
        {
            if ( mCommRank == 0 )
            {
                BELFEM_ASSERT( aRhs.length() == mMatrix->n_rows(),
                    "RHS size mismatch ( is %lu, expect %lu )",
                    ( long unsigned int ) aRhs.length(),
                    ( long unsigned int ) mMatrix->n_rows() );

                mRhs.set_size( mMatrix->n_rows(), 0 );
                mMyRhs.set_size( mMyNumRows );
                if ( mPermutationSwitch )
                {
                    index_t tCount = 0;
                    for ( index_t k : mForwardPermutation )
                    {
                        mRhs( tCount++ ) = aRhs( k );
                    }
                    for ( index_t i = 0; i < mMyNumRows; ++i )
                    {
                        mMyRhs( i ) = aRhs( mForwardPermutation( i ) );
                    }

                    comm_barrier();
                    distribute( mRhs.data(), mDist );
                }
                else
                {
                    for ( index_t k=0; k<mMyNumRows; ++k )
                    {
                        mMyRhs( k ) = aRhs( k );
                    }
                    comm_barrier();
                    distribute( aRhs.data(), mDist );
                }
            }
            else
            {
                comm_barrier();
                receive( mMyRhs );
            }
        }

        void
        DistMatrix::distribute_lhs( const Vector< real > & aLhs )
        {
            if ( mCommRank == 0 )
            {
                BELFEM_ASSERT( aLhs.length() == mMatrix->n_rows(),
                    "RHS size mismatch ( is %lu, expect %lu )",
                    ( long unsigned int ) aLhs.length(),
                    ( long unsigned int ) mMatrix->n_rows() );

                mLhs.set_size( mMatrix->n_rows(), 0 );
                mMyLhs.set_size( mMyNumRows );

                if ( mPermutationSwitch )
                {
                    index_t tCount = 0;
                    for ( index_t k : mForwardPermutation )
                    {
                        mLhs( tCount++ ) = aLhs( k );
                    }
                    for ( index_t i = 0; i < mMyNumRows; ++i )
                    {
                        mMyLhs( i ) = aLhs( mForwardPermutation( i ) );
                    }
                    comm_barrier();
                    distribute( mLhs.data(), mDist );
                }
                else
                {
                    for ( index_t k=0; k<mMyNumRows; ++k )
                    {
                        mMyLhs( k ) = aLhs( k );
                    }
                    comm_barrier();
                    distribute( aLhs.data(), mDist );
                }
            }
            else
            {
                comm_barrier();
                receive( mMyLhs );
            }
        }

        void
        DistMatrix::collect_lhs( Vector< real > & aLhs )
        {
            if ( mCommRank == 0 )
            {
                if ( ! mPermutationSwitch )
                {
                    aLhs.set_size( mMatrix->n_rows() );

                    // copy rank 0's local slice to full vector before collecting
                    for ( index_t i=0; i<mMyNumRows; ++i )
                    {
                        aLhs( i ) = mMyLhs( i );
                    }
                    comm_barrier();
                    collect( aLhs.data(), mDist );
                }
                else
                {
                    mLhs.set_size( mMatrix->n_rows() );
                    aLhs.set_size( mMatrix->n_rows() );
                    for ( index_t i=0; i<mMyNumRows; ++i )
                    {
                        mLhs( i ) = mMyLhs( i );
                    }
                    comm_barrier();
                    collect( mLhs.data(), mDist );

                    index_t tCount = 0 ;
                    for (  index_t k : mBackwardPermutation )
                    {
                        aLhs( tCount++ ) = mLhs( k );
                    }
                }
            }
            else
            {
                comm_barrier();
                send( mMyLhs );
            }
        }

    }
}