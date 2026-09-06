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

#ifndef BELFEM_GRAPHTOOLS_HPP
#define BELFEM_GRAPHTOOLS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Vector.hpp"
#include "cl_Graph_Vertex.hpp"
#include "op_Graph_Vertex_Index.hpp"
#include "commtools.hpp"


namespace belfem
{
    namespace graph
    {

//------------------------------------------------------------------------------

        template< typename T >
        void
        apply_graph_permutation(
                 Graph & aGraph,
           const Vector< T >  & aPermutation )
        {
            // Apply the permutation to the graphx
            index_t tCount = 0;
            for( Vertex * tVertex : aGraph )
            {
                tVertex->set_index( static_cast< index_t> ( aPermutation( tCount++ ) ) );
            }

            // Sort the graph by the new indices
            sort( aGraph, opVertexIndex );

            DynamicBitset tBitset( aGraph.size() );

            Cell< index_t > tIndices ;
            for ( Vertex * tVertex : aGraph )
            {
                tBitset.reset();
                for ( uint k=0; k<tVertex->number_of_vertices(); ++k )
                {
                    tBitset.set( tVertex->vertex( k )->index() );
                }
                if ( tVertex->is_flagged() )
                {
                    tBitset.set( tVertex->index() );
                }
                tBitset.where( tIndices );
                tVertex->reset_vertex_container();
                tVertex->init_vertex_container( tIndices.size() );
                for ( index_t k : tIndices )
                {
                    tVertex->insert_vertex( aGraph( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        template< typename T >
        void
        build_graph_adjacency(
            Graph & aGraph,
            Vector< T >  & aVertices,
            Vector< T >  & aEdges )
        {
            T tNumVertices = static_cast< T >( aGraph.size() );

            // First pass: ensure continuous indices
            T tCount = 0;

            for( Vertex * tVertex : aGraph )
            {
                tVertex->set_index( tCount++ );
            }

            // Second pass: count total edges. Self-loops are dropped here and
            // in the fill below: a graph built from a matrix carries one per
            // diagonal entry, and METIS_NodeNDP corrupts its heap on them
            // ( reproduced 2026-09-04, FM_2WayNodeRefine1Sided → rpqDestroy ).
            // METIS and SCOTCH define their input as loop-free; ParMETIS passes
            // loops through unchecked. The Graph itself keeps them: DistMatrix
            // builds the permuted SpMatrix from the Graph and needs the diagonal
            tCount = 0;
            for( Vertex * tVertex : aGraph )
            {
                for( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                {
                    if( tVertex->vertex( k ) != tVertex )
                    {
                        ++tCount;
                    }
                }
            }

            // Allocate METIS arrays
            aVertices.set_size( tNumVertices + 1, 0 );
            // Always allocate at least 1 element to avoid NULL pointers
            if( tCount > 0 )
            {
                aEdges.set_size( tCount, 0 );
            }
            else
            {
                aEdges.set_size( 1, 0 );
            }

            // Build the adjacency structure (CSR format)
            tCount = 0;

            for( Vertex * tVertex : aGraph )
            {
                aVertices( tVertex->index() ) = tCount;

                for( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                {
                    if( tVertex->vertex( k ) != tVertex )
                    {
                        aEdges( tCount++ ) = static_cast< T >( tVertex->vertex( k )->index() );
                    }
                }
            }
            aVertices( tNumVertices ) = tCount;
        }

        /**
         * assign owners to an ownerless graph in contiguous blocks of graph
         * order: the first aCommSize - ( N mod aCommSize ) ranks take
         * N / aCommSize vertices, the remaining ranks one more ( the split
         * DistMatrix uses for its rows ). A working distribution for the
         * parallel nested-dissection wrappers, which need every vertex owned
         * before build_pargraph_adjacency(); not a partition. Two loops,
         * never k / ( N / P ): that quotient is zero when N < P.
         */
        inline void
        block_distribution( Graph & aGraph, const proc_t aCommSize )
        {
            BELFEM_ASSERT( aCommSize > 0, "block_distribution: comm size %d", ( int ) aCommSize );

            const index_t tN     = aGraph.size();
            const index_t tP     = static_cast< index_t >( aCommSize );
            const index_t tDiv   = tN / tP;
            const index_t tSplit = tP - ( tN % tP );

            index_t tCount = 0;
            for ( index_t p = 0; p < tSplit; ++p )
            {
                for ( index_t k = 0; k < tDiv; ++k )
                {
                    aGraph( tCount++ )->set_owner( static_cast< proc_t >( p ) );
                }
            }
            for ( index_t p = tSplit; p < tP; ++p )
            {
                for ( index_t k = 0; k < tDiv + 1; ++k )
                {
                    aGraph( tCount++ )->set_owner( static_cast< proc_t >( p ) );
                }
            }
            BELFEM_ASSERT( tCount == tN, "block_distribution: assigned %lu of %lu vertices",
                ( long unsigned int ) tCount, ( long unsigned int ) tN );
        }

        template < typename  T >
        void
        build_pargraph_adjacency(
            Graph & aGraph,
                  Vector< T > & aDistribution,
            Cell< Vector< T > > & aVertices,
            Cell< Vector< T > > & aEdges )
        {
                proc_t tCommSize = comm_size();
                BELFEM_ASSERT( comm_rank() == 0, "This build_pargraph_adjacency() can only be called by root" );

                // first we need to reorder the graph based on its owner
                Vector< T > tCount( tCommSize, 0 );
                for ( Vertex * tVertex : aGraph )
                {
                    tVertex->unflag();

                    // unowned vertices (e.g. circuit dofs) go to root. Two sentinels are in
                    // use: the Kernel's comm_size() marker and the gNoOwner default.
                    proc_t tOwner = tVertex->owner();
                    if ( tOwner == tCommSize || tOwner == gNoOwner )
                    {
                        tOwner = 0;
                        tVertex->set_owner( tOwner );
                    }
                    BELFEM_ERROR( tOwner >= 0 && tOwner < tCommSize,
                        "Vertex %lu has owner %d outside [0, %d)",
                        ( long unsigned int ) tVertex->id(), ( int ) tOwner, ( int ) tCommSize );

                    ++tCount( tOwner );
                }

                aDistribution.set_size( comm_size()+1, 0 );
                for ( proc_t p = 0; p < comm_size(); ++p )
                {
                    aDistribution( p+1 ) = aDistribution( p ) + tCount( p );
                }

                // now we need to reorder the graph based on its owner
                Graph tGraph( aGraph.size(), nullptr );
                tCount.fill( 0 );
                for ( Vertex * tVertex : aGraph )
                {
                    tGraph( aDistribution( tVertex->owner() ) + tCount( tVertex->owner() )++ ) = tVertex;
                }
                aGraph = std::move( tGraph );
                tGraph.clear();
#if defined( DEBUG )
                for ( Vertex * tVertex : aGraph )
                {
                    for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                    {
                        tVertex->vertex( k )->set_index( gNoIndex );
                    }
                    tVertex->set_index( gNoIndex );
                }
#endif

                // next, we make sure that the indices are consecutive
                index_t tIndex = 0;
                for ( Vertex * tVertex : aGraph )
                {
                    tVertex->set_index( tIndex++ );
                }

#if defined( DEBUG )
                for ( Vertex * tVertex : aGraph )
                {
                    for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                    {
                        BELFEM_ERROR( tVertex->vertex( k )->index() != gNoIndex,
                            "Vertex index of %lu (owner %u) is not set",
                            ( long unsigned int ) tVertex->vertex( k )->id() ,
                            ( uint ) tVertex->vertex( k)->owner() );
                    }
                }
#endif

                // Count edges per processor, self-loops dropped ( see
                // build_graph_adjacency for why )
                Vector< T > tNNZ( tCommSize, 0 );
                for ( Vertex * tVertex : aGraph )
                {
                    for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                    {
                        if ( tVertex->vertex( k ) != tVertex )
                        {
                            ++tNNZ( tVertex->owner() );
                        }
                    }
                }

            aVertices.set_size( tCommSize, {} );
            aEdges.set_size( tCommSize, {} );

            for ( proc_t p = 0; p < tCommSize; ++p )
            {
                aVertices( p ).set_size( tCount( p ) + 1 );

                // Always allocate at least 1 element to avoid NULL pointers in
                // ParMETIS. Branch on the EDGE count: a slice whose vertices
                // carry only self-loops has vertices but, after the loop drop,
                // no edges, and must get the placeholder too
                if ( tNNZ( p ) > 0 )
                {
                    aEdges( p ).set_size( tNNZ( p ) );
                }
                else
                {
                    aEdges( p ).set_size( 1, 0 );
                }
            }

            Vector< T > tXcount( tCommSize, 0 );
            Vector< T > tAcount( tCommSize, 0 );

            for( Vertex * tVertex : aGraph )
            {
                proc_t tOwner = tVertex->owner();
                aVertices( tOwner )( tXcount( tOwner )++ ) = tAcount( tOwner );

                for ( uint k=0; k<tVertex->number_of_vertices(); ++k )
                {
                    if ( tVertex->vertex( k ) != tVertex )
                    {
                        aEdges( tOwner )( tAcount( tOwner )++ ) = static_cast< T >( tVertex->vertex( k )->index() );
                    }
                }
            }

            BELFEM_ASSERT( tAcount == tNNZ, "The number of non-zeros in the adjacency matrix is not correct" );

            for ( proc_t p = 0; p < tCommSize; ++p )
            {
                aVertices( p )( tCount( p ) ) = tNNZ( p );
            }
        }

    }
}
#endif //BELFEM_GRAPHTOOLS_HPP