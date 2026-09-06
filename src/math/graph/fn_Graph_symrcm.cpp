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

#include "typedefs.hpp"

#include "fn_Graph_find_pseudo_peripheral_vertex.hpp"
#include "fn_Graph_symrcm.hpp"
#include "fn_Graph_bfs.hpp"
#include "op_Graph_Vertex_Index.hpp"
#include "op_Graph_Vertex_Degree.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Queue.hpp"
namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        void
        symrcm( Graph & aGraph, Vertex * aStart )
        {
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return;
            }

            index_t tCount = 0 ;
            for ( Vertex * tVertex: aGraph )
            {
                tVertex->set_index( tCount++ );
            }

            // BFS over all components: sets the vertex levels. Degrees are
            // the stored adjacency counts read below, not a BFS product
            bfs( aGraph );

            // Find maximum degree using the stored values
            uint tMaxDegree = 0;
            for( index_t i = 0; i < aGraph.size(); ++i )
            {
                uint tDegree = aGraph( i )->number_of_vertices();
                if( tDegree > tMaxDegree )
                {
                    tMaxDegree = tDegree;
                }
            }

            // Find starting vertex if not provided
            if( aStart == nullptr )
            {
                // First, find vertices with minimum degree as potential starts
                uint tMinDegree = tMaxDegree;
                Graph tMinDegreeVertices;

                for( index_t i = 0; i < aGraph.size(); ++i )
                {
                    uint tDegree = aGraph( i )->number_of_vertices();
                    if( tDegree < tMinDegree )
                    {
                        tMinDegree = tDegree;
                        tMinDegreeVertices.clear();
                        tMinDegreeVertices.push( aGraph( i ) );
                    }
                    else if( tDegree == tMinDegree )
                    {
                        tMinDegreeVertices.push( aGraph( i ) );
                    }
                }

                // Among minimum degree vertices, find pseudo-peripheral
                if( tMinDegreeVertices.size() == 1 )
                {
                    aStart = find_pseudo_peripheral_vertex( aGraph, tMinDegreeVertices( 0 ) );
                }
                else
                {
                    aStart = find_pseudo_peripheral_vertex( aGraph );
                }
            }

            // Initialize result and tracking structures
            Cell< index_t > tPermutation;
            tPermutation.reserve( aGraph.size() );

            DynamicBitset tVisited( aGraph.size() );

            // Pre-allocate queue with reasonable size (typically ~sqrt(V) for sparse graphs)
            Queue< Vertex * > tQueue;

            // Reusable container for unvisited neighbors - sized for maximum degree
            Graph tUnvisitedNeighbors;
            tUnvisitedNeighbors.reserve( tMaxDegree );

            // Start Cuthill-McKee from the starting vertex
            tQueue.push( aStart );
            tVisited.set( aStart->index() );

            // Cuthill-McKee traversal
            while( !tQueue.empty() )
            {
                Vertex * tCurrent = tQueue.pop();
                tPermutation.push( tCurrent->index() );

                // Clear and reuse the neighbors container
                tUnvisitedNeighbors.clear();

                for( uint k = 0; k < tCurrent->number_of_vertices(); ++k )
                {
                    Vertex * tNeighbor = tCurrent->vertex( k );
                    if( !tVisited.test( tNeighbor->index() ) )
                    {
                        tUnvisitedNeighbors.push( tNeighbor );
                        tVisited.set( tNeighbor->index() );
                    }
                }

                // Sort neighbors by degree (ascending order for better bandwidth reduction)
                if( tUnvisitedNeighbors.size() > 1 )
                {
                    sort( tUnvisitedNeighbors, opVertexDegree );
                }

                // Add sorted neighbors to queue
                for( index_t i = 0; i < tUnvisitedNeighbors.size(); ++i )
                {
                    tQueue.push( tUnvisitedNeighbors( i ) );
                }
            }
            // Handle disconnected components
            if( tPermutation.size() < aGraph.size() )
            {
                // Find unvisited vertices and process each component
                for( index_t i = 0; i < aGraph.size(); ++i )
                {
                    if( !tVisited.test( i ) )
                    {
                        // Start new component with minimum degree vertex
                        Vertex * tComponentStart = aGraph( i );
                        uint tMinComponentDegree = tComponentStart->number_of_vertices();  // Using stored degree

                        // Quick scan for minimum degree vertex in this component
                        for( index_t j = i + 1; j < aGraph.size(); ++j )
                        {
                            if( !tVisited.test( j ) && aGraph( j )->number_of_vertices() < tMinComponentDegree )
                            {
                                tComponentStart = aGraph( j );
                                tMinComponentDegree = tComponentStart->number_of_vertices();
                            }
                        }

                        tQueue.push( tComponentStart );
                        tVisited.set( tComponentStart->index() );

                        // Process component with Cuthill-McKee
                        while( !tQueue.empty() )
                        {
                            Vertex * tCurrent = tQueue.pop();
                            tPermutation.push( tCurrent->index() );

                            // Clear and reuse the neighbors container
                            tUnvisitedNeighbors.clear();

                            for( uint k = 0; k < tCurrent->number_of_vertices(); ++k )
                            {
                                Vertex * tNeighbor = tCurrent->vertex( k );
                                if( !tVisited.test( tNeighbor->index() ) )
                                {
                                    tUnvisitedNeighbors.push( tNeighbor );
                                    tVisited.set( tNeighbor->index() );
                                }
                            }

                            // Sort neighbors by degree
                            if( tUnvisitedNeighbors.size() > 1 )
                            {
                                sort( tUnvisitedNeighbors, opVertexDegree );
                            }

                            // Add sorted neighbors to queue
                            for( index_t k = 0; k < tUnvisitedNeighbors.size(); ++k )
                            {
                                tQueue.push( tUnvisitedNeighbors( k ) );
                            }
                        }
                    }
                }
            }

            // Reverse the permutation for RCM
            reverse( tPermutation );

            // Apply the RCM permutation to the graph
            // First, update vertex indices according to RCM order
            tCount = 0 ;
            for ( index_t k : tPermutation )
            {
                aGraph( k )->set_index( tCount++ );
            }

            // Now sort the graph by the new indices
            sort( aGraph, opVertexIndex );
        }

        /**
         * Applies the RCM permutation to reorder the graph vertices.
         * Updates the index of each vertex according to the new ordering.
         *
         * @param aGraph        Cell containing all vertices in the graph
         * @param aPermutation  The permutation vector from symrcm
         */
        void apply_rcm_permutation( Graph & aGraph, const Cell< index_t > & aPermutation )
        {
            // Create inverse permutation for updating indices
            Cell< index_t > tInversePermutation( aPermutation.size(), 0 );
            for( index_t i = 0; i < aPermutation.size(); ++i )
            {
                tInversePermutation( aPermutation( i ) ) = i;
            }

            // Update vertex indices
            for( index_t i = 0; i < aGraph.size(); ++i )
            {
                aGraph( i )->set_index( tInversePermutation( i ) );
            }

            // Sort the graph array by the new indices using std::sort
            sort( aGraph, opVertexIndex );
        }
    }

//------------------------------------------------------------------------------
}
