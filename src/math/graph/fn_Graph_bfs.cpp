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

#include "fn_Graph_bfs.hpp"
#include "cl_Queue.hpp"
#include "cl_DynamicBitset.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        index_t
        bfs( Graph & aGraph, Vertex * aStart = nullptr )
        {
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return 0;
            }

            // Use first vertex if no start vertex is provided
            if( aStart == nullptr )
            {
                aStart = aGraph( 0 );
            }

            index_t tCount = 0 ;
            for ( Vertex * tVertex: aGraph )
            {
                tVertex->set_index( tCount++ );
                tVertex->set_level( gNoIndex );
            }
            aStart->set_level( 0 );

            // Initialize visited bitset
            DynamicBitset tVisited( aGraph.size() );

            // Queues for current and next level
            Queue< Vertex * > tCurrentLevel;
            Queue< Vertex * > tNextLevel;

            // Start BFS
            tCurrentLevel.push( aStart );
            tVisited.set( aStart->index() );

            // Track maximum width
            index_t aMaxWidth = 1;  // At least one vertex (the start)

            // BFS traversal
            while( !tCurrentLevel.empty() )
            {
                // Process all vertices at current level
                while( !tCurrentLevel.empty() )
                {
                    Vertex * tVertex = tCurrentLevel.pop();

                    // Visit all neighbors
                    for( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                    {
                        Vertex * tNeighbor = tVertex->vertex( k );

                        // Check if neighbor has been visited
                        if( !tVisited.test( tNeighbor->index() ) )
                        {
                            tVisited.set( tNeighbor->index() );
                            tNeighbor->set_level( tVertex->level() + 1 );
                            tNextLevel.push( tNeighbor );
                        }
                    }
                }

                // Update maximum width if next level has more vertices
                if( tNextLevel.size() > aMaxWidth )
                {
                    aMaxWidth = tNextLevel.size();
                }

                // Move to next level
                std::swap( tCurrentLevel, tNextLevel );
            }

            return aMaxWidth;
        }

        /**
         * Alternative implementation that also handles disconnected graphs
         * by finding the maximum width across all connected components.
         *
         * @param aGraph     Cell containing all vertices in the graph
         * @return           The maximum width across all connected components
         */
        index_t bfs( Graph & aGraph )
        {
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return 0;
            }

            // Initialize visited bitset for entire graph
            DynamicBitset tGlobalVisited( aGraph.size() );

            index_t aMaxWidth = 0;

            // Process each connected component
            for( index_t i = 0; i < aGraph.size(); ++i )
            {
                if( !tGlobalVisited.test( i ) )
                {
                    // Found unvisited vertex - new component
                    Vertex * tStart = aGraph( i );

                    tStart->set_level( 0 );

                    // Queues for current and next level
                    Queue< Vertex * > tCurrentLevel;
                    Queue< Vertex * > tNextLevel;

                    // Start BFS for this component
                    tCurrentLevel.push( tStart );
                    tGlobalVisited.set( tStart->index() );

                    // Track maximum width for this component
                    index_t tComponentMaxWidth = 1;

                    // BFS traversal
                    while( !tCurrentLevel.empty() )
                    {
                        // Process all vertices at current level
                        while( !tCurrentLevel.empty() )
                        {
                            Vertex * tVertex = tCurrentLevel.pop();

                            // Visit all neighbors
                            for( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                            {
                                Vertex * tNeighbor = tVertex->vertex( k );

                                // Check if neighbor has been visited
                                if( !tGlobalVisited.test( tNeighbor->index() ) )
                                {
                                    tGlobalVisited.set( tNeighbor->index() );
                                    tNeighbor->set_level( tVertex->level() + 1 );
                                    tNextLevel.push( tNeighbor );
                                }
                            }
                        }

                        // Update component maximum width if next level has more vertices
                        if( tNextLevel.size() > tComponentMaxWidth )
                        {
                            tComponentMaxWidth = tNextLevel.size();
                        }

                        // Move to next level
                        std::swap( tCurrentLevel, tNextLevel );
                    }

                    // Update global maximum width
                    if( tComponentMaxWidth > aMaxWidth )
                    {
                        aMaxWidth = tComponentMaxWidth;
                    }
                }
            }

            return aMaxWidth;
        }

//------------------------------------------------------------------------------
    }
}
