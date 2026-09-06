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

#include "cl_Logger.hpp"
#include "fn_Graph_find_pseudo_peripheral_vertex.hpp"
#include "cl_Queue.hpp"
#include "cl_DynamicBitset.hpp"
namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        Vertex *
        find_pseudo_peripheral_vertex( Graph & aGraph, Vertex * aStart )
        {
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return nullptr;
            }

            // Use first vertex if no start vertex is provided
            if( aStart == nullptr )
            {
                aStart = aGraph( 0 );
            }

            Vertex * aCurrent = aStart;
            index_t tPreviousEccentricity = 0;

            // Iterate until we can't find a vertex with larger eccentricity
            while( true )
            {
                // Perform BFS from current vertex to find farthest vertices
                DynamicBitset tVisited( aGraph.size() );
                Queue< Vertex * > tQueue;

                // Start BFS
                tQueue.push( aCurrent );
                tVisited.set( aCurrent->index() );

                // Track vertices at the farthest level
                Graph tFarthestVertices;
                index_t tCurrentLevel = 0;

                // BFS traversal tracking levels
                while( !tQueue.empty() )
                {
                    // Clear previous farthest vertices as we found a new level
                    tFarthestVertices.clear();

                    // Get size of current level
                    size_t tLevelSize = tQueue.size();

                    // Process all vertices at current level
                    for( size_t i = 0; i < tLevelSize; ++i )
                    {
                        Vertex * tVertex = tQueue.pop();
                        tFarthestVertices.push( tVertex );

                        // Visit all neighbors
                        for( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                        {
                            Vertex * tNeighbor = tVertex->vertex( k );

                            // Check if neighbor has been visited
                            if( !tVisited.test( tNeighbor->index() ) )
                            {
                                tVisited.set( tNeighbor->index() );
                                tQueue.push( tNeighbor );
                            }
                        }
                    }

                    if( !tQueue.empty() )
                    {
                        ++tCurrentLevel;
                    }
                }

                // Check if we've improved the eccentricity
                if( tCurrentLevel <= tPreviousEccentricity )
                {
                    // No improvement, we've found a pseudo-peripheral vertex
                    return aCurrent;
                }

                // Update eccentricity
                tPreviousEccentricity = tCurrentLevel;

                // Among the farthest vertices, choose one with minimum degree
                // This heuristic often leads to better pseudo-peripheral vertices
                Vertex * tNextCandidate = tFarthestVertices( 0 );
                uint tMinDegree = tNextCandidate->number_of_vertices();

                for( index_t i = 1; i < tFarthestVertices.size(); ++i )
                {
                    uint tDegree = tFarthestVertices( i )->number_of_vertices();
                    if( tDegree < tMinDegree )
                    {
                        tMinDegree = tDegree;
                        tNextCandidate = tFarthestVertices( i );
                    }
                }

                aCurrent = tNextCandidate;
            }
        }

    }
}