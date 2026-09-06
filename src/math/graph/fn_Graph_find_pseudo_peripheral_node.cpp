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

#include "assert.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "fn_Graph_find_pseudo_peripheral_node.hpp"
#include "fn_Graph_bfs.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        Vertex *
        find_pseudo_peripheral_node(
                Graph & aGraph,
                Vertex *     aStart  )
        {
            // perform a breadth-first search
            graph::bfs( aGraph, aStart );

            // compute max level (eccentricity) — NOT the BFS return value,
            // which is the max width (largest number of vertices at any level)
            index_t tMaxLevel = 0;
            for ( Vertex * tVertex : aGraph )
            {
                if ( tVertex->level() != gNoIndex )
                {
                    tMaxLevel = std::max( tMaxLevel, tVertex->level() );
                }
            }

            // loop counter
            index_t tCount = 0;

            while( true )
            {
                // remember last max level
                index_t tOldMaxLevel = tMaxLevel;

                // reset max degree
                uint tMaxDegree = 0;

                // find max degree among vertices at the farthest level
                for ( Vertex * tVertex : aGraph )
                {
                    if ( tVertex->level() == tMaxLevel )
                    {
                        tMaxDegree = std::max( tMaxDegree, tVertex->number_of_vertices() );
                    }
                }

                // find any vertex with max degree at the farthest level
                for ( graph::Vertex * tVertex : aGraph )
                {
                    if ( tVertex->level() == tMaxLevel )
                    {
                        if ( tVertex->number_of_vertices() == tMaxDegree )
                        {
                            aStart = tVertex;
                            break;
                        }
                    }
                }

                // perform BFS from the new candidate
                graph::bfs( aGraph, aStart );

                // recompute max level
                tMaxLevel = 0;
                for ( Vertex * tVertex : aGraph )
                {
                    if ( tVertex->level() != gNoIndex )
                    {
                        tMaxLevel = std::max( tMaxLevel, tVertex->level() );
                    }
                }

                // break if eccentricity did not increase
                if( tMaxLevel <= tOldMaxLevel )
                {
                    break;
                }
                else
                {
                    tCount++;
                    BELFEM_ERROR( tCount < aGraph.size(),
                               "find_pseudo_peripheral_node failed" );
                }
            }

            return aStart;
        }
    }
}