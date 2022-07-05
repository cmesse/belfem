//
// Created by Christian Messe on 2019-01-20.
//
#include "assert.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "fn_Graph_find_pseudo_peripheral_node.hpp"
#include "fn_Graph_BFS.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        Vertex *
        find_pseudo_peripheral_node(
                Cell< Vertex * > & aGraph,
                Vertex *     aStart  )
        {

            // start timer
            //Timer tTimer;

            // width of graph
            uint tWidth;

            // perform a breath width search
            graph::bfs( aGraph, tWidth, aStart );

            // loop counter
            index_t tCount = 0;

            while( true )
            {
                // remember last width
                uint tOldWidth = tWidth;

                // reset max degree = 0;
                uint tMaxDegree = 0;

                // count nodes with max level
                for ( Vertex * tVertex : aGraph )
                {

                    if ( tVertex->level() == tWidth )
                    {
                        // increment counter
                        tMaxDegree = std::max( tMaxDegree, tVertex->number_of_vertices() );
                    }
                }

                // find any node with max degree
                for ( graph::Vertex * tVertex : aGraph )
                {
                    if ( tVertex->level() == tWidth )
                    {
                        if ( tVertex->number_of_vertices() == tMaxDegree )
                        {
                            aStart = tVertex;
                            break;
                        }
                    }
                }

                // perform a breath width search
                graph::bfs( aGraph, tWidth, aStart );

                // break the while loop
                if( tWidth <= tOldWidth )
                {
                    break;
                }
                else
                {
                    // increment loop counter
                    tCount++;

                    BELFEM_ERROR( tCount < aGraph.size(),
                               "find_pseudo_peripheral_node failed" );

                }
            }

            //message( 4, "    Time for finding pseudo peripheral node : %u ms\n",
            //         ( unsigned int ) tTimer.stop() );

            // return node
            return aStart;
        }
    }
}