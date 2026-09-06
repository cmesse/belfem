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

#include <stack>
#include <tuple>

#include "fn_Graph_dfs.hpp"

namespace belfem
{
    namespace graph
    {


        proc_t
        dfs( Graph & aGraph )
        {
            // early exit
            if ( aGraph.size() == 0 ) return 0;

            // first, let's make sure that the graph is reset
            for ( Vertex * tVertex: aGraph )
            {
                tVertex->set_owner( gNoOwner );
                tVertex->set_level( gNoIndex );
                tVertex->unflag();
            }

            proc_t aNumSubGraphs = 0 ;

            // next, we loop over all vertices
            for ( Vertex * tVertex: aGraph )
            {
                if ( ! tVertex->is_flagged() )
                {
                    tVertex->set_owner( aNumSubGraphs++ );
                    dfs_from_start( aGraph, tVertex );
                }
            }

            return aNumSubGraphs;
        }

        void
        dfs_from_start(  Graph & aGraph, Vertex * aStart  )
        {
            std::stack< Vertex * >  tStack;
            tStack.push( aStart );
            aStart->set_level( 0 );
            aStart->flag();

            while ( ! tStack.empty() )
            {
                Vertex * tVertex  = tStack.top(); tStack.pop();

                for ( uint k=0; k<tVertex->number_of_vertices(); ++k )
                {
                    Vertex * tNeighbor = tVertex->vertex( k );

                    if ( ! tNeighbor->is_flagged() )
                    {
                        // set visited flag
                        tNeighbor->flag();

                        // set the level of the neigbor
                        tNeighbor->set_level( tVertex->level() + 1 );

                        // set the owner
                        tNeighbor->set_owner( aStart->owner() );

                        // add neighbor to stack
                        tStack.push( tNeighbor );
                    }
                }
            }
        }

    }
}
