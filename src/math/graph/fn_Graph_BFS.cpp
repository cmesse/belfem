//
// Created by Christian Messe on 2019-01-20.
//

#include "fn_Graph_BFS.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        void
        bfs(    Cell< Vertex * > & aGraph,
                   uint & aWidth,
                   Vertex * aStart )
        {
            // reset all node flags and levels
            for( Vertex * tVertex : aGraph )
            {
                tVertex->unflag();
                tVertex->level() = 0;
            }

            // get start if not prescribed
            if( aStart == nullptr )
            {
                aStart = aGraph( 0 );
            }

            // Queue
            Cell< Vertex * > tQueue;

            // put start into Queue
            tQueue.push( aStart );

            // init level of start node
            aStart->level() = 1;

            // reset width
            aWidth = 0;

            while( tQueue.size() > 0 )
            {
                // get node from queue
                Vertex * tVertex = tQueue.pop();

                // compare level with width
                aWidth = std::max( tVertex->level(), aWidth );

                // get number of neighbors
                index_t tNumberOfNeighbors = tVertex->number_of_vertices();

                // loop over all neighbors of this node
                for( index_t k=0; k<tNumberOfNeighbors; ++k )
                {
                    // get pointer to neighbor
                    Vertex * tNeighbor = tVertex->vertex( k );

                    // test if node was visited
                    if( ! tNeighbor->is_flagged() )
                    {
                        // put neighbor on queue
                        tQueue.push( tNeighbor );

                        // flag neighbor
                        tNeighbor->flag();

                        // set level of this node
                        tNeighbor->level() = tVertex->level() + 1;
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    }
}
