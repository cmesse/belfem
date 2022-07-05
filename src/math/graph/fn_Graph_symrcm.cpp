//
// Created by Christian Messe on 2019-01-20.
//
#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "fn_Graph_find_pseudo_peripheral_node.hpp"
#include "fn_Graph_symcrm.hpp"
#include "op_Graph_Vertex_Degree.hpp"
#include "cl_Progressbar.hpp"
namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        void
        symrcm( Cell< Vertex * > & aGraph,
                      Vertex *     aStart,
                const size_t       aNumberOfVertices )
        {
            // start timer
            Timer tTimer;



            // only look for node if no start was passed
            if( aStart == nullptr )
            {
                aStart = aGraph( 0 );
                graph::find_pseudo_peripheral_node( aGraph, aStart );
            }

            index_t tCount = 0;
            for( Vertex * tVertex : aGraph )
            {
                tVertex->unflag();
                tVertex->set_index( tCount++ );
            }

            // number of nodes
            lsint tNumberOfVertices;
            if( aNumberOfVertices == 0 )
            {
                tNumberOfVertices = aGraph.size();
            }
            else
            {
                tNumberOfVertices = aNumberOfVertices;
            }

            message( 2, "\nperforming reverse Cuthill–McKee reordering on %lu vertices ... ",
                     ( long unsigned int ) tNumberOfVertices );

            // the new sorting
            Vector< index_t > tP( tNumberOfVertices, tNumberOfVertices );

            // index
            lsint k = tNumberOfVertices-1;

            // backup index for infinite loop detection
            lsint k_old = tNumberOfVertices;

            // start node is last
            tP( k ) = aStart->index();
            aStart->flag();

            //Vector< uint > tNodesPerStep( tNumberOfVertices );
            //Vector< uint > tDegreePerStep( tNumberOfVertices );

            Progressbar tProgess( tNumberOfVertices );

            tCount = 0;

            while( true )
            {
                // backup index
                k_old = k;

                // reverse loop
                for ( lsint i = tNumberOfVertices - 1; i >= 0; --i )
                {
                    // a node queue
                    Cell< graph::Vertex * > tQueue;
                    for ( lsint j = 0; j < tNumberOfVertices; ++j )
                    {
                        // get pointer to node
                        Vertex * tVertex = aGraph( j );

                        // test if node was visited
                        if ( ! tVertex->is_flagged() )
                        {

                            lsint tNumberOfNeighbors = tVertex->number_of_vertices();

                            // loop over all neighbors of node
                            for ( lsint l = 0; l < tNumberOfNeighbors; ++l )
                            {
                                // get pointer to neighbor
                                Vertex * tNeighbor = tVertex->vertex( l );

                                // test if node is connected to last node on queue
                                if ( tNeighbor->index() == tP( i ) )
                                {
                                    tProgess.step( tCount++ );

                                    // flag this node
                                    tVertex->flag();

                                    // add node to queue
                                    tQueue.push( tVertex );

                                    // break the neighbor loop ( index k )
                                    break;
                                }
                            }
                        }
                    }

                    // sort queue according to degree
                    sort( tQueue, opVertexDegree );

                    // debug parameters
                    uint tOldDegree = 0;

                    // loop over all nodes in queue
                    for ( Vertex * tVertex : tQueue )
                    {
                        // get degree of node
                        uint tDegree = tVertex->number_of_vertices();

                        BELFEM_ERROR ( tDegree >= tOldDegree, "error in sorting" );
                        tOldDegree = tDegree;

                        // add index to new order
                        tP( --k ) = tVertex->index();
                    }


                } // end loop over all nodes

                // test if all nodes have been used
                if( k <= 0 )
                {
                    // break the i-loop
                    break;
                }
                else
                {
                    // detect infinite loop
                    BELFEM_ERROR( k_old != k, "Infinite loop at symcrm function detected. Value of k is %i", (int) k );

                    // find a vertex that is not flagged and connected
                    for( graph::Vertex * tVertex : aGraph )
                    {
                        if( ! tVertex->is_flagged() && tVertex->number_of_vertices() > 0 )
                        {
                            tP( k ) = tVertex->index();

                            // break the for loop and start a new i-loop
                            break;
                        }
                    }
                }

            }

            tProgess.finish();

            // change IDs of connected nodes
            for( lsint k=0; k<tNumberOfVertices; ++k )
            {
                aGraph( tP( k ) )->set_index( k );
            }

            lsint tNumberOfAllVertices = aGraph.size();
            for( lsint k=tNumberOfVertices; k<tNumberOfAllVertices; ++k )
            {
                aGraph( k )->set_index( k );
            }

            Cell< Vertex * > tTempGraph( aGraph );

            aGraph.set_size( tNumberOfVertices, nullptr );

            // rearrange nodes in target array
            for( Vertex * tVertex: tTempGraph )
            {
                aGraph( tVertex->index() ) = tVertex;
            }

            // write output message
            message( 4, "    ....  Time: %i\n",
            tTimer.stop() );
        }

//------------------------------------------------------------------------------
    }
}