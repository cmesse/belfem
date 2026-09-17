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

#include "fn_Graph_multibfs.hpp"


#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Queue.hpp"

namespace belfem
{
    namespace graph
    {
        Vertex *
        multibfs(
            Cell< Vertex * > & aGraph,
            const proc_t aOwnerSubset )
        {
            if ( aGraph.size() == 0 )
            {
                return nullptr ;
            }

            // distance transform: seeds ( flagged vertices ) start at
            // level 1, unvisited at 0, the BFS writes parent + 1 — so
            // level - 1 is the graph distance to the nearest seed
            index_t tCount = 0 ;

            // set levels
            for ( Vertex * tVertex : aGraph )
            {
                if ( tVertex->is_flagged() )
                {
                    ++tCount ;
                    tVertex->set_level( 1 );
                }
                else
                {
                    tVertex->set_level( 0 );
                }
            }

            // collect sources
            Cell< Vertex * > tSources( tCount, nullptr );
            tCount = 0 ;
            for ( Vertex * tVertex : aGraph )
            {
                if ( tVertex->is_flagged() )
                {
                    tSources( tCount++ ) = tVertex ;
                }
            }

            Queue tQueue( tSources );

            if ( aOwnerSubset == gNoOwner )
            {
                while (  ! tQueue.empty() )
                {
                    Vertex * tU = tQueue.pop() ;
                    for ( uint k=0; k<tU->number_of_vertices(); ++k )
                    {
                        Vertex * tV = tU->vertex( k );

                        if ( tV->level() == 0 )
                        {
                            tV->set_level( tU->level() + 1 );
                            tQueue.push( tV );
                        }
                    }
                }
            }
            else
            {
                while (  ! tQueue.empty() )
                {
                    Vertex * tU = tQueue.pop() ;
                    for ( uint k=0; k<tU->number_of_vertices(); ++k )
                    {
                        Vertex * tV = tU->vertex( k );

                        if ( tV->level() == 0 && tV->owner() == aOwnerSubset )
                        {
                            tV->set_level( tU->level() + 1 );
                            tQueue.push( tV );
                        }
                    }
                }
            }

            // running best RESTRICTED to the requested subset: deepest
            // level first, ties broken by smallest id — deterministic
            // across runs and stdlibs, and never a vertex of another
            // component ( the whole-graph sort + last() this replaces
            // could steal a foreign seed when the subset had none ).
            // Semantics of the winner's level for the caller:
            //   0  — the subset holds no seed at all; the tie-break has
            //        returned its smallest id ( an enclosed component )
            //   1  — seeds but no interior: every subset vertex is a seed
            //  >=2 — the regular case, an interior vertex
            Vertex * tResult = nullptr ;
            for ( Vertex * tVertex : aGraph )
            {
                if ( aOwnerSubset != gNoOwner
                     && tVertex->owner() != aOwnerSubset )
                {
                    continue ;
                }

                if (    tResult == nullptr
                     || tVertex->level() > tResult->level()
                     || (    tVertex->level() == tResult->level()
                          && tVertex->id() < tResult->id() ) )
                {
                    tResult = tVertex ;
                }
            }

            return tResult ;
        }
    }
}
