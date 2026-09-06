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
#include "graph_typedefs.hpp"
#include "graphtools.hpp"
#include "fn_Graph_PTSCOTCH.hpp"

#ifdef BELFEM_PTSCOTCH
#include <ptscotch.h>
#endif
#include "commtools.hpp"
#include "fn_Graph_SCOTCH.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"

namespace belfem
{
    namespace graph
    {
        void
        ptscotch_nd( Graph & aGraph )
        {
#ifdef BELFEM_PTSCOTCH
            if ( comm_size() < 2 )
            {
              scotch_nd( aGraph );
              return;
            }

            Vector< scotch_t > tDistribution;
            Vector< scotch_t > tVertices;
            Vector< scotch_t > tEdges;

            Timer tTimer;

            proc_t tCommRank = comm_rank();
            proc_t tCommSize = comm_size();

            scotch_t tNumVertices = 0;

            if ( tCommRank == 0 )
            {
              BELFEM_ERROR( aGraph.size() > 0, "Graph is empty" );

              Cell< Vector< scotch_t > > tAllVertices ;
              Cell< Vector< scotch_t > > tAllEdges;

              build_pargraph_adjacency(
                  aGraph, tDistribution,
                  tAllVertices,
                  tAllEdges );

              tNumVertices = static_cast< scotch_t >( aGraph.size() );
              comm_barrier();
              broadcast( tNumVertices );
              share( tDistribution );
              distribute( tAllVertices );
              distribute( tAllEdges );
              tVertices = std::move( tAllVertices( 0 ) );
              tEdges = std::move( tAllEdges( 0 ) );
            }
            else
            {
              comm_barrier();
              broadcast( tNumVertices );
              receive( tDistribution );
              receive( tVertices );
              receive( tEdges );
            }

            // Local number of vertices
            scotch_t tLocalNumVerts = static_cast< scotch_t >( tVertices.length() - 1 );

            // Global vertex labels (based on distribution). Sized at least one
            // so vlblloctab is never null on an empty rank: Armadillo hands out
            // a null memptr for a zero-length vector ( Blaze does not ), and
            // SCOTCH_dgraphBuild requires every rank to agree on whether labels
            // are supplied -- same reason the edge buffer keeps a placeholder.
            Vector< scotch_t > tLocalVertexIndices( tLocalNumVerts > 0 ? tLocalNumVerts : 1, 0 );
            scotch_t tGlobalOffset = static_cast< scotch_t >( tDistribution( tCommRank ) );
            for ( scotch_t i = 0; i < tLocalNumVerts; ++i )
            {
              tLocalVertexIndices( i ) = tGlobalOffset + i;
            }

            // the logical edge count is the CSR terminal, not the buffer length:
            // build_pargraph_adjacency() hands an empty rank a one-element
            // placeholder buffer so the pointer is never null
            scotch_t tLocalNumEdges = tVertices( tLocalNumVerts );
            scotch_t tEdgeBufferSize = static_cast< scotch_t >( tEdges.length() );

            SCOTCH_Dgraph tDGraph;
            SCOTCH_dgraphInit( &tDGraph, gComm.world() );

            comm_barrier();

            // grafptr, baseval, vertlocnbr, vertlocmax, vertloctab, vendloctab,
            // veloloctab, vlblloctab, edgelocnbr, edgelocsiz, edgeloctab, edgegsttab, edloloctab
            int tStatus = SCOTCH_dgraphBuild(
              &tDGraph,                      // Distributed graph structure
              0,                             // baseval (0-based indexing)
              tLocalNumVerts,                // vertlocnbr (number of local vertices)
              tLocalNumVerts,                // vertlocmax (max local vertices)
              tVertices.data(),              // vertloctab (CSR vertex start array)
              nullptr,                       // vendloctab (vertex end array, NULL means vertloctab+1)
              nullptr,                       // veloloctab (vertex weights, NULL = no weights)
              tLocalVertexIndices.data(),    // vlblloctab (global vertex labels)
              tLocalNumEdges,                // edgelocnbr (number of local edges)
              tEdgeBufferSize,               // edgelocsiz (size of edge array)
              tEdges.data(),                 // edgeloctab (edge array)
              nullptr,                       // edgegsttab (ghost edges, NULL = no ghosts)
              nullptr );                     // edloloctab (edge weights, NULL = no weights)

            BELFEM_ERROR( tStatus == 0, "SCOTCH_dgraphBuild failed with status %d", tStatus );

            SCOTCH_Strat tStrat;
            SCOTCH_stratInit( &tStrat );

            // Strategy for distributed ND
            tStatus = SCOTCH_stratDgraphOrder( &tStrat, "n" );
            BELFEM_ERROR( tStatus == 0, "SCOTCH_stratDgraphOrder failed with status %d", tStatus );

            // Type is SCOTCH_Dordering, not SCOTCH_Dorder
            SCOTCH_Dordering tDOrder;
            tStatus = SCOTCH_dgraphOrderInit( &tDGraph, &tDOrder );
            BELFEM_ERROR( tStatus == 0, "SCOTCH_dgraphOrderInit failed with status %d", tStatus );

            tStatus = SCOTCH_dgraphOrderCompute( &tDGraph, &tDOrder, &tStrat );
            BELFEM_ERROR( tStatus == 0, "SCOTCH_dgraphOrderCompute failed with status %d", tStatus );

            // Get local permutation
            Vector< scotch_t > tLocalPerm( tLocalNumVerts );
            tStatus = SCOTCH_dgraphOrderPerm( &tDGraph, &tDOrder, tLocalPerm.data() );
            BELFEM_ERROR( tStatus == 0, "SCOTCH_dgraphOrderPerm failed with status %d", tStatus );

            if ( tCommRank == 0 )
            {
              message( InfoLevel::Verbose, "    ... time for reordering             : %u ms",
                       ( unsigned int ) tTimer.stop() );

              Cell< Vector< scotch_t > > tAllPerms( tCommSize, {} );
              comm_barrier();
              collect( tAllPerms );
              tAllPerms( 0 ) = std::move( tLocalPerm );

              Vector< scotch_t > tGlobalPerm( tNumVertices );

              Vector< index_t > tCount( tCommSize, 0 );
              for ( Vertex * tVertex : aGraph )
              {
                  proc_t tOwner = tVertex->owner();
                  tGlobalPerm( tVertex->index() ) = tAllPerms( tOwner )( tCount( tOwner )++ );
              }

              apply_graph_permutation( aGraph, tGlobalPerm );
            }
            else
            {
              comm_barrier();
              send( tLocalPerm );
            }

            // Cleanup
            SCOTCH_dgraphOrderExit( &tDGraph, &tDOrder );
            SCOTCH_stratExit( &tStrat );
            SCOTCH_dgraphExit( &tDGraph );

#else
            BELFEM_ERROR( false, "trying to call PTSCOTCH_dgraphOrderCompute, but we are not linked against PTSCOTCH!" );
#endif
        }
    }
}