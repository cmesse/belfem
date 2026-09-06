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
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "graph_typedefs.hpp"
#include "graphtools.hpp"
#include "fn_Graph_SCOTCH.hpp"




namespace belfem
{
    namespace graph
    {
        void
        scotch_nd( Graph & aGraph )
        {

#ifdef BELFEM_SCOTCH
             // Handle empty graph
              if ( aGraph.size() == 0 )
              {
                  return;
              }

              scotch_t tNumVertices = static_cast< scotch_t >( aGraph.size() );

              message( InfoLevel::Verbose, "    Starting SCOTCH Nested Dissection ..." );

              // Start the timer
              Timer tTimer;

              // Build adjacency structure (reuses METIS format)
              Vector< scotch_t > tVertices;
              Vector< scotch_t > tEdges;
              build_graph_adjacency( aGraph, tVertices, tEdges );

              // logical edge count from the CSR terminal, not the buffer: the
              // builder keeps a one-element placeholder when there are no edges
              scotch_t tNumEdges = tVertices( tNumVertices );

              SCOTCH_Graph tGraph;
              SCOTCH_graphInit( &tGraph );

              scotch_t tBaseVal = 0;  // 0-based indexing

              int tStatus = SCOTCH_graphBuild( &tGraph, tBaseVal, tNumVertices,
                                               tVertices.data(), nullptr, nullptr, nullptr,
                                               tNumEdges, tEdges.data(), nullptr );

              BELFEM_ERROR( tStatus == 0, "SCOTCH_graphBuild failed with status %d", tStatus );

              SCOTCH_Strat tStrat;
              SCOTCH_stratInit( &tStrat );

              // Set strategy for nested dissection (simple "n" for ND; customize as needed)
              tStatus = SCOTCH_stratGraphOrder( &tStrat, "n" );
              BELFEM_ERROR( tStatus == 0, "SCOTCH_stratGraphOrder failed with status %d", tStatus );

              // Allocate permutation array that SCOTCH will fill
              Vector< scotch_t > tPerm( tNumVertices );

              // Initialize ordering structure with our permutation array
              SCOTCH_Ordering tOrder;
              tStatus = SCOTCH_graphOrderInit( &tGraph, &tOrder, tPerm.data(), nullptr,
                                               nullptr, nullptr, nullptr );
              BELFEM_ERROR( tStatus == 0, "SCOTCH_graphOrderInit failed with status %d", tStatus );

              // Compute the ordering - this fills tPerm automatically
              tStatus = SCOTCH_graphOrderCompute( &tGraph, &tOrder, &tStrat );
              BELFEM_ERROR( tStatus == 0, "SCOTCH_graphOrderCompute failed with status %d", tStatus );

              // Apply the permutation (tPerm is now filled by SCOTCH)
              apply_graph_permutation( aGraph, tPerm );

              // Clean up
              SCOTCH_graphOrderExit( &tGraph, &tOrder );
              SCOTCH_stratExit( &tStrat );
              SCOTCH_graphExit( &tGraph );

              message( InfoLevel::Verbose, "    ... time for reordering             : %u ms",
                       ( unsigned int ) tTimer.stop() );

#else
            BELFEM_ERROR( false, "trying to call SCOTCH_graphOrderCompute, but we are not linked against SCOTCH!" );
#endif
        }

        void
        scotch_partition( Graph & aGraph, const uint aNumPartitions )
        {
#ifdef BELFEM_SCOTCH
            if ( aGraph.size() == 0 )
            {
              return;
            }

            if ( aNumPartitions < 2 )
            {
              for ( Vertex * tVertex : aGraph )
              {
                  tVertex->set_owner( 0 );
              }
              return;
            }

            scotch_t tNumVertices = static_cast< scotch_t >( aGraph.size() );
            scotch_t tNumPartitions = static_cast< scotch_t >( aNumPartitions );

            message( InfoLevel::Verbose, "    Starting SCOTCH Partitioning ..." );

            Timer tTimer;

            Vector< scotch_t > tVertices;
            Vector< scotch_t > tEdges;
            build_graph_adjacency( aGraph, tVertices, tEdges );

            // CSR terminal, not buffer length ( placeholder, see scotch_nd )
            scotch_t tNumEdges = tVertices( tNumVertices );

            SCOTCH_Graph tGraph;
            SCOTCH_graphInit( &tGraph );

            scotch_t tBaseVal = 0;

            int tStatus = SCOTCH_graphBuild( &tGraph, tBaseVal, tNumVertices,
                                           tVertices.data(), nullptr, nullptr, nullptr,
                                           tNumEdges, tEdges.data(), nullptr );

            BELFEM_ERROR( tStatus == 0, "SCOTCH_graphBuild failed with status %d", tStatus );

            SCOTCH_Strat tStrat;
            SCOTCH_stratInit( &tStrat );

            // Strategy for mapping/partitioning - use SCOTCH_stratGraphMap
            tStatus = SCOTCH_stratGraphMap( &tStrat, "r" );
            BELFEM_ERROR( tStatus == 0, "SCOTCH_stratGraphMap failed with status %d", tStatus );

            Vector< scotch_t > tPart( tNumVertices );

            tStatus = SCOTCH_graphPart( &tGraph, tNumPartitions, &tStrat, tPart.data() );
            BELFEM_ERROR( tStatus == 0, "SCOTCH_graphPart failed with status %d", tStatus );

            // Apply owners
            index_t tCount = 0;
            for ( Vertex * tVertex : aGraph )
            {
              tVertex->set_owner( static_cast< proc_t >( tPart( tCount++ ) ) );
            }

            // Cleanup
            SCOTCH_stratExit( &tStrat );
            SCOTCH_graphExit( &tGraph );

            message( InfoLevel::Verbose, "    ... time for partitioning           : %u ms",
                   ( unsigned int ) tTimer.stop() );


#else
            BELFEM_ERROR( false, "trying to call SCOTCH_graphPart, but we are not linked against SCOTCH!" );
#endif
        }

    }
}
