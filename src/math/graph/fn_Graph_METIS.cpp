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

#include "typedefs.hpp"
#include "assert.hpp"
#include "graph_typedefs.hpp"
#include "graphtools.hpp"
#include "fn_Graph_METIS.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "cl_Cell.hpp"
#include "fn_to_vector.hpp"


namespace belfem
{
    namespace graph
    {

//------------------------------------------------------------------------------

        void
        metis_nd( Graph & aGraph )
        {
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return;
            }

#ifdef BELFEM_METIS

            metis_t tNumVertices = static_cast< metis_t >( aGraph.size() );

            message( InfoLevel::Verbose, "    Starting METIS NodeND ..." );

            // start the timer
            Timer tTimer;

            // Build adjacency structure
            Vector< metis_t > tVertices ;
            Vector< metis_t > tEdges ;
            build_graph_adjacency( aGraph, tVertices, tEdges );

            // Allocate output arrays
            Vector< metis_t > tPerm( tNumVertices, 0 );
            Vector< metis_t > tIPerm( tNumVertices, 0 );

            // Set METIS options
            Cell< metis_t > tOptions( METIS_NOPTIONS, 0 );
            METIS_SetDefaultOptions( tOptions.data() );

            tOptions( METIS_OPTION_NUMBERING ) = 0;
            tOptions( METIS_OPTION_COMPRESS )  = 0;
            tOptions( METIS_OPTION_CTYPE )     = METIS_CTYPE_SHEM;
            tOptions( METIS_OPTION_RTYPE )     = METIS_RTYPE_SEP1SIDED;

            if( gLog.info_level() >= static_cast< uint >( InfoLevel::Everything ) )
            {
                tOptions( METIS_OPTION_DBGLVL ) = METIS_DBG_INFO;
            }

            // Call METIS nested dissection
            int tStatus = METIS_NodeND(
                &tNumVertices,  // Number of vertices
                tVertices.data(),          // Row pointers
                tEdges.data(),        // Column indices
                nullptr,                // Vertex weights (nullptr = unweighted)
                tOptions.data(),       // Options
                tPerm.data(),          // Permutation array (output): new_index -> old_index
                tIPerm.data()          // Inverse permutation (output): old_index -> new_index
            );

            BELFEM_ERROR( tStatus == METIS_OK,
                         "METIS_NodeND failed with status %d : %s",
                         ( unsigned int ) tStatus,
                         metis_status( tStatus ).c_str() );

            message( InfoLevel::Verbose, "    ... time for reordering             : %u ms",
                     ( unsigned int ) tTimer.stop() );

            // Apply permutation to graph
            apply_graph_permutation( aGraph, tIPerm );

#else
          BELFEM_ERROR( false, "trying to call METIS_NodeND, but we are not linked against METIS!" );
#endif
        }

//------------------------------------------------------------------------------

        void
        metis_ndp( Graph & aGraph, const uint aNumPartitions )
        {
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return;
            }

            // METIS_NodeNDP requires at least 2 partitions
            // For single partition, use standard nested dissection instead
            /*if( aNumPartitions < 2 )
            {
                metis_nd( aGraph );
                return;
            }*/

#ifdef BELFEM_METIS

            metis_t tNumVertices = static_cast< metis_t >( aGraph.size() );
            metis_t tNumPartitions = static_cast< metis_t >( aNumPartitions );

            message( InfoLevel::Verbose, "    Starting METIS NodeNDP ..." );

            // start the timer
            Timer tTimer;

            // Build adjacency structure
            Vector< metis_t > tVertices ;
            Vector< metis_t > tEdges ;
            build_graph_adjacency( aGraph, tVertices, tEdges );

            // Allocate output arrays
            Vector< metis_t > tPerm( tNumVertices, 0 );
            Vector< metis_t > tIPerm( tNumVertices, 0 );
            
            // Sizes array for NodeNDP: contains info about separator sizes
            // Size is 2*npes - 1 according to METIS documentation
            Cell< metis_t > tSizes( 2 * tNumPartitions - 1, 0);

            // Set METIS options
            Cell< metis_t > tOptions( METIS_NOPTIONS, 0 );
            METIS_SetDefaultOptions( tOptions.data() );

            tOptions( METIS_OPTION_NUMBERING ) = 0;
            tOptions( METIS_OPTION_COMPRESS )  = 0;
            tOptions( METIS_OPTION_CTYPE )     = METIS_CTYPE_SHEM;
            tOptions( METIS_OPTION_RTYPE )     = METIS_RTYPE_SEP1SIDED;

            if( gLog.info_level() >= static_cast< uint >( InfoLevel::Everything ) )
            {
                tOptions( METIS_OPTION_DBGLVL ) = METIS_DBG_INFO;
            }

            // Call METIS nested dissection with partitioning
            int tStatus = METIS_NodeNDP(
                tNumVertices,          // Number of vertices
                tVertices.data(),      // Row pointers
                tEdges.data(),         // Column indices
                nullptr,               // Vertex weights (nullptr = unweighted)
                tNumPartitions,        // Number of partitions
                tOptions.data(),       // Options
                tPerm.data(),          // Permutation array (output): new_index -> old_index
                tIPerm.data(),         // Inverse permutation (output): old_index -> new_index
                tSizes.data()          // Separator sizes (output)
            );

            // interpret status from METIS
            BELFEM_ERROR( tStatus == METIS_OK,
                         "METIS_NodeND failed with status %d : %s",
                         ( unsigned int ) tStatus,
                         metis_status( tStatus ).c_str() );

            // Apply permutation to graph
            apply_graph_permutation( aGraph, tIPerm );

            message( InfoLevel::Verbose, "    ... time for reordering             : %u ms",
            ( unsigned int ) tTimer.stop() );

#else
            BELFEM_ERROR( false, "trying to call METIS_NodeND, but we are not linked against METIS!" );
#endif
        }

//------------------------------------------------------------------------------

        void
        metis_partition(
            Graph & aGraph,
            const uint aNumPartitions,
            const bool aForceContinuousPartitions,
            Vector< proc_t > * aPartitions )
        {
#ifdef BELFEM_METIS
            // Handle empty graph
            if( aGraph.size() == 0 )
            {
                return;
            }
            if ( aNumPartitions < 2 )
            {
                if ( aPartitions == nullptr )
                {
                    for ( Vertex * tVertex : aGraph )
                    {
                        tVertex->set_owner( 0 );
                    }
                    return;
                }
                else
                {
                    aPartitions->set_size( aGraph.size(), 0 );
                }
            }

            metis_t tNumVertices = static_cast< metis_t >( aGraph.size() );
            metis_t tNumPartitions = static_cast< metis_t >( aNumPartitions );
            metis_t tNumberOfConstraints = 1;
            message( InfoLevel::Verbose, "    Starting METIS METIS_PartGraphKway ..." );

            // start the timer
            Timer tTimer;

            // Build adjacency structure
            Vector< metis_t > tVertices ;
            Vector< metis_t > tEdges ;
            build_graph_adjacency( aGraph, tVertices, tEdges );

            // Set METIS options
            Cell< metis_t > tOptions( METIS_NOPTIONS, 0 );
            METIS_SetDefaultOptions( tOptions.data() );

            // graph is zero based
            tOptions( METIS_OPTION_NUMBERING ) = 0;

            // no compression of the graph
            tOptions( METIS_OPTION_COMPRESS )  = 0;

            // force contiguous partitions if requested (default: yes)
            tOptions( METIS_OPTION_CONTIG )    = static_cast< metis_t >( aForceContinuousPartitions );

            tOptions( METIS_OPTION_CTYPE )     = METIS_CTYPE_SHEM;
            tOptions( METIS_OPTION_RTYPE )     = METIS_RTYPE_SEP1SIDED;
            tOptions( METIS_OPTION_OBJTYPE )   = METIS_OBJTYPE_VOL;

            // verbisity flag
            if( gLog.info_level() >= static_cast< uint >( InfoLevel::Everything ) )
            {
                tOptions( METIS_OPTION_DBGLVL ) = METIS_DBG_INFO;
            }

            Cell< metis_t > tPart( tNumVertices, 0 );

            // the result
            metis_t tEdgeCut;

            int tStatus =
                   METIS_PartGraphKway(
                           &tNumVertices,            // The number of vertices in the graph.
                           &tNumberOfConstraints,  // The number of balancing constraints. It should be at least 1.
                           tVertices.data(),       // The adjacency structure of the graph
                           tEdges.data(),             // The adjacency structure of the graph
                           nullptr,                   // The weights of the vertices
                           nullptr,                   // The size of the vertices for computing the total communication volume
                           nullptr,                   // The weights of the edges
                           &tNumPartitions,           // The number of parts to partition the graph.
                           nullptr,                   // Target Partition weights
                           nullptr,                   // the allowed load imbalance tolerance for each constraint
                           tOptions.data(),               // Options for METIS
                           &tEdgeCut,
                           tPart.data() );


            // interpret status from METIS
            BELFEM_ERROR( tStatus == METIS_OK,
                         "METIS_PartGraphKway failed with status %d : %s",
                         ( unsigned int ) tStatus,
                         metis_status( tStatus ).c_str() );


            if ( aPartitions == nullptr )
            {
                // Apply partitioning to the graph
                index_t tCount = 0;

                for ( Vertex * tVertex : aGraph )
                {
                    tVertex->set_owner( tPart( tCount++ ) );
                }
            }
            else
            {
                // copy partitions to the output
                *aPartitions = to_vector( tPart );
            }
            message( InfoLevel::Verbose, "    ... time for partitioning           : %u ms",
            ( unsigned int ) tTimer.stop() );

#ifdef DEBUG
            Vector< index_t > tCounters( aNumPartitions, 0 );
            for ( Vertex * tVertex : aGraph )
            {
                ++tCounters( tVertex->owner() );
            }
#endif

#else
            BELFEM_ERROR( false, "trying to call METIS_PartGraphKway, but we are not linked against METIS!" );
#endif


        }
//------------------------------------------------------------------------------
        string
        metis_status( const int aStatus )
        {
#ifdef BELFEM_METIS
            switch ( aStatus )
            {
                case( METIS_OK ) :
                {
                    return "Success";
                }
                case( METIS_ERROR_INPUT ) :
                {
                    return "Invalid input";
                }
                case( METIS_ERROR_MEMORY ) :
                {
                    return "Memory Error";
                }
                default :
                {
                    return "Unknown Error";
                }
            }
#else
            return "We are not linked against METIS";
#endif
        }

//------------------------------------------------------------------------------
    }
}
