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

#include "cl_Vector.hpp"
#include "fn_Graph_find_connected_partitions.hpp"
#include "op_Graph_Vertex_Owner.hpp"

namespace belfem
{
    namespace graph
    {
        index_t
        find_connected_partitions( Graph & aGraph )
        {
            if( aGraph.size() == 0 ) return 0;

            // unflag all vertices and make sure that
            // indices are consistent
            index_t tIndex = 0 ;
            for( Vertex * tVertex : aGraph )
            {
                tVertex->set_index( tIndex++ );
            }

            // perform dfs
            proc_t tNumPartitions = dfs( aGraph );

            // count entries per partition
            Vector< index_t > tNumVertices ( tNumPartitions, 0 );
            for( Vertex * tVertex : aGraph )
            {
                ++tNumVertices( tVertex->owner() );
            }

            // Create a vector of pairs (size, partition) for sorting
            std::vector<std::pair<index_t, proc_t>> tSizePartitionPairs;
            for( proc_t p = 0; p < tNumPartitions; ++p )
            {
                tSizePartitionPairs.emplace_back(tNumVertices( p ), p);
            }

            // Sort partitions by size in descending order
            std::sort(tSizePartitionPairs.begin(), tSizePartitionPairs.end(), std::greater<>());

            // Rearrange the indices based on their size, starting with the largest
            Vector< proc_t > tGroup( tNumPartitions );
            for( proc_t q = 0; q < tNumPartitions; ++q )
            {
                // Get old partition index and assign new one
                tGroup( tSizePartitionPairs[q].second ) = q;
            }

            // Update group numbers in the graph
            for( Vertex * tVertex : aGraph )
            {
                tVertex->set_owner( tGroup( tVertex->owner() ) );
            }

            // Finally, rearrange the graph based on the owner index
            sort( aGraph, opVertexOwner );

            // Return the number of vertices in the largest group
            return tSizePartitionPairs.front().first;
        }
    }
}