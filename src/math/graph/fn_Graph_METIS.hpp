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

#ifndef BELFEM_FN_GRAPH_METIS_HPP
#define BELFEM_FN_GRAPH_METIS_HPP


#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
          * Performs nested dissection ordering on the graph using METIS_NodeND.
          * This algorithm reorders vertices to minimize fill-in during
          * sparse matrix factorization, producing a shallow elimination tree
          * that is optimal for direct solvers like STRUMPACK.
          *
          * The graph is reordered in-place. After calling this function:
          * - vertex->index() contains the new (permuted) index
          * - The graph Cell is sorted by the new indices
          * - vertex->id() is not touched; record the original index yourself
          *   before the call if you need the permutation
          *
          * @param aGraph  Cell containing all vertices in the graph (will be reordered)
          * 
          * @note Requires BELFEM_METIS to be defined ( set by USE_METIS )
          * @note Without METIS this raises BELFEM_ERROR. There is no symrcm fallback.
          */
        void
        metis_nd( Graph & aGraph );

//------------------------------------------------------------------------------

        /**
          * Performs nested dissection ordering with partitioning using METIS_NodeNDP.
          * This variant creates a specified number of top-level partitions before
          * applying nested dissection, which can produce better orderings for
          * parallel factorization with many processors.
          *
          * STRUMPACK documentation suggests this often works better than standard
          * NodeND for certain problem types (use --sp_enable_METIS_NodeNDP).
          *
          * The graph is reordered in-place. After calling this function:
          * - vertex->index() contains the new (permuted) index
          * - The graph Cell is sorted by the new indices
          * - vertex->id() is not touched; record the original index yourself
          *   before the call if you need the permutation
          *
          * @param aGraph          Cell containing all vertices (will be reordered)
          * @param aNumPartitions  Number of top-level partitions
          * 
          * @note Requires BELFEM_METIS to be defined ( set by USE_METIS ); without it,
          *       BELFEM_ERROR
          */
        void
        metis_ndp( Graph & aGraph, const uint aNumPartitions );

//------------------------------------------------------------------------------

        void
        metis_partition(
            Graph & aGraph,
            const uint aNumPartitions,
            const bool aForceContinuousPartitions = true,
            Vector< proc_t > * aPartitions = nullptr );

//------------------------------------------------------------------------------

        string
        metis_status( const int aStatus );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_GRAPH_METIS_HPP
