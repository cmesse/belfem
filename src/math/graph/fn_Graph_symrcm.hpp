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

#ifndef BELFEM_FN_GRAPH_SYMCRM_HPP
#define BELFEM_FN_GRAPH_SYMCRM_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
          * Performs Reverse Cuthill-McKee (RCM) ordering on the graph.
          * This algorithm reorders vertices to reduce the bandwidth of the adjacency matrix.
          * The graph is reordered in-place.
          *
          * @param aGraph     Cell containing all vertices in the graph (will be reordered)
          * @param aStart     Starting vertex for the algorithm (if nullptr, uses pseudo-peripheral vertex)
          */

        void
        symrcm( Graph & aGraph,
                      Vertex *     aStart = nullptr );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GRAPH_SYMCRM_HPP
