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

#ifndef BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_NODE_HPP
#define BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_NODE_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        Vertex *
        find_pseudo_peripheral_node(
                Graph & aGraph,
                      Vertex *     aStart = nullptr );

//------------------------------------------------------------------------------
    }
}
#endif // BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_NODE_HPP
