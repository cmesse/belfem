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

#ifndef BELFEM_FN_GRAPH_BFS_HPP
#define BELFEM_FN_GRAPH_BFS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        index_t
        bfs( Graph & aGraph, Vertex * aStart );

        // can also handle non-connected graphs
        index_t
        bfs( Graph & aGraph );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GRAPH_BFS_HPP
