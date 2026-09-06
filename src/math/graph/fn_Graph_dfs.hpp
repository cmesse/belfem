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

#ifndef FN_GRAPH_DFS_HPP
#define FN_GRAPH_DFS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {


        proc_t
        dfs( Graph & aGraph );

        void
        dfs_from_start(  Graph & aGraph, Vertex * aStart );

    }
}
#endif //FN_GRAPH_DFS_HPP
