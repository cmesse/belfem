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

#ifndef FN_GRAPH_FIND_CONNECTED_PARTITIONS_HPP
#define FN_GRAPH_FIND_CONNECTED_PARTITIONS_HPP

#include "fn_Graph_dfs.hpp"

namespace belfem
{
    namespace graph
    {
        index_t
        find_connected_partitions( Graph & aGraph );

    }
}

#endif //FN_GRAPH_FIND_CONNECTED_PARTITIONS_HPP
