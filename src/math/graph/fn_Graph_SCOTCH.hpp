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

#ifndef BELFEM_FN_GRAPH_SCOTCH_HPP
#define BELFEM_FN_GRAPH_SCOTCH_HPP

#ifdef BELFEM_SCOTCH
#include <scotch.h>
#endif

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {

        void
        scotch_nd( Graph & aGraph );

        void
        scotch_partition( Graph & aGraph, const uint aNumPartitions );

    }
}
#endif //BELFEM_FN_GRAPH_SCOTCH_HPP