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

#ifndef BELFEM_OP_GRAPH_VERTEX_INDEX_HPP
#define BELFEM_OP_GRAPH_VERTEX_INDEX_HPP

#include "cl_Graph_Vertex.hpp"

namespace belfem
{

//------------------------------------------------------------------------------

    // comparison object
    inline struct OpVertexIndex
    {
        inline bool
        operator()( const graph::Vertex * aA, const graph::Vertex * aB )
        {
            return aA->index() < aB->index();
        }
    } opVertexIndex;

//------------------------------------------------------------------------------
}

#endif //BELFEM_OP_GRAPH_VERTEX_INDEX_HPP
