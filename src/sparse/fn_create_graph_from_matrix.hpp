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

#ifndef BELFEM_FN_CREATE_GRAPH_FROM_MATRIX_HPP
#define BELFEM_FN_CREATE_GRAPH_FROM_MATRIX_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_SpMatrix.hpp"

namespace belfem
{
    namespace sparse
    {

        void
        create_graph_from_matrix(
        const SpMatrix & aMatrix,
            Graph & aGraph );

    }
}

#endif //BELFEM_FN_CREATE_GRAPH_FROM_MATRIX_HPP