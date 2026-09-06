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

#ifndef BELFEM_FN_COMPUTE_PERMUTATION_HPP
#define BELFEM_FN_COMPUTE_PERMUTATION_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_SpMatrix.hpp"

namespace belfem
{
    namespace sparse
    {
        void
        compute_permutation(
            const SpMatrix & aMatrix ,
            Graph & aGraph ,
            Cell< int_t > & aForwardPermutation ,
            Cell< int_t > & aBackwardPermutation ,
            Cell< int_t > & aIndexPermutation );
    }
}

#endif //BELFEM_FN_COMPUTE_PERMUTATION_HPP