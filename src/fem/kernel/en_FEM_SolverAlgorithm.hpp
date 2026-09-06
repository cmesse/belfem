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

#ifndef BELFEM_EN_FEM_SOLVERALGORITHM_HPP
#define BELFEM_EN_FEM_SOLVERALGORITHM_HPP

namespace belfem
{
    namespace fem
    {
        enum class SolverAlgorithm
        {
            Direct        = 0,
            NewtonRaphson = 1,
            Picard        = 2,
            UNDEFINED     = 3
        };
    }
}

#endif //BELFEM_EN_FEM_SOLVERALGORITHM_HPP
