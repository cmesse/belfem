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

#ifndef BELFEM_EN_OPT_ALGORITHM_HPP
#define BELFEM_EN_OPT_ALGORITHM_HPP

namespace belfem
{
    namespace opt
    {
        /**
         * optimization algorithm selector. The names are library-agnostic; the
         * mapping to the concrete nlopt algorithm lives in cl_Optimizer.cpp so
         * that no nlopt header leaks into the public interface.
         *
         * Naming follows nlopt's own convention where L = local, N = no
         * derivatives, D = uses derivatives.
         */
        enum class Algorithm
        {
            // --- local, derivative-free ----------------------------------
            BOBYQA     = 0, // LN_BOBYQA   - quadratic model, bound constrained
            COBYLA     = 1, // LN_COBYLA   - linear model, supports constraints
            NELDERMEAD = 2, // LN_NELDERMEAD - simplex
            SBPLX      = 3, // LN_SBPLX    - Rowan's subplex
            PRAXIS     = 4, // LN_PRAXIS   - principal-axis

            // --- local, gradient-based -----------------------------------
            MMA        = 5, // LD_MMA      - method of moving asymptotes
            SLSQP      = 6, // LD_SLSQP    - sequential quadratic programming
            LBFGS      = 7  // LD_LBFGS    - low-storage BFGS
        };
    }
}
#endif //BELFEM_EN_OPT_ALGORITHM_HPP
