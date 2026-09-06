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

#ifndef BELFEM_EN_ODE_TYPE_HPP
#define BELFEM_EN_ODE_TYPE_HPP

namespace belfem
{
    namespace ode
    {
        enum class Type
        {
            RK45 = 0,
            DOP853 = 1
        };
    }
}

#endif //BELFEM_EN_ODE_TYPE_HPP
