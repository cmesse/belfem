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

#ifndef BELFEM_FN_ODE_DOP853_HPP
#define BELFEM_FN_ODE_DOP853_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_ODE.hpp"
#include "en_ODE_Status.hpp"

namespace belfem
{
    namespace ode
    {
//------------------------------------------------------------------------------

        /**
         * Dormand-Prince 8(5,3) explicit Runge-Kutta method (DOP853).
         *
         * An 8th order method with embedded 5th and 3rd order error
         * estimators. The step size controller uses the 5th order
         * estimator. Recommended for smooth, non-stiff problems
         * requiring high accuracy (tolerances tighter than ~1e-8).
         * For moderate accuracy, RK45 is more efficient.
         *
         * Reference: Hairer, Norsett, Wanner: Solving Ordinary
         * Differential Equations I, 2nd ed., Springer (1993),
         * Section II.6.
         */
        Status
        DOP853( ODE         & aODE,
                real           & aT,
                Vector< real > & aY,
                real           & aStep,
                Cell< Vector< real > > & aWork,
                const real aEpsilon = 1e-7,
                const uint aMaxIterations = 1000,
                const real aTmax = BELFEM_REAL_MAX,
                const bool aAutoTimestep = true );

//------------------------------------------------------------------------------

        void
        DOP853_init( ODE & aODE, Cell< Vector< real > > & aWork );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_ODE_DOP853_HPP
