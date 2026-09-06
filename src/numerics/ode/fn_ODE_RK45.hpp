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

#ifndef BELFEM_FN_ODE_RK45_HPP
#define BELFEM_FN_ODE_RK45_HPP

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

        // doi: 10.1007/BF02241732
        Status
        RK45( ODE         & aODE,
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
        RK45_init( ODE & aODE, Cell< Vector< real > > & aWork );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_ODE_RK45_HPP
