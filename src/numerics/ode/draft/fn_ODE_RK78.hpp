//
// Created by Christian Messe on 27.12.19.
//

#ifndef BELFEM_FN_ODE_RK78_HPP
#define BELFEM_FN_ODE_RK78_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_ODE.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    int
    RK78(    ODE             & aODE,
             real                   & aT,
             Vector< real >         & aY,
             real                   & aStep,
             Cell< Vector< real > > & aWork,
             const real               aEpsilon=1e-7,
             const uint               aMaxIterations = 1000,
             const real               aXmax = BELFEM_REAL_MAX );

//------------------------------------------------------------------------------
}


#endif //BELFEM_FN_ODE_RK78_HPP
