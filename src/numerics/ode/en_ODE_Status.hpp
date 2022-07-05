//
// Created by Christian Messe on 27.12.19.
//

#ifndef BELFEM_EN_ODE_STATUS_HPP
#define BELFEM_EN_ODE_STATUS_HPP

namespace belfem
{
    namespace ode
    {
        enum class Status
        {
            OK      = 0 , // all OK
            TRAPPED = 1,  // all OK, but step was trapped
            MAXIT   = 2,  // too many iterations
        };
    }
}
#endif //BELFEM_EN_ODE_STATUS_HPP
