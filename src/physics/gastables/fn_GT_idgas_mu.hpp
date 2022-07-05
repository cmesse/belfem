//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_FN_GT_IDGAS_MU_HPP
#define BELFEM_FN_GT_IDGAS_MU_HPP

#include "typedefs.hpp"
#include "cl_GT_RefGas.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        /**
         * calculate the viscosity in case there is no data in the database
         */
        real
        idgas_mu( RefGas * aGas, const real aT );

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_FN_GT_IDGAS_MU_HPP
