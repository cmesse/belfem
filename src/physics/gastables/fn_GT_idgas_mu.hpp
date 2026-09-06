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
        idgas_mu( RefGas * aGas, const real T );

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_FN_GT_IDGAS_MU_HPP
