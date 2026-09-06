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

#ifndef BELFEM_FN_GM_HELMHOLZ_DERIVTEST_HPP
#define BELFEM_FN_GM_HELMHOLZ_DERIVTEST_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GM_Helmholtz.hpp"

namespace belfem
{
    namespace gasmodels
    {
        void
        deriv_test( Helmholtz & aGas, Vector< real > & aR2 );

    }
}
#endif //BELFEM_FN_GM_HELMHOLZ_DERIVTEST_HPP
