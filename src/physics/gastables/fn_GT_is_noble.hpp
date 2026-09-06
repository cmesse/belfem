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

#ifndef BELFEM_FN_GT_IS_NOBLE
#define BELFEM_FN_GT_IS_NOBLE

#include "typedefs.hpp"

namespace belfem
{
    namespace gastables
    {
        inline bool
        is_noble( const string & aLabel )
        {
            return aLabel == "He" || aLabel == "Ne" ||
                   aLabel == "Ar" || aLabel == "Kr" ||
                   aLabel == "Xe" || aLabel == "Rn";
        }
    }
}
#endif // BELFEM_FN_GT_IS_NOBLE
