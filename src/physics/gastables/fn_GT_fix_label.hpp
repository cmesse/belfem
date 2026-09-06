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

#ifndef BELFEM_FN_GT_FIX_LABEL_HPP
#define BELFEM_FN_GT_FIX_LABEL_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
         *
         * fix illogical naming in CEA, eg AL --> Al, CL --> Cl
         */
        string
        fix_label( const string & aLabel );

    }
}

#endif //BELFEM_FN_GT_FIX_LABEL_HPP
