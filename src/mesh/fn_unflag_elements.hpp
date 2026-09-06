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

#ifndef BELFEM_FN_UNFLAG_ELEMENTS_HPP
#define BELFEM_FN_UNFLAG_ELEMENTS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
        inline void
        unflag_elements( Cell< Element * > & aElements )
        {
            for ( Element * tElement : aElements )
            {
                tElement->unflag();
            }
        }
    }
}
#endif //BELFEM_FN_UNFLAG_ELEMENTS_HPP
