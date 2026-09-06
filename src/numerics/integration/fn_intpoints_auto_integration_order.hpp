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

#ifndef BELFEM_FN_INTPOINTS_AUTO_ORDER_HPP
#define BELFEM_FN_INTPOINTS_AUTO_ORDER_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * tries to automatically detect the appropriate integration order
     * @param aType
     * @return
     */
    uint
    auto_integration_order( const ElementType aType );

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_INTPOINTS_AUTO_ORDER_HPP
