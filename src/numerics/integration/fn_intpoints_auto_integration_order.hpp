//
// Created by Christian Messe on 16.06.20.
//

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
