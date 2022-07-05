//
// Created by Christian Messe on 04.05.20.
//

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
