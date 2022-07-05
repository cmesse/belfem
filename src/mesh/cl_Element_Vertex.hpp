//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_ELEMENT_VERTEX_HPP
#define BELFEM_CL_ELEMENT_VERTEX_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_LagrangeElement.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        template <>
        ElementType
        LagrangeElement< 1, 1, 0, 0, 0 >::type() const
        {
            return ElementType::VERTEX;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_ELEMENT_VERTEX_HPP
