//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_ELEMENT_FACTORY_HPP
#define BELFEM_CL_ELEMENT_FACTORY_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"

#include "cl_Element.hpp"
#include "cl_ElementTemplate.hpp"

namespace belfem
{
    namespace mesh
    {
        class ElementFactory
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ElementFactory() = default;

            ~ElementFactory() = default;

//------------------------------------------------------------------------------

            Element *
            create_element(
                    const ElementType aType,
                    const id_t aID ) const;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_ELEMENT_FACTORY_HPP
