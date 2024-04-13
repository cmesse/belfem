//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_ELEMENT_LINE2_HPP
#define BELFEM_CL_ELEMENT_LINE2_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_ElementTemplate.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        template <>
        ElementType
        ElementTemplate< 2, 2, 1, 0, 0 >::type() const
        {
            return ElementType::LINE2;
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 2, 2, 1, 0, 0 >::get_nodes_of_edge(
                const uint aFacetIndex,
                Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            aNodes( 0 ) = mNodes[ 0 ];
            aNodes( 1 ) = mNodes[ 1 ];
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_ELEMENT_LINE2_HPP
