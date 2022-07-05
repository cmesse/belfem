//
// Created by Christian Messe on 2019-08-04.
//

#ifndef BELFEM_CL_ELEMENT_LINE3_HPP
#define BELFEM_CL_ELEMENT_LINE3_HPP

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
        LagrangeElement< 3, 2, 1, 0, 0 >::type() const
        {
            return ElementType::LINE3;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 3, 2, 1, 0, 0 >::get_nodes_of_edge(
                const uint aFacetIndex,
                Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 3, nullptr );

            aNodes( 0 ) = mNodes[ 0 ];
            aNodes( 1 ) = mNodes[ 1 ];
            aNodes( 2 ) = mNodes[ 2 ];
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_ELEMENT_LINE3_HPP
