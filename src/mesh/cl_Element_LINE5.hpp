//
// Created by Christian Messe on 25.07.20.
//

#ifndef BELFEM_CL_ELEMENT_LINE5_HPP
#define BELFEM_CL_ELEMENT_LINE5_HPP

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
        LagrangeElement< 5, 2, 1, 0, 0 >::type() const
        {
            return ElementType::LINE5;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 5, 2, 1, 0, 0 >::get_nodes_of_edge(
                const uint aFacetIndex,
                Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 5, nullptr );

            aNodes( 0 ) = mNodes[ 0 ];
            aNodes( 1 ) = mNodes[ 1 ];
            aNodes( 2 ) = mNodes[ 2 ];
            aNodes( 3 ) = mNodes[ 3 ];
            aNodes( 4 ) = mNodes[ 4 ];
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_ELEMENT_LINE5_HPP
