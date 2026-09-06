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

#ifndef BELFEM_CL_ELEMENT_LINE4_HPP
#define BELFEM_CL_ELEMENT_LINE4_HPP

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
        ElementTemplate< 4, 2, 1, 0, 0 >::type() const
        {
            return ElementType::LINE4;
        }

//------------------------------------------------------------------------------

        template <>
        uint
        ElementTemplate< 4, 2, 1, 0, 0 >::dimension() const
        {
            return 1 ;
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 4, 2, 1, 0, 0 >::get_nodes_of_edge(
                const uint aFacetIndex,
                Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 4, nullptr );

            aNodes( 0 ) = mNodes[ 0 ];
            aNodes( 1 ) = mNodes[ 1 ];
            aNodes( 2 ) = mNodes[ 2 ];
            aNodes( 3 ) = mNodes[ 3 ];
       }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_ELEMENT_LINE4_HPP
