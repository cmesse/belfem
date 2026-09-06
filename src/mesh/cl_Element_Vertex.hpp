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

#ifndef BELFEM_CL_ELEMENT_VERTEX_HPP
#define BELFEM_CL_ELEMENT_VERTEX_HPP

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
        ElementTemplate< 1, 1, 0, 0, 0 >::type() const
        {
            return ElementType::VERTEX;
        }

//------------------------------------------------------------------------------

        template <>
        uint
        ElementTemplate< 1, 1, 0, 0, 0 >::dimension() const
        {
            return 0 ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_ELEMENT_VERTEX_HPP
