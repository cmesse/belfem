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

#ifndef BELFEM_CL_ELEMENT_FACTORY_HPP
#define BELFEM_CL_ELEMENT_FACTORY_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"

#include "cl_Element.hpp"
#include "cl_Matrix.hpp"
#include "cl_ReferenceElement.hpp"

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

            ReferenceElement *
            create_reference_element(
                    const ElementType aType ) const;

            void
            create_orientation_table ( const ElementType aType, Matrix< uint > & aTable ) const;

//------------------------------------------------------------------------------
        private:

            void
            create_unity_nodes( const ElementType aType, Matrix< real > & aNodes ) const;


        };
    }
}
#endif //BELFEM_CL_ELEMENT_FACTORY_HPP
