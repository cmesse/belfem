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

#ifndef BELFEM_VTKTOOLS_HPP
#define BELFEM_VTKTOOLS_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace vtk
    {
//------------------------------------------------------------------------------

        /**
          * This function inverts little endian to big endian and vice versa.
          * Needed for VTK output files.
          */
        template <typename T> T swap_byte_endian(T aValue)
        {
            T aOutValue;
            auto *tPointer = (char*) & aValue;
            auto *tOutPointer = (char*) & aOutValue;
            int size = sizeof(T);
            for( int i=0; i<size; i++ )
            {
                tOutPointer[size - 1 - i] = tPointer[i];
            }
            return aOutValue;
        }

//------------------------------------------------------------------------------

        /**
         * returns the VTK type of an element
         */
        int
        vtk_type( const ElementType & aElementType );

//------------------------------------------------------------------------------

        /**
         * get node IDs of this element in VTK order
         */
         void
         get_node_ids( mesh::Element * aElement, Vector< id_t > & aNodeIDs );

//-----------------------------------------------------------------------------
    }
}
#endif //BELFEM_VTKTOOLS_HPP
