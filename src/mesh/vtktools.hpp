//
// Created by Christian Messe on 22.06.20.
//

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
