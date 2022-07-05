//
// Created by Christian Messe on 14.06.20.
//

#ifndef BELFEM_FN_MESH_COMPUTE_EDGE_LENGTHS_HPP
#define BELFEM_FN_MESH_COMPUTE_EDGE_LENGTHS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        void
        compute_edge_lengths(
                const uint                aDimension,
                Cell< mesh::Element * > & aElements,
                Vector< real >          & aEdgeLength );

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_FN_MESH_COMPUTE_EDGE_LENGTHS_HPP
