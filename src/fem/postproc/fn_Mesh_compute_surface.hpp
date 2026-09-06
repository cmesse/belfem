//
// Created by Christian Messe on 27.06.20.
//

#ifndef BELFEM_FN_MESH_COMPUTE_SURFACE_HPP
#define BELFEM_FN_MESH_COMPUTE_SURFACE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        real
        compute_surface( Mesh * aMesh, const Vector< id_t > & aSideSetIDs );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_MESH_COMPUTE_SURFACE_HPP
