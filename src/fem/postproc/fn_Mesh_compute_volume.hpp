//
// Created by Christian Messe on 27.06.20.
//

#ifndef BELFEM_FN_MESH_COMPUTE_VOLUME_HPP
#define BELFEM_FN_MESH_COMPUTE_VOLUME_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        real
        compute_volume( Mesh * aMesh, const id_t aBlockIDs );

//------------------------------------------------------------------------------

        real
        compute_volume( Mesh * aMesh, const Vector< id_t > & aBlockIDs );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_MESH_COMPUTE_VOLUME_HPP
