//
// Created by Christian Messe on 16.12.20.
//

#ifndef BELFEM_FN_MESH_INTEGRATE_SCALAR_OVER_SIDESETS_HPP
#define BELFEM_FN_MESH_INTEGRATE_SCALAR_OVER_SIDESETS_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        real
        integrate_scalar_over_sidesets(
                Mesh * aMesh,
                const string         & aFieldLabel,
                const Vector< id_t > & aSidesetIDs,
                const proc_t           aMasterRank );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_MESH_INTEGRATE_SCALAR_OVER_SIDESETS_HPP
