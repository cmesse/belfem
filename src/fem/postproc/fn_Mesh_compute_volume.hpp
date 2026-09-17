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
