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
