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

#ifndef BELFEM_FN_IF_INITIALIZE_INTEGRATION_POINTS_HPP
#define BELFEM_FN_IF_INITIALIZE_INTEGRATION_POINTS_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "en_IntegrationScheme.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        // aIntegrationOrder == 0 selects auto_integration_order() ( 4 / 7 / 10 for linear / quadratic / cubic elements )
        void
        initialize_integration_points(
                const ElementType    & aElementType,
                      Vector< real > & aWeights,
                      Matrix< real > & aPoints,
                const uint             aIntegrationOrder=0,
                const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

        // aIntegrationOrder == 0 selects auto_integration_order() ( 4 / 7 / 10 for linear / quadratic / cubic elements )
        void
        initialize_integration_points(
                const GeometryType  & aGeometryType,
                Vector< real >      & aWeights,
                Matrix< real >      & aPoints,
                const uint            aIntegrationOrder,
                const IntegrationScheme aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_IF_INITIALIZE_INTEGRATION_POINTS_HPP
