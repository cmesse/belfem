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

#include "fn_IF_initialize_integration_points.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_intpoints.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        void
        initialize_integration_points(
                const ElementType    & aElementType,
                Vector< real >       & aWeights,
                Matrix< real >       & aPoints,
                const uint             aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )

        {
            // aIntegrationOrder == 0 selects auto_integration_order() ( 4 / 7 / 10 for linear / quadratic / cubic elements )
            uint tIntegrationOrder = aIntegrationOrder == 0 ?
                    auto_integration_order( aElementType ) :
                    aIntegrationOrder ;

            // find out geometry type of this element
            GeometryType tGeometryType = mesh::geometry_type( aElementType );

            initialize_integration_points(
                tGeometryType,
                aWeights,
                aPoints,
                tIntegrationOrder,
                aIntegrationScheme ) ;
        }

//------------------------------------------------------------------------------

        void
        initialize_integration_points(
                const GeometryType   & aGeometryType,
                Vector< real >       & aWeights,
                Matrix< real >       & aPoints,
                const uint             aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            intpoints(
                    aIntegrationScheme,
                    aGeometryType,
                    aIntegrationOrder,
                    aWeights,
                    aPoints );
        }

//------------------------------------------------------------------------------

    }
}