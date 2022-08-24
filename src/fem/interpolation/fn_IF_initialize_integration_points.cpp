//
// Created by Christian Messe on 03.11.19.
//

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
            // by default we want twice the order of the element
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