//
// Created by Christian Messe on 03.11.19.
//

#ifndef BELFEM_FN_FEM_INITIALIZE_INTEGRATION_POINTS_HPP
#define BELFEM_FN_FEM_INITIALIZE_INTEGRATION_POINTS_HPP

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

        // by default we want twice the order of the element
        void
        initialize_integration_points(
                const ElementType    & aElementType,
                      Vector< real > & aWeights,
                      Matrix< real > & aPoints,
                const uint             aIntegrationOrder=0,
                const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

        // by default we want twice the order of the element
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
#endif //BELFEM_FN_FEM_INITIALIZE_INTEGRATION_POINTS_HPP
