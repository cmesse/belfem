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

#ifndef BELFEM_FN_IF_INITIALIZE_INTEGRATION_POINTS_ON_FACET_HPP
#define BELFEM_FN_IF_INITIALIZE_INTEGRATION_POINTS_ON_FACET_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"
#include "en_IntegrationScheme.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        void
        initialize_integration_points_on_facet(
                const ElementType    aElementType,
                const uint           aSideIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint             aIntegrationOrder=0,
                const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

        void
        initialize_integration_points_on_facet(
            const ElementType    aElementType,
            const uint           aSideIndex,
            const uint           aOrientation,
            Vector< real > & aWeights,
            Matrix< real > & aPoints,
            const uint             aIntegrationOrder=0,
            const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

        namespace facetintpoints
        {
//------------------------------------------------------------------------------

            void
            intpoints_tri(
                    const uint       aMasterIndex,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            void
            intpoints_quad(
                    const uint       aMasterIndex,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            void
            intpoints_tet(
                    const uint       aMasterIndex,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            void
            intpoints_penta(
                    const uint       aMasterIndex,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            void
            intpoints_hex(
                    const uint       aMasterIndex,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            void
            intpoints_pyra(
                    const uint       aMasterIndex,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            /**
             * special function for slave side
             */
            void
            intpoints_quad(
                    const uint       aSlaveIndex,
                    const uint       aOrientation,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            /**
             * special function for slave side
             */
            void
            intpoints_pyra(
                    const uint       aSlaveIndex,
                    const uint       aOrientation,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            /**
             * special function for slave side
             */
            void
            intpoints_penta(
                    const uint       aSlaveIndex,
                    const uint       aOrientation,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            /**
             * special function for slave side
             */
            void
            intpoints_tet(
                    const uint       aSlaveIndex,
                    const uint       aOrientation,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            /**
             * special function for slave side
             */
            void
            intpoints_hex(
                    const uint       aSlaveIndex,
                    const uint       aOrientation,
                    Vector< real > & aWeights,
                    Matrix< real > & aPoints,
                    const uint       aIntegrationOrder,
                    const IntegrationScheme  aIntegrationScheme=IntegrationScheme::GAUSS );


//------------------------------------------------------------------------------
        }
    }
}
#endif //BELFEM_FN_IF_INITIALIZE_INTEGRATION_POINTS_ON_FACET_HPP
