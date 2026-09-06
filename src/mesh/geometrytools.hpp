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

#ifndef BELFEM_GEOMETRYTOOLS_HPP
#define BELFEM_GEOMETRYTOOLS_HPP

#include "typedefs.hpp"
#include "cl_Element.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace mesh
    {
//------------------------------------------------------------------------------

        void
        collect_node_coords(
                Element        * aElement,
                Matrix< real > & aNodeCoords,
                const uint aNumberOfDimensions = 0,
                const bool aCornerNodesOnly = false );

//------------------------------------------------------------------------------

        void
        collect_node_data(
                Element              * aElement,
                const Vector< real > & aFieldOnMesh ,
                Vector< real >       & aNodalField,
                const bool aCornerNodesOnly = false );

//------------------------------------------------------------------------------

        // Bronstein, Eq. 8.152c
        real
        compute_surface_increment(
                const Matrix< real > & adNdXi,
                const Matrix< real > & aNodeCoords,
                const uint             aNumSpatialDimensions=0 );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_GEOMETRYTOOLS_HPP
