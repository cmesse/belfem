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

#ifndef BELFEM_FN_INTPOINTS_HPP
#define BELFEM_FN_INTPOINTS_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "assert.hpp"
#include "meshtools.hpp"
#include "en_IntegrationScheme.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    void
    intpoints( const enum IntegrationScheme aIntegrationScheme,
               const enum GeometryType      aGeometryType,
               const               uint     aOrder,
                          Vector< real >  & aWeights,
                          Matrix< real >  & aPoints );

//------------------------------------------------------------------------------
    namespace integration
    {
//------------------------------------------------------------------------------
        void
        gauss_line(
                const int       aOrder,
                Vector <real> & aWeights,
                Matrix <real> & aPoints );

//------------------------------------------------------------------------------

        void
        gauss_quad(
                const int       aOrder,
                Vector <real> & aWeights,
                Matrix <real> & aPoints );

 //------------------------------------------------------------------------------

        void
        gauss_hex(
                const int       aOrder,
                Vector <real> & aWeights,
                Matrix <real> & aPoints );

//------------------------------------------------------------------------------
    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_INTPOINTS_HPP

#include "fn_intpoints_gauss_tri3.hpp"