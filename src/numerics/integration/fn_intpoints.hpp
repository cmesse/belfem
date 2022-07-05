//
// Created by Christian Messe on 2019-01-15.
//

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