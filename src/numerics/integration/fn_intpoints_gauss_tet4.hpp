//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET4_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET4_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------
        inline void
        gauss_tet4(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            real tAlpha = ( 5.0 - std::sqrt( 5.0 ) )/20.0;
            real tBeta = 1.0 - 3.0*tAlpha;

            aPoints.set_size( 4, 4 );

            aPoints( 0, 0 ) = tBeta;
            aPoints( 1, 0 ) = tAlpha;
            aPoints( 2, 0 ) = tAlpha;
            aPoints( 3, 0 ) = tAlpha;

            aPoints( 0, 1 ) = tAlpha;
            aPoints( 1, 1 ) = tBeta;
            aPoints( 2, 1 ) = tAlpha;
            aPoints( 3, 1 ) = tAlpha;

            aPoints( 0, 2 ) = tAlpha;
            aPoints( 1, 2 ) = tAlpha;
            aPoints( 2, 2 ) = tBeta;
            aPoints( 3, 2 ) = tAlpha;

            aPoints( 0, 3 ) = tAlpha;
            aPoints( 1, 3 ) = tAlpha;
            aPoints( 2, 3 ) = tAlpha;
            aPoints( 3, 3 ) = tBeta;

            aWeights.set_size( 4, 1.0/24.0 );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET4_HPP
