//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET5_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET5_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        inline void
        gauss_tet5(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            real tBeta = 1.0 / 6.0;

            aPoints.set_size( 3, 5 );

            aPoints( 0, 0 ) = 0.25;
            aPoints( 1, 0 ) = 0.25;
            aPoints( 2, 0 ) = 0.25;

            aPoints( 0, 1 ) = 0.5;
            aPoints( 1, 1 ) = tBeta;
            aPoints( 2, 1 ) = tBeta;

            aPoints( 0, 2 ) = tBeta;
            aPoints( 1, 2 ) = 0.5;
            aPoints( 2, 2 ) = tBeta;

            aPoints( 0, 3 ) = tBeta;
            aPoints( 1, 3 ) = tBeta;
            aPoints( 2, 3 ) = 0.5;

            aPoints( 0, 4 ) = tBeta;
            aPoints( 1, 4 ) = tBeta;
            aPoints( 2, 4 ) = tBeta;;

            aWeights.set_size( 5 );


            aWeights( 0 ) = -4.0/30;
            aWeights( 1 ) = 0.075;
            aWeights( 2 ) = 0.075;
            aWeights( 3 ) = 0.075;
            aWeights( 4 ) = 0.075;
        }
// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET5_HPP
