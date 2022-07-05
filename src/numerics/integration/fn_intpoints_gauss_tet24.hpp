//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET24_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET24_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        inline void
        gauss_tet24(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            // source Keast: 10.1016/0045-7825(86)90059-9

            real tAlpha1 = 0.214602871259151684;
            real tBeta1  = 1.0 - 3.0*tAlpha1;

            real tAlpha2 = 0.0406739585346113397;
            real tBeta2  = 1.0 - 3.0*tAlpha2;

            real tAlpha3 = 0.322337890142275646;
            real tBeta3  = 1.0 - 3.0*tAlpha3;

            real tAlpha4 = 0.0636610018750175299;
            real tBeta4  = 1.0/3.0 - tAlpha4;
            real tDelta4 = 1.0 - 2.0*tAlpha4 - tBeta4;

            aPoints.set_size( 3, 24 );

            aPoints( 0, 0 ) = tAlpha1;
            aPoints( 1, 0 ) = tAlpha1;
            aPoints( 2, 0 ) = tAlpha1;

            aPoints( 0, 1 ) = tAlpha1;
            aPoints( 1, 1 ) = tAlpha1;
            aPoints( 2, 1 ) = tBeta1;

            aPoints( 0, 2 ) = tAlpha1;
            aPoints( 1, 2 ) = tBeta1;
            aPoints( 2, 2 ) = tAlpha1;

            aPoints( 0, 3 ) = tBeta1;
            aPoints( 1, 3 ) = tAlpha1;
            aPoints( 2, 3 ) = tAlpha1;

            aPoints( 0, 4 ) = tAlpha2;
            aPoints( 1, 4 ) = tAlpha2;
            aPoints( 2, 4 ) = tAlpha2;

            aPoints( 0, 5 ) = tAlpha2;
            aPoints( 1, 5 ) = tAlpha2;
            aPoints( 2, 5 ) = tBeta2;

            aPoints( 0, 6 ) = tAlpha2;
            aPoints( 1, 6 ) = tBeta2;
            aPoints( 2, 6 ) = tAlpha2;

            aPoints( 0, 7 ) = tBeta2;
            aPoints( 1, 7 ) = tAlpha2;
            aPoints( 2, 7 ) = tAlpha2;

            aPoints( 0, 8 ) = tAlpha3;
            aPoints( 1, 8 ) = tAlpha3;
            aPoints( 2, 8 ) = tAlpha3;

            aPoints( 0, 9 ) = tAlpha3;
            aPoints( 1, 9 ) = tAlpha3;
            aPoints( 2, 9 ) = tBeta3;

            aPoints( 0, 10 ) = tAlpha3;
            aPoints( 1, 10 ) = tBeta3;
            aPoints( 2, 10 ) = tAlpha3;

            aPoints( 0, 11 ) = tBeta3;
            aPoints( 1, 11 ) = tAlpha3;
            aPoints( 2, 11 ) = tAlpha3;

            aPoints( 0, 12 ) = tAlpha4;
            aPoints( 1, 12 ) = tAlpha4;
            aPoints( 2, 12 ) = tBeta4;

            aPoints( 0, 13 ) = tAlpha4;
            aPoints( 1, 13 ) = tAlpha4;
            aPoints( 2, 13 ) = tDelta4;

            aPoints( 0, 14 ) = tAlpha4;
            aPoints( 1, 14 ) = tBeta4;
            aPoints( 2, 14 ) = tAlpha4;

            aPoints( 0, 15 ) = tAlpha4;
            aPoints( 1, 15 ) = tBeta4;
            aPoints( 2, 15 ) = tDelta4;

            aPoints( 0, 16 ) = tAlpha4;
            aPoints( 1, 16 ) = tDelta4;
            aPoints( 2, 16 ) = tAlpha4;

            aPoints( 0, 17 ) = tAlpha4;
            aPoints( 1, 17 ) = tDelta4;
            aPoints( 2, 17 ) = tBeta4;

            aPoints( 0, 18 ) = tBeta4;
            aPoints( 1, 18 ) = tAlpha4;
            aPoints( 2, 18 ) = tAlpha4;

            aPoints( 0, 19 ) = tBeta4;
            aPoints( 1, 19 ) = tAlpha4;
            aPoints( 2, 19 ) = tDelta4;

            aPoints( 0, 20 ) = tBeta4;
            aPoints( 1, 20 ) = tDelta4;
            aPoints( 2, 20 ) = tAlpha4;

            aPoints( 0, 21 ) = tDelta4;
            aPoints( 1, 21 ) = tAlpha4;
            aPoints( 2, 21 ) = tAlpha4;

            aPoints( 0, 22 ) = tDelta4;
            aPoints( 1, 22 ) = tAlpha4;
            aPoints( 2, 22 ) = tBeta4;

            aPoints( 0, 23 ) = tDelta4;
            aPoints( 1, 23 ) = tBeta4;
            aPoints( 2, 23 ) = tAlpha4;

            aWeights.set_size( 24 );

            aWeights(  0 ) = 0.00665379170969464506;
            aWeights(  1 ) = 0.00665379170969464506;
            aWeights(  2 ) = 0.00665379170969464506;
            aWeights(  3 ) = 0.00665379170969464506;

            aWeights(  4 ) = 0.00167953517588677620;
            aWeights(  5 ) = 0.00167953517588677620;
            aWeights(  6 ) = 0.00167953517588677620;
            aWeights(  7 ) = 0.00167953517588677620;

            aWeights(  8 ) = 0.00922619692394239843;
            aWeights(  9 ) = 0.00922619692394239843;
            aWeights( 10 ) = 0.00922619692394239843;
            aWeights( 11 ) = 0.00922619692394239843;

            aWeights( 12 ) = 0.00803571428571428248;
            aWeights( 13) = 0.00803571428571428248;
            aWeights( 14 ) = 0.00803571428571428248;
            aWeights( 15 ) = 0.00803571428571428248;

            aWeights( 16 ) = 0.00803571428571428248;
            aWeights( 17 ) = 0.00803571428571428248;
            aWeights( 18 ) = 0.00803571428571428248;
            aWeights( 19 ) = 0.00803571428571428248;

            aWeights( 20 ) = 0.00803571428571428248;
            aWeights( 21 ) = 0.00803571428571428248;
            aWeights( 22 ) = 0.00803571428571428248;
            aWeights( 23 ) = 0.00803571428571428248;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET24_HPP
