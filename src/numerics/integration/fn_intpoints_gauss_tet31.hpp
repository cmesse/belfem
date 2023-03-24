//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET31_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET31_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------
        inline void
        gauss_tet31(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {


            // source Keast: 10.1016/0045-7825(86)90059-9

            aPoints.set_size( 4, 31 );

            real tAlpha1 = 0.0782131923303186549;
            real tBeta1  = 1.0 - 3.0*tAlpha1;

            real tAlpha2 = 0.121843216663904411;
            real tBeta2  = 1.0 - 3.0*tAlpha2;

            real tAlpha3 = 0.332539164446420554;
            real tBeta3  = 1.0 - 3.0*tAlpha3;

            real tAlpha4 = 0.5;
            real tBeta4  = 0.0;

            real tAlpha5 = 0.1;
            real tBeta5  = 0.2;
            real tDelta5 = 0.6;

            aPoints( 0, 0 ) = 0.25;
            aPoints( 1, 0 ) = 0.25;
            aPoints( 2, 0 ) = 0.25;

            aPoints( 0, 1 ) = tAlpha1;
            aPoints( 1, 1 ) = tAlpha1;
            aPoints( 2, 1 ) = tAlpha1;

            aPoints( 0, 2 ) = tAlpha1;
            aPoints( 1, 2 ) = tAlpha1;
            aPoints( 2, 2 ) = tBeta1;

            aPoints( 0, 3 ) = tAlpha1;
            aPoints( 1, 3 ) = tBeta1;
            aPoints( 2, 3 ) = tAlpha1;

            aPoints( 0, 4 ) = tBeta1;
            aPoints( 1, 4 ) = tAlpha1;
            aPoints( 2, 4 ) = tAlpha1;

            aPoints( 0, 5 ) = tAlpha2;
            aPoints( 1, 5 ) = tAlpha2;
            aPoints( 2, 5 ) = tAlpha2;

            aPoints( 0, 6 ) = tAlpha2;
            aPoints( 1, 6 ) = tAlpha2;
            aPoints( 2, 6 ) = tBeta2;

            aPoints( 0, 7 ) = tAlpha2;
            aPoints( 1, 7 ) = tBeta2;
            aPoints( 2, 7 ) = tAlpha2;

            aPoints( 0, 8 ) = tBeta2;
            aPoints( 1, 8 ) = tAlpha2;
            aPoints( 2, 8 ) = tAlpha2;

            aPoints( 0, 9 ) = tAlpha3;
            aPoints( 1, 9 ) = tAlpha3;
            aPoints( 2, 9 ) = tAlpha3;

            aPoints( 0, 10 ) = tAlpha3;
            aPoints( 1, 10 ) = tAlpha3;
            aPoints( 2, 10 ) = tBeta3;

            aPoints( 0, 11 ) = tAlpha3;
            aPoints( 1, 11 ) = tBeta3;
            aPoints( 2, 11 ) = tAlpha3;

            aPoints( 0, 12 ) = tBeta3;
            aPoints( 1, 12 ) = tAlpha3;
            aPoints( 2, 12 ) = tAlpha3;

            aPoints( 0, 13 ) = tAlpha4;
            aPoints( 1, 13 ) = tAlpha4;
            aPoints( 2, 13 ) = tBeta4;

            aPoints( 0, 14 ) = tAlpha4;
            aPoints( 1, 14 ) = tBeta4;
            aPoints( 2, 14 ) = tAlpha4;

            aPoints( 0, 15 ) = tAlpha4;
            aPoints( 1, 15 ) = tBeta4;
            aPoints( 2, 15 ) = tBeta4;

            aPoints( 0, 16 ) = tBeta4;
            aPoints( 1, 16 ) = tAlpha4;
            aPoints( 2, 16 ) = tAlpha4;

            aPoints( 0, 17 ) = tBeta4;
            aPoints( 1, 17 ) = tAlpha4;
            aPoints( 2, 17 ) = tBeta4;

            aPoints( 0, 18 ) = tBeta4;
            aPoints( 1, 18 ) = tBeta4;
            aPoints( 2, 18 ) = tAlpha4;

            aPoints( 0, 19 ) = tAlpha5;
            aPoints( 1, 19 ) = tAlpha5;
            aPoints( 2, 19 ) = tBeta5;

            aPoints( 0, 20 ) = tAlpha5;
            aPoints( 1, 20 ) = tAlpha5;
            aPoints( 2, 20 ) = tDelta5;

            aPoints( 0, 21 ) = tAlpha5;
            aPoints( 1, 21 ) = tBeta5;
            aPoints( 2, 21 ) = tAlpha5;

            aPoints( 0, 22 ) = tAlpha5;
            aPoints( 1, 22 ) = tBeta5;
            aPoints( 2, 22 ) = tDelta5;

            aPoints( 0, 23 ) = tAlpha5;
            aPoints( 1, 23 ) = tDelta5;
            aPoints( 2, 23 ) = tAlpha5;

            aPoints( 0, 24 ) = tAlpha5;
            aPoints( 1, 24 ) = tDelta5;
            aPoints( 2, 24 ) = tBeta5;

            aPoints( 0, 25 ) = tBeta5;
            aPoints( 1, 25 ) = tAlpha5;
            aPoints( 2, 25 ) = tAlpha5;

            aPoints( 0, 26 ) = tBeta5;
            aPoints( 1, 26 ) = tAlpha5;
            aPoints( 2, 26 ) = tDelta5;

            aPoints( 0, 27 ) = tBeta5;
            aPoints( 1, 27 ) = tDelta5;
            aPoints( 2, 27 ) = tAlpha5;

            aPoints( 0, 28 ) = tDelta5;
            aPoints( 1, 28 ) = tAlpha5;
            aPoints( 2, 28 ) = tAlpha5;

            aPoints( 0, 29 ) = tDelta5;
            aPoints( 1, 29 ) = tAlpha5;
            aPoints( 2, 29 ) = tBeta5;

            aPoints( 0, 30 ) = tDelta5;
            aPoints( 1, 30 ) = tBeta5;
            aPoints( 2, 30 ) = tAlpha5;

            for( uint k=0; k<31; ++k )
            {
                aPoints( 3, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k )
                                   - aPoints( 2, k );
            }

            aWeights.set_size( 31 );

            aWeights(  0 ) =  0.0182642234661087939;

            aWeights(  1 ) =  0.0105999415244141609;
            aWeights(  2 ) =  0.0105999415244141609;
            aWeights(  3 ) =  0.0105999415244141609;
            aWeights(  4 ) =  0.0105999415244141609;

            aWeights(  5 ) = -0.0625177401143299494;
            aWeights(  6 ) = -0.0625177401143299494;
            aWeights(  7 ) = -0.0625177401143299494;
            aWeights(  8 ) = -0.0625177401143299494;

            aWeights(  9 ) =  0.00489142526307353653;
            aWeights( 10 ) =  0.00489142526307353653;
            aWeights( 11 ) =  0.00489142526307353653;
            aWeights( 12 ) =  0.00489142526307353653;

            aWeights( 13 ) =  0.000970017636684296702;
            aWeights( 14 ) =  0.000970017636684296702;
            aWeights( 15 ) =  0.000970017636684296702;
            aWeights( 16 ) =  0.000970017636684296702;
            aWeights( 17 ) =  0.000970017636684296702;
            aWeights( 18 ) =  0.000970017636684296702;

            aWeights( 19 ) =  0.0275573192239850917;
            aWeights( 20 ) =  0.0275573192239850917;
            aWeights( 21 ) =  0.0275573192239850917;
            aWeights( 22 ) =  0.0275573192239850917;
            aWeights( 23 ) =  0.0275573192239850917;
            aWeights( 24 ) =  0.0275573192239850917;

            aWeights( 25 ) =  0.0275573192239850917;
            aWeights( 26 ) =  0.0275573192239850917;
            aWeights( 27 ) =  0.0275573192239850917;
            aWeights( 28 ) =  0.0275573192239850917;
            aWeights( 29 ) =  0.0275573192239850917;
            aWeights( 30 ) =  0.0275573192239850917;

            aWeights( 0 ) = 0.0;
            aWeights( 0 ) = 1./6. - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET31_HPP
