//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET15_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET15_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        inline void
        gauss_tet15(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            // source 10.1016/0045-7825(86)90059-9

            real tAlpha1 = 1.0/3.0;
            real tBeta1  = 0.0;

            real tAlpha2 = 90.0/990.0;
            real tBeta2  = 1.0 - 3.0*tAlpha2;

            real tAlpha3 = 0.0665501535736642813;
            real tBeta3 = 0.5 - tAlpha3;

            aPoints.set_size( 4, 15 );

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
            aPoints( 2, 9 ) = tBeta3;

            aPoints( 0, 10 ) = tAlpha3;
            aPoints( 1, 10 ) = tBeta3;
            aPoints( 2, 10 ) = tAlpha3;

            aPoints( 0, 11 ) = tAlpha3;
            aPoints( 1, 11 ) = tBeta3;
            aPoints( 2, 11 ) = tBeta3;

            aPoints( 0, 12 ) = tBeta3;
            aPoints( 1, 12 ) = tAlpha3;
            aPoints( 2, 12 ) = tAlpha3;

            aPoints( 0, 13 ) = tBeta3;
            aPoints( 1, 13 ) = tAlpha3;
            aPoints( 2, 13 ) = tBeta3;

            aPoints( 0, 14 ) = tBeta3;
            aPoints( 1, 14 ) = tBeta3;
            aPoints( 2, 14 ) = tAlpha3;

            for( uint k=0; k<15; ++k )
            {
                aPoints( 3, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k )
                                   - aPoints( 2, k );
            }

            aWeights.set_size( 15 );

            aWeights(  0 ) = 0.0302836780970891856;

            aWeights(  1 ) = 0.00602678571428571597;
            aWeights(  2 ) = 0.00602678571428571597;
            aWeights(  3 ) = 0.00602678571428571597;
            aWeights(  4 ) = 0.00602678571428571597;

            aWeights(  5 ) = 0.0116452490860289742;
            aWeights(  6 ) = 0.0116452490860289742;
            aWeights(  7 ) = 0.0116452490860289742;
            aWeights(  8 ) = 0.0116452490860289742;

            aWeights(  9 ) = 0.0109491415613864534;
            aWeights( 10 ) = 0.0109491415613864534;
            aWeights( 11 ) = 0.0109491415613864534;
            aWeights( 12 ) = 0.0109491415613864534;
            aWeights( 13 ) = 0.0109491415613864534;
            aWeights( 14 ) = 0.0109491415613864534;

            aWeights( 0 ) = 0.0 ;
            aWeights( 0 ) = 1./6. - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET15_HPP
