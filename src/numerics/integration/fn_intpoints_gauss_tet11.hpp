//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET11_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET11_HPP

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
        gauss_tet11(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            // source Keast: 10.1016/0045-7825(86)90059-9

            aPoints.set_size( 4, 11 );

            real tAlpha1 = 1.0/14.0;
            real tBeta1  = 1.0 - 3.0*tAlpha1;

            real tAlpha2 = 0.25*( 1.0 - std::sqrt( 5.0/14.0 ) );
            real tBeta2  = 0.5 - tAlpha2;

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
            aPoints( 2, 5 ) = tBeta2;

            aPoints( 0, 6 ) = tAlpha2;
            aPoints( 1, 6 ) = tBeta2;
            aPoints( 2, 6 ) = tAlpha2;

            aPoints( 0, 7 ) = tAlpha2;
            aPoints( 1, 7 ) = tBeta2;
            aPoints( 2, 7 ) = tBeta2;

            aPoints( 0, 8 ) = tBeta2;
            aPoints( 1, 8 ) = tAlpha2;
            aPoints( 2, 8 ) = tAlpha2;

            aPoints( 0, 9 ) = tBeta2;
            aPoints( 1, 9 ) = tAlpha2;
            aPoints( 2, 9 ) = tBeta2;

            aPoints( 0, 10 ) = tBeta2;
            aPoints( 1, 10 ) = tBeta2;
            aPoints( 2, 10 ) = tAlpha2;

            for( uint k=0; k<11; ++k )
            {
                aPoints( 3, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k )
                                   - aPoints( 2, k );
            }


            aWeights.set_size( 11 );

            aWeights(  0 ) = -74.0/5625.0;
            aWeights(  1 ) = 343.0/45000.0;
            aWeights(  2 ) = aWeights( 1 );
            aWeights(  3 ) = aWeights( 1 );
            aWeights(  4 ) = aWeights( 1 );
            aWeights(  5 ) = 56.0/2250.0;
            aWeights(  6 ) = aWeights( 5 );
            aWeights(  7 ) = aWeights( 5 );
            aWeights(  8 ) = aWeights( 5 );
            aWeights(  9 ) = aWeights( 5 );
            aWeights( 10 ) = aWeights( 5 );

            aWeights( 0 ) = 0.0 ;
            aWeights( 0 ) = 1./6. - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET11_HPP
