//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PENTA16_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PENTA16_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        // source: 10.1016/j.compfluid.2013.01.002
        inline void
        gauss_penta16(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 16 );

            aPoints( 0, 0 ) =  0.33333333333333;
            aPoints( 1, 0 ) =  0.33333333333333;
            aPoints( 2, 0 ) =  0.00000000000000;

            aPoints( 0, 1 ) =  0.02540007089951;
            aPoints( 1, 1 ) =  0.48729996455025;
            aPoints( 2, 1 ) =  0.00000000000000;

            aPoints( 0, 2 ) =  0.48729996455025;
            aPoints( 1, 2 ) =  0.48729996455025;
            aPoints( 2, 2 ) =  0.00000000000000;

            aPoints( 0, 3 ) =  0.48729996455025;
            aPoints( 1, 3 ) =  0.02540007089951;
            aPoints( 2, 3 ) =  0.00000000000000;

            aPoints( 0, 4 ) =  0.10880379065926;
            aPoints( 1, 4 ) =  0.44559810467037;
            aPoints( 2, 4 ) = -0.87100293486544;

            aPoints( 0, 5 ) =  0.10880379065926;
            aPoints( 1, 5 ) =  0.44559810467037;
            aPoints( 2, 5 ) =  0.87100293486544;

            aPoints( 0, 6 ) =  0.44559810467037;
            aPoints( 1, 6 ) =  0.44559810467037;
            aPoints( 2, 6 ) = -0.87100293486544;

            aPoints( 0, 7 ) =  0.44559810467037;
            aPoints( 1, 7 ) =  0.44559810467037;
            aPoints( 2, 7 ) =  0.87100293486544;

            aPoints( 0, 8 ) =  0.44559810467037;
            aPoints( 1, 8 ) =  0.10880379065926;
            aPoints( 2, 8 ) = -0.87100293486544;

            aPoints( 0, 9 ) =  0.44559810467037;
            aPoints( 1, 9 ) =  0.10880379065926;
            aPoints( 2, 9 ) =  0.87100293486544;

            aPoints( 0, 10 ) =  0.79828210803458;
            aPoints( 1, 10 ) =  0.10085894598271;
            aPoints( 2, 10 ) = -0.57042698070516;

            aPoints( 0, 11 ) =  0.79828210803458;
            aPoints( 1, 11 ) =  0.10085894598271;
            aPoints( 2, 11 ) =  0.57042698070516;

            aPoints( 0, 12 ) =  0.10085894598271;
            aPoints( 1, 12 ) =  0.10085894598271;
            aPoints( 2, 12 ) = -0.57042698070516;

            aPoints( 0, 13 ) =  0.10085894598271;
            aPoints( 1, 13 ) =  0.10085894598271;
            aPoints( 2, 13 ) =  0.57042698070516;

            aPoints( 0, 14 ) =  0.10085894598271;
            aPoints( 1, 14 ) =  0.79828210803458;
            aPoints( 2, 14 ) = -0.57042698070516;

            aPoints( 0, 15 ) =  0.10085894598271;
            aPoints( 1, 15 ) =  0.79828210803458;
            aPoints( 2, 15 ) =  0.57042698070516;

            aWeights.set_size( 16 );
            aWeights( 0 ) =  0.17786388898287;
            aWeights( 1 ) =  0.05617751680707;
            aWeights( 2 ) =  0.05617751680707;
            aWeights( 3 ) =  0.05617751680707;
            aWeights( 4 ) =  0.04641535532904;
            aWeights( 5 ) =  0.04641535532904;
            aWeights( 6 ) =  0.04641535532904;
            aWeights( 7 ) =  0.04641535532904;
            aWeights( 8 ) =  0.04641535532904;
            aWeights( 9 ) =  0.04641535532904;
            aWeights( 10 ) =  0.06251857143695;
            aWeights( 11 ) =  0.06251857143695;
            aWeights( 12 ) =  0.06251857143695;
            aWeights( 13 ) =  0.06251857143695;
            aWeights( 14 ) =  0.06251857143695;
            aWeights( 15 ) =  0.06251857143695;
        }
// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PENTA16_HPP
