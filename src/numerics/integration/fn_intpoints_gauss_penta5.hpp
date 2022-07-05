//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PENTA5_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PENTA5_HPP

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
        gauss_penta5(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 5 );

            aPoints( 0, 0 ) =  1.000000000000000;
            aPoints( 1, 0 ) =  0.000000000000000;
            aPoints( 2, 0 ) =  0.000000000000000;

            aPoints( 0, 1 ) =  0.000000000000000;
            aPoints( 1, 1 ) =  0.000000000000000;
            aPoints( 2, 1 ) =  0.000000000000000;

            aPoints( 0, 2 ) =  0.000000000000000;
            aPoints( 1, 2 ) =  1.000000000000000;
            aPoints( 2, 2 ) =  0.000000000000000;

            aPoints( 0, 3 ) =  0.33333333333333333;
            aPoints( 1, 3 ) =  0.33333333333333333;
            aPoints( 2, 3 ) =  0.66666666666666667;

            aPoints( 0, 4 ) =  0.33333333333333333;
            aPoints( 1, 4 ) =  0.33333333333333333;
            aPoints( 2, 4 ) = -0.66666666666666667;

            aWeights.set_size( 5 );
            aWeights( 0 ) =  0.08333333333333333;
            aWeights( 1 ) =  0.08333333333333333;
            aWeights( 2 ) =  0.08333333333333333;
            aWeights( 3 ) =  0.375000000000000;
            aWeights( 4 ) =  0.375000000000000;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PENTA5_HPP
