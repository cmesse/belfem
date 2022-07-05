//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PENTA11_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PENTA11_HPP

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
        gauss_penta11(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 11 );

            aPoints( 0, 0 ) =  0.06268838027601;
            aPoints( 1, 0 ) =  0.46865580986200;
            aPoints( 2, 0 ) =  0.00000000000000;

            aPoints( 0, 1 ) =  0.46865580986200;
            aPoints( 1, 1 ) =  0.46865580986200;
            aPoints( 2, 1 ) =  0.00000000000000;

            aPoints( 0, 2 ) =  0.46865580986200;
            aPoints( 1, 2 ) =  0.06268838027601;
            aPoints( 2, 2 ) =  0.00000000000000;

            aPoints( 0, 3 ) =  0.33333333333333;
            aPoints( 1, 3 ) =  0.33333333333333;
            aPoints( 2, 3 ) = -0.86686197400903;

            aPoints( 0, 4 ) =  0.33333333333333;
            aPoints( 1, 4 ) =  0.33333333333333;
            aPoints( 2, 4 ) =  0.86686197400903;

            aPoints( 0, 5 ) =  0.79851918840218;
            aPoints( 1, 5 ) =  0.10074040579891;
            aPoints( 2, 5 ) = -0.67563982368227;

            aPoints( 0, 6 ) =  0.79851918840218;
            aPoints( 1, 6 ) =  0.10074040579891;
            aPoints( 2, 6 ) =  0.67563982368227;

            aPoints( 0, 7 ) =  0.10074040579891;
            aPoints( 1, 7 ) =  0.10074040579891;
            aPoints( 2, 7 ) = -0.67563982368227;

            aPoints( 0, 8 ) =  0.10074040579891;
            aPoints( 1, 8 ) =  0.10074040579891;
            aPoints( 2, 8 ) =  0.67563982368227;

            aPoints( 0, 9 ) =  0.10074040579891;
            aPoints( 1, 9 ) =  0.79851918840218;
            aPoints( 2, 9 ) = -0.67563982368227;

            aPoints( 0, 10 ) =  0.10074040579891;
            aPoints( 1, 10 ) =  0.79851918840218;
            aPoints( 2, 10 ) =  0.67563982368227;

            aWeights.set_size( 11 );
            aWeights( 0 ) =  0.13641461260548;
            aWeights( 1 ) =  0.13641461260548;
            aWeights( 2 ) =  0.13641461260548;
            aWeights( 3 ) =  0.10791197481553;
            aWeights( 4 ) =  0.10791197481553;
            aWeights( 5 ) =  0.06248870209208;
            aWeights( 6 ) =  0.06248870209208;
            aWeights( 7 ) =  0.06248870209208;
            aWeights( 8 ) =  0.06248870209208;
            aWeights( 9 ) =  0.06248870209208;
            aWeights( 10 ) =  0.06248870209208;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PENTA11_HPP
