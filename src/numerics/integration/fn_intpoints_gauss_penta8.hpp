//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PENTA8_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PENTA8_HPP

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
        gauss_penta8(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 8 );

            aPoints( 0, 0 ) =  0.72496028726946;
            aPoints( 1, 0 ) =  0.23737598643765;
            aPoints( 2, 0 ) =  0.00000000000000;

            aPoints( 0, 1 ) =  0.23737598643765;
            aPoints( 1, 1 ) =  0.72496028726946;
            aPoints( 2, 1 ) =  0.00000000000000;

            aPoints( 0, 2 ) =  0.03766372629289;
            aPoints( 1, 2 ) =  0.72496028726946;
            aPoints( 2, 2 ) =  0.00000000000000;

            aPoints( 0, 3 ) =  0.03766372629289;
            aPoints( 1, 3 ) =  0.23737598643765;
            aPoints( 2, 3 ) =  0.00000000000000;

            aPoints( 0, 4 ) =  0.23737598643765;
            aPoints( 1, 4 ) =  0.03766372629289;
            aPoints( 2, 4 ) =  0.00000000000000;

            aPoints( 0, 5 ) =  0.72496028726946;
            aPoints( 1, 5 ) =  0.03766372629289;
            aPoints( 2, 5 ) =  0.00000000000000;

            aPoints( 0, 6 ) =  0.33333333333333;
            aPoints( 1, 6 ) =  0.33333333333333;
            aPoints( 2, 6 ) =  1.00000000000000;

            aPoints( 0, 7 ) =  0.33333333333333;
            aPoints( 1, 7 ) =  0.33333333333333;
            aPoints( 2, 7 ) = -1.00000000000000;

            aWeights.set_size( 8 );
            aWeights( 0 ) =  0.11111111111111;
            aWeights( 1 ) =  0.11111111111111;
            aWeights( 2 ) =  0.11111111111111;
            aWeights( 3 ) =  0.11111111111111;
            aWeights( 4 ) =  0.11111111111111;
            aWeights( 5 ) =  0.11111111111111;
            aWeights( 6 ) =  0.16666666666667;
            aWeights( 7 ) =  0.16666666666667;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */
#endif //BELFEM_FN_INTPOINTS_GAUSS_PENTA8_HPP
