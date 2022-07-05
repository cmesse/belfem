//
// Created by Christian Messe on 2019-01-15.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI4_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI4_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        inline void
        gauss_tri4(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            // source : 10.1002/nme.1620070316
            aPoints.set_size( 2, 4 );

            aPoints( 0, 0 ) = 1.0 / 3.0;
            aPoints( 1, 0 ) = aPoints( 0, 0 );

            aPoints( 0, 1 ) = 0.6;
            aPoints( 1, 1 ) = 0.2;

            aPoints( 0, 2 ) = 0.2;
            aPoints( 1, 2 ) = 0.6;

            aPoints( 0, 3 ) = 0.2;
            aPoints( 1, 3 ) = 0.2;

            aWeights.set_size( 4 );

            aWeights( 0 ) = -27.0 / 96.0;
            aWeights( 1 ) = 25.0 / 96.0;
            aWeights( 2 ) = aWeights( 1 );
            aWeights( 3 ) = aWeights( 1 );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TRI3_HPP
