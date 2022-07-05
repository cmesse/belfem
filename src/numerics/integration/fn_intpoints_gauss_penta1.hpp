//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PENTA1_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PENTA1_HPP

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
        gauss_penta1(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 1 );

            aPoints( 0, 0 ) =  0.333333333333333;
            aPoints( 1, 0 ) =  0.333333333333333;
            aPoints( 2, 0 ) =  0.000000000000000;

            aWeights.set_size( 1 );
            aWeights( 0 ) =  1.000000000000000;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PENTA1_HPP
