// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI7_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI7_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 5th order interpolation
         *
         * P.C. Hammer, O.J. Marlowe, A.H. Stroud :
         * Numerical Integration Over Simplexes and Cones
         * Mathematical Tables and Other Aids to Computation, vol. 10, no. 55, pp. 130-137, 1956
         * https://doi.org/10.1090/S0025-5718-1956-0086389-6
         */
        inline void
        gauss_tri7(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 7 );

            aPoints( 0, 0 ) = 1./3.;
            aPoints( 1, 0 ) = 1./3.;

            aPoints( 0, 1 ) = 0.10128650732345634;
            aPoints( 1, 1 ) = 0.10128650732345634;

            aPoints( 0, 2 ) = 0.4701420641051151;
            aPoints( 1, 2 ) = 0.4701420641051151;

            aPoints( 0, 3 ) = 0.10128650732345634;
            aPoints( 1, 3 ) = 0.7974269853530873;

            aPoints( 0, 4 ) = 0.4701420641051151;
            aPoints( 1, 4 ) = 0.05971587178976982;

            aPoints( 0, 5 ) = 0.7974269853530873;
            aPoints( 1, 5 ) = 0.10128650732345634;

            aPoints( 0, 6 ) = 0.05971587178976982;
            aPoints( 1, 6 ) = 0.4701420641051151;

            for( uint k=0; k<7; ++k )
            {
                aPoints( 2, k ) =  1.0
                        - aPoints( 0, k )
                        - aPoints( 1, k );
            }

            aWeights.set_size( 7 );

            aWeights( 0 ) = 0.1125;
            aWeights( 0 ) = 0.0 ;
            aWeights( 1 ) = 0.06296959027241358;
            aWeights( 2 ) = 0.0661970763942531;
            aWeights( 3 ) = 0.06296959027241358;
            aWeights( 4 ) = 0.0661970763942531;
            aWeights( 5 ) = 0.06296959027241358;
            aWeights( 6 ) = 0.0661970763942531;

            aWeights( 0 ) = 0.0;
            aWeights( 0 ) = 0.5 - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI7_HPP
