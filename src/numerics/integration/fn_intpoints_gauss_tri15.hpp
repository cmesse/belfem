// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI15_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI15_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 7th order interpolation
         *
         * M.E. Laursen, M. Gellert :
         * Some criteria for numerically integrated matrices and quadrature formulas for triangles
         * International Journal for Numerical Methods in Engineering, vol. 12, no. 1, pp. 67–76, 1978
         * https://doi.org/10.1002/nme.1620120107
         */
        inline void
        gauss_tri15(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 15 );

            aPoints( 0, 0 ) = 0.06493051315916486;
            aPoints( 1, 0 ) = 0.06493051315916486;

            aPoints( 0, 1 ) = 0.06493051315916486;
            aPoints( 1, 1 ) = 0.8701389736816703;

            aPoints( 0, 2 ) = 0.8701389736816703;
            aPoints( 1, 2 ) = 0.06493051315916486;

            aPoints( 0, 3 ) = 0.2845755842491703;
            aPoints( 1, 3 ) = 0.517039939069323;

            aPoints( 0, 4 ) = 0.3135591843849315;
            aPoints( 1, 4 ) = 0.04386347179237249;

            aPoints( 0, 5 ) = 0.19838447668150672;
            aPoints( 1, 5 ) = 0.2845755842491703;

            aPoints( 0, 6 ) = 0.642577343822696;
            aPoints( 1, 6 ) = 0.3135591843849315;

            aPoints( 0, 7 ) = 0.517039939069323;
            aPoints( 1, 7 ) = 0.19838447668150672;

            aPoints( 0, 8 ) = 0.04386347179237249;
            aPoints( 1, 8 ) = 0.642577343822696;

            aPoints( 0, 9 ) = 0.517039939069323;
            aPoints( 1, 9 ) = 0.2845755842491703;

            aPoints( 0, 10 ) = 0.04386347179237249;
            aPoints( 1, 10 ) = 0.3135591843849315;

            aPoints( 0, 11 ) = 0.19838447668150672;
            aPoints( 1, 11 ) = 0.517039939069323;

            aPoints( 0, 12 ) = 0.642577343822696;
            aPoints( 1, 12 ) = 0.04386347179237249;

            aPoints( 0, 13 ) = 0.2845755842491703;
            aPoints( 1, 13 ) = 0.19838447668150672;

            aPoints( 0, 14 ) = 0.3135591843849315;
            aPoints( 1, 14 ) = 0.642577343822696;

            for( uint k=0; k<15; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 15 );

            aWeights( 0 ) = 0.026538900895116242;
            aWeights( 1 ) = 0.026538900895116242;
            aWeights( 2 ) = 0.026538900895116242;
            aWeights( 3 ) = 0.03542654184606683;
            aWeights( 4 ) = 0.0346373410397085;
            aWeights( 5 ) = 0.03542654184606683;
            aWeights( 6 ) = 0.0346373410397085;
            aWeights( 7 ) = 0.03542654184606683;
            aWeights( 8 ) = 0.0346373410397085;
            aWeights( 9 ) = 0.03542654184606683;
            aWeights( 10 ) = 0.0346373410397085;
            aWeights( 11 ) = 0.03542654184606683;
            aWeights( 12 ) = 0.0346373410397085;
            aWeights( 13 ) = 0.03542654184606683;
            aWeights( 14 ) = 0.0346373410397085;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI15_HPP
