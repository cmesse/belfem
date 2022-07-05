// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI19_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI19_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 9th order interpolation
         *
         * J.N. Lyness, D. Jespersen :
         * Moderate Degree Symmetric Quadrature Rules for the Triangle
         * IMA Journal of Applied Mathematics, vol. 15, no. 1, pp. 19-32, 1975
         * https://doi.org/10.1093/imamat/15.1.19
         */
        inline void
        gauss_tri19(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 2, 19 );

            aPoints( 0, 0 ) = 0.3333333333333333;
            aPoints( 1, 0 ) = 0.3333333333333333;

            aPoints( 0, 1 ) = 0.489682519198737;
            aPoints( 1, 1 ) = 0.489682519198737;

            aPoints( 0, 2 ) = 0.4370895914929355;
            aPoints( 1, 2 ) = 0.4370895914929355;

            aPoints( 0, 3 ) = 0.1882035356190322;
            aPoints( 1, 3 ) = 0.1882035356190322;

            aPoints( 0, 4 ) = 0.04472951339445297;
            aPoints( 1, 4 ) = 0.04472951339445297;

            aPoints( 0, 5 ) = 0.489682519198737;
            aPoints( 1, 5 ) = 0.02063496160252598;

            aPoints( 0, 6 ) = 0.4370895914929355;
            aPoints( 1, 6 ) = 0.12582081701412895;

            aPoints( 0, 7 ) = 0.1882035356190322;
            aPoints( 1, 7 ) = 0.6235929287619356;

            aPoints( 0, 8 ) = 0.04472951339445297;
            aPoints( 1, 8 ) = 0.9105409732110941;

            aPoints( 0, 9 ) = 0.02063496160252598;
            aPoints( 1, 9 ) = 0.489682519198737;

            aPoints( 0, 10 ) = 0.12582081701412895;
            aPoints( 1, 10 ) = 0.4370895914929355;

            aPoints( 0, 11 ) = 0.6235929287619356;
            aPoints( 1, 11 ) = 0.1882035356190322;

            aPoints( 0, 12 ) = 0.9105409732110941;
            aPoints( 1, 12 ) = 0.04472951339445297;

            aPoints( 0, 13 ) = 0.03683841205473626;
            aPoints( 1, 13 ) = 0.741198598784498;

            aPoints( 0, 14 ) = 0.22196298916076573;
            aPoints( 1, 14 ) = 0.03683841205473626;

            aPoints( 0, 15 ) = 0.741198598784498;
            aPoints( 1, 15 ) = 0.22196298916076573;

            aPoints( 0, 16 ) = 0.741198598784498;
            aPoints( 1, 16 ) = 0.03683841205473626;

            aPoints( 0, 17 ) = 0.22196298916076573;
            aPoints( 1, 17 ) = 0.741198598784498;

            aPoints( 0, 18 ) = 0.03683841205473626;
            aPoints( 1, 18 ) = 0.22196298916076573;

            aWeights.set_size( 19 );

            aWeights( 0 ) = 0.04856789814139805;
            aWeights( 1 ) = 0.015667350113569917;
            aWeights( 2 ) = 0.038913770502387715;
            aWeights( 3 ) = 0.03982386946360452;
            aWeights( 4 ) = 0.01278883782934905;
            aWeights( 5 ) = 0.015667350113569917;
            aWeights( 6 ) = 0.038913770502387715;
            aWeights( 7 ) = 0.03982386946360452;
            aWeights( 8 ) = 0.01278883782934905;
            aWeights( 9 ) = 0.015667350113569917;
            aWeights( 10 ) = 0.038913770502387715;
            aWeights( 11 ) = 0.03982386946360452;
            aWeights( 12 ) = 0.01278883782934905;
            aWeights( 13 ) = 0.021641769688644702;
            aWeights( 14 ) = 0.021641769688644702;
            aWeights( 15 ) = 0.021641769688644702;
            aWeights( 16 ) = 0.021641769688644702;
            aWeights( 17 ) = 0.021641769688644702;
            aWeights( 18 ) = 0.021641769688644702;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI19_HPP
