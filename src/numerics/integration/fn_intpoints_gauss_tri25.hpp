// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI25_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI25_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 10th order interpolation
         *
         * F.D. Witherden, P.E. Vincent :
         * On the identification of symmetric quadrature rules for finite element methods
         * Computers & Mathematics with Applications, vol. 69, no. 10, pp. 1232–1241, 2015
         * https://doi.org/10.1016/j.camwa.2015.03.017
         */
        inline void
        gauss_tri25(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 2, 25 );

            aPoints( 0, 0 ) = 0.3333333333333333;
            aPoints( 1, 0 ) = 0.3333333333333333;

            aPoints( 0, 1 ) = 0.03205537321694352;
            aPoints( 1, 1 ) = 0.03205537321694352;

            aPoints( 0, 2 ) = 0.14216110105656438;
            aPoints( 1, 2 ) = 0.14216110105656438;

            aPoints( 0, 3 ) = 0.03205537321694352;
            aPoints( 1, 3 ) = 0.935889253566113;

            aPoints( 0, 4 ) = 0.14216110105656438;
            aPoints( 1, 4 ) = 0.7156777978868712;

            aPoints( 0, 5 ) = 0.935889253566113;
            aPoints( 1, 5 ) = 0.03205537321694352;

            aPoints( 0, 6 ) = 0.7156777978868712;
            aPoints( 1, 6 ) = 0.14216110105656438;

            aPoints( 0, 7 ) = 0.028367665339938453;
            aPoints( 1, 7 ) = 0.1637017337371825;

            aPoints( 0, 8 ) = 0.029619889488729734;
            aPoints( 1, 8 ) = 0.369146781827811;

            aPoints( 0, 9 ) = 0.14813288578382056;
            aPoints( 1, 9 ) = 0.32181299528883545;

            aPoints( 0, 10 ) = 0.807930600922879;
            aPoints( 1, 10 ) = 0.028367665339938453;

            aPoints( 0, 11 ) = 0.6012333286834592;
            aPoints( 1, 11 ) = 0.029619889488729734;

            aPoints( 0, 12 ) = 0.530054118927344;
            aPoints( 1, 12 ) = 0.14813288578382056;

            aPoints( 0, 13 ) = 0.1637017337371825;
            aPoints( 1, 13 ) = 0.807930600922879;

            aPoints( 0, 14 ) = 0.369146781827811;
            aPoints( 1, 14 ) = 0.6012333286834592;

            aPoints( 0, 15 ) = 0.32181299528883545;
            aPoints( 1, 15 ) = 0.530054118927344;

            aPoints( 0, 16 ) = 0.1637017337371825;
            aPoints( 1, 16 ) = 0.028367665339938453;

            aPoints( 0, 17 ) = 0.369146781827811;
            aPoints( 1, 17 ) = 0.029619889488729734;

            aPoints( 0, 18 ) = 0.32181299528883545;
            aPoints( 1, 18 ) = 0.14813288578382056;

            aPoints( 0, 19 ) = 0.807930600922879;
            aPoints( 1, 19 ) = 0.1637017337371825;

            aPoints( 0, 20 ) = 0.6012333286834592;
            aPoints( 1, 20 ) = 0.369146781827811;

            aPoints( 0, 21 ) = 0.530054118927344;
            aPoints( 1, 21 ) = 0.32181299528883545;

            aPoints( 0, 22 ) = 0.028367665339938453;
            aPoints( 1, 22 ) = 0.807930600922879;

            aPoints( 0, 23 ) = 0.029619889488729734;
            aPoints( 1, 23 ) = 0.6012333286834592;

            aPoints( 0, 24 ) = 0.14813288578382056;
            aPoints( 1, 24 ) = 0.530054118927344;

            aWeights.set_size( 25 );

            aWeights( 0 ) = 0.040871664573142986;
            aWeights( 1 ) = 0.006676484406574783;
            aWeights( 2 ) = 0.022978981802372365;
            aWeights( 3 ) = 0.006676484406574783;
            aWeights( 4 ) = 0.022978981802372365;
            aWeights( 5 ) = 0.006676484406574783;
            aWeights( 6 ) = 0.022978981802372365;
            aWeights( 7 ) = 0.012648878853644192;
            aWeights( 8 ) = 0.017092324081479714;
            aWeights( 9 ) = 0.03195245319821202;
            aWeights( 10 ) = 0.012648878853644192;
            aWeights( 11 ) = 0.017092324081479714;
            aWeights( 12 ) = 0.03195245319821202;
            aWeights( 13 ) = 0.012648878853644192;
            aWeights( 14 ) = 0.017092324081479714;
            aWeights( 15 ) = 0.03195245319821202;
            aWeights( 16 ) = 0.012648878853644192;
            aWeights( 17 ) = 0.017092324081479714;
            aWeights( 18 ) = 0.03195245319821202;
            aWeights( 19 ) = 0.012648878853644192;
            aWeights( 20 ) = 0.017092324081479714;
            aWeights( 21 ) = 0.03195245319821202;
            aWeights( 22 ) = 0.012648878853644192;
            aWeights( 23 ) = 0.017092324081479714;
            aWeights( 24 ) = 0.03195245319821202;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI25_HPP
