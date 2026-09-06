/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_IF_PYRA14_HPP
#define BELFEM_CL_IF_PYRA14_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >
        ::element_type() const
        {
            return ElementType::PYRA14;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 14 );

            aXiHat( 0,  0 ) = -1.0;
            aXiHat( 1,  0 ) = -1.0;
            aXiHat( 2,  0 ) =  0.0;

            aXiHat( 0,  1 ) =  1.0;
            aXiHat( 1,  1 ) = -1.0;
            aXiHat( 2,  1 ) =  0.0;

            aXiHat( 0,  2 ) =  1.0;
            aXiHat( 1,  2 ) =  1.0;
            aXiHat( 2,  2 ) =  0.0;

            aXiHat( 0,  3 ) = -1.0;
            aXiHat( 1,  3 ) =  1.0;
            aXiHat( 2,  3 ) =  0.0;

            aXiHat( 0,  4 ) =  0.0;
            aXiHat( 1,  4 ) =  0.0;
            aXiHat( 2,  4 ) =  1.0;

            aXiHat( 0,  5 ) =  0.0 ;
            aXiHat( 1,  5 ) = -1.0 ;
            aXiHat( 2,  5 ) =  0.0 ;

            aXiHat( 0,  6 ) =  1.0 ;
            aXiHat( 1,  6 ) =  0.0 ;
            aXiHat( 2,  6 ) =  0.0 ;

            aXiHat( 0,  7 ) =  0.0 ;
            aXiHat( 1,  7 ) =  1.0 ;
            aXiHat( 2,  7 ) =  0.0 ;

            aXiHat( 0,  8 ) = -1.0 ;
            aXiHat( 1,  8 ) =  0.0 ;
            aXiHat( 2,  8 ) =  0.0 ;

            aXiHat( 0,  9 ) = -0.5 ;
            aXiHat( 1,  9 ) = -0.5 ;
            aXiHat( 2,  9 ) =  0.5 ;

            aXiHat( 0, 10 ) =  0.5 ;
            aXiHat( 1, 10 ) = -0.5 ;
            aXiHat( 2, 10 ) =  0.5 ;

            aXiHat( 0, 11 ) =  0.5 ;
            aXiHat( 1, 11 ) =  0.5 ;
            aXiHat( 2, 11 ) =  0.5 ;

            aXiHat( 0, 12 ) = -0.5 ;
            aXiHat( 1, 12 ) =  0.5 ;
            aXiHat( 2, 12 ) =  0.5 ;

            aXiHat( 0, 13 ) =  0.0;
            aXiHat( 1, 13 ) =  0.0;
            aXiHat( 2, 13 ) =  0.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >
        ::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            real xi2 = xi * xi;
            real eta2 = eta * eta;
            real zeta2 = zeta * zeta;
            real phi = eta * zeta;
            real psi = xi * zeta;
            real chi = xi * eta;
            real alpha = eta * xi2;
            real beta = xi * eta2;
            real gamma  = chi * zeta;
            real gamma2 = gamma + gamma ;

            real delta = xi2 * eta2;

            aN.set_size( 1, 14 );

            aN( 0,  0 ) = 0.25 * ( delta - alpha - beta + chi ) - 0.5 * gamma + 0.125 * ( phi + psi ) + 0.0625 * ( zeta2 - zeta );
            aN( 0,  1 ) = 0.25 * ( delta - alpha + beta - chi ) + 0.5 * gamma + 0.125 * ( phi - psi ) + 0.0625 * ( zeta2 - zeta );
            aN( 0,  2 ) = 0.25 * ( delta + alpha + beta + chi ) - 0.5 * gamma - 0.125 * ( phi + psi ) + 0.0625 * ( zeta2 - zeta );
            aN( 0,  3 ) = 0.25 * ( delta + alpha - beta - chi ) + 0.5 * gamma - 0.125 * ( phi - psi ) + 0.0625 * ( zeta2 - zeta );

            aN( 0,  4 ) =  zeta2 + zeta2 - zeta ;

            aN( 0,  5 ) = 0.5 * ( eta2 - eta + alpha - delta) + 0.75 * phi + 0.375 * ( zeta2 - zeta ) ;
            aN( 0,  6 ) = 0.5 * ( xi2 + xi - beta  - delta )  - 0.75 * psi + 0.375 * ( zeta2 - zeta ) ;
            aN( 0,  7 ) = 0.5 * ( eta2 + eta - alpha - delta) - 0.75 * phi + 0.375 * ( zeta2 - zeta ) ;
            aN( 0,  8 ) = 0.5 * ( xi2 - xi + beta - delta )   + 0.75 * psi + 0.375 * ( zeta2 - zeta ) ;

            aN( 0,  9 ) = zeta - phi - psi - zeta2 + gamma2 ;
            aN( 0, 10 ) = zeta - phi + psi - zeta2 - gamma2 ;
            aN( 0, 11 ) = zeta + phi + psi - zeta2 + gamma2 ;
            aN( 0, 12 ) = zeta + phi - psi - zeta2 - gamma2 ;

            aN( 0, 13 ) = delta - xi2-eta2 - 1.25*zeta + 0.25*zeta2 + 1.0 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >
        ::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real xi   = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            real xi2 = xi * xi;
            real eta2 = eta * eta;

            real phi = eta * zeta;
            real psi = xi * zeta;
            real chi = xi * eta;

            real phi2 = phi + phi ;
            real psi2 = psi + psi ;

            real alpha = eta * xi2;
            real beta = xi * eta2;
            real zeta8 = 0.125 * zeta ;

            adNdXi.set_size( 3, 14 );

            adNdXi( 0,  0 ) = 0.5 * ( beta - phi - chi )  + 0.25 * ( eta - eta2 ) + zeta8 ;
            adNdXi( 1,  0 ) = 0.5 * ( alpha - psi - chi ) + 0.25 * ( xi - xi2 )   + zeta8 ;
            adNdXi( 2,  0 ) = 0.125 * ( xi + eta + zeta ) - 0.5 * chi - 0.0625 ;

            adNdXi( 0,  1 ) = 0.5 * ( beta + phi - chi )  + 0.25 * ( eta2 - eta ) - zeta8 ;
            adNdXi( 1,  1 ) = 0.5 * ( alpha + psi + chi ) - 0.25 * ( xi + xi2)    + zeta8 ;
            adNdXi( 2,  1 ) = 0.125 * ( eta - xi + zeta ) + 0.5 * chi - 0.0625 ;

            adNdXi( 0,  2 ) = 0.5 * ( beta - phi + chi  ) + 0.25 * ( eta2 + eta ) - zeta8 ;
            adNdXi( 1,  2 ) = 0.5 * ( alpha - psi + chi ) + 0.25 * ( xi + xi2 )   - zeta8 ;
            adNdXi( 2,  2 ) = 0.125 * ( zeta - xi - eta ) - 0.5 * chi - 0.0625 ;

            adNdXi( 0,  3 ) = 0.5 * ( beta + phi + chi )  - 0.25 * ( eta + eta2 ) + zeta8 ;
            adNdXi( 1,  3 ) = 0.5 * ( alpha + psi - chi ) + 0.25 * ( xi2 - xi )   - zeta8 ;
            adNdXi( 2,  3 ) = 0.125 * ( xi - eta + zeta ) + 0.5 * chi - 0.0625 ;

            adNdXi( 0,  4 ) = 0.0 ;
            adNdXi( 1,  4 ) = 0.0 ;
            adNdXi( 2,  4 ) = 4.0 * zeta - 1.0 ;

            adNdXi( 0,  5 ) = chi - beta ;
            adNdXi( 1,  5 ) = eta - alpha + 0.75 * zeta + 0.5 * xi2 - 0.5 ;
            adNdXi( 2,  5 ) = 0.75 * ( zeta + eta ) - 0.375 ;

            adNdXi( 0,  6 ) = xi - beta - 0.75 * zeta - 0.5 * eta2 + 0.5 ;
            adNdXi( 1,  6 ) = -chi-alpha ;
            adNdXi( 2,  6 ) = 0.75 * ( zeta - xi )  - 0.375 ;

            adNdXi( 0,  7 ) = -chi-beta ;
            adNdXi( 1,  7 ) = eta - alpha - 0.75 * zeta - 0.5 * xi2 + 0.5 ;
            adNdXi( 2,  7 ) = 0.75 * ( zeta - eta ) - 0.375 ;

            adNdXi( 0,  8 ) = xi - beta + 0.75 * zeta + 0.5 * eta2 - 0.5 ;
            adNdXi( 1,  8 ) = chi - alpha ;
            adNdXi( 2,  8 ) = 0.75 * ( zeta + xi )  - 0.375 ;

            adNdXi( 0,  9 ) = phi2 - zeta ;
            adNdXi( 1,  9 ) = psi2 - zeta ;
            adNdXi( 2,  9 ) = 1.0 - xi - eta + 2.0 * ( chi - zeta ) ;

            adNdXi( 0, 10 ) =  zeta - phi2 ;
            adNdXi( 1, 10 ) = -zeta - psi2 ;
            adNdXi( 2, 10 ) =  1.0 + xi - eta - 2.0 * ( chi + zeta );

            adNdXi( 0, 11 ) =  zeta + phi2 ;
            adNdXi( 1, 11 ) =  zeta + psi2 ;
            adNdXi( 2, 11 ) =  1.0 + xi + eta + 2.0 * ( chi - zeta ) ;

            adNdXi( 0, 12 ) = -zeta - phi2 ;
            adNdXi( 1, 12 ) =  zeta - psi2 ;
            adNdXi( 2, 12 ) =  1.0 - xi + eta - 2.0 * ( chi + zeta ) ;

            adNdXi( 0, 13 ) = 2.0 * ( beta - xi );
            adNdXi( 1, 13 ) = 2.0 * ( alpha - eta );
            adNdXi( 2, 13 ) = 0.5 * zeta - 1.25 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            // rows: d2/dxi2, d2/deta2, d2/dzeta2, d2/(deta*dzeta), d2/(dxi*dzeta), d2/(dxi*deta)
            const real xi   = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            const real hxi    = 0.5 * xi ;
            const real heta   = 0.5 * eta ;
            const real hzeta  = 0.5 * zeta ;

            const real xi2    = xi * xi ;
            const real eta2   = eta * eta ;
            const real chi    = xi * eta ;


            ad2NdXi2.set_size( 6, 14, 0.0 );
            ad2NdXi2( 0,   0 ) = heta * ( eta - 1.0 ) ;
            ad2NdXi2( 1,   0 ) = hxi  * ( xi  - 1.0 ) ;
            ad2NdXi2( 2,   0 ) = 0.125 ;
            ad2NdXi2( 3,   0 ) = 0.125 - hxi ;
            ad2NdXi2( 4,   0 ) = 0.125 - heta ;
            ad2NdXi2( 5,   0 ) = chi - heta - hxi - hzeta + 0.25 ;

            ad2NdXi2( 0,   1 ) = heta * ( eta - 1.0 ) ;
            ad2NdXi2( 1,   1 ) = hxi  * ( xi  + 1.0 ) ;
            ad2NdXi2( 2,   1 ) = 0.125 ;
            ad2NdXi2( 3,   1 ) = hxi + 0.125 ;
            ad2NdXi2( 4,   1 ) = heta - 0.125 ;
            ad2NdXi2( 5,   1 ) = chi + heta - hxi + hzeta - 0.25 ;

            ad2NdXi2( 0,   2 ) = heta * ( eta + 1.0 ) ;
            ad2NdXi2( 1,   2 ) = hxi  * ( xi  + 1.0 ) ;
            ad2NdXi2( 2,   2 ) = 0.125 ;
            ad2NdXi2( 3,   2 ) = - hxi - 0.125 ;
            ad2NdXi2( 4,   2 ) = - heta - 0.125 ;
            ad2NdXi2( 5,   2 ) = chi + heta + hxi - hzeta + 0.25 ;

            ad2NdXi2( 0,   3 ) = heta * ( eta + 1.0 ) ;
            ad2NdXi2( 1,   3 ) = hxi  * ( xi  - 1.0 ) ;
            ad2NdXi2( 2,   3 ) = 0.125 ;
            ad2NdXi2( 3,   3 ) = hxi - 0.125 ;
            ad2NdXi2( 4,   3 ) = heta + 0.125 ;
            ad2NdXi2( 5,   3 ) = chi - heta + hxi + hzeta - 0.25 ;

            ad2NdXi2( 2,   4 ) = 4.0 ;

            ad2NdXi2( 0,   5 ) = eta * ( 1.0 - eta ) ;
            ad2NdXi2( 1,   5 ) = 1.0 - xi2 ;
            ad2NdXi2( 2,   5 ) = 0.75 ;
            ad2NdXi2( 3,   5 ) = 0.75 ;
            ad2NdXi2( 5,   5 ) = xi * ( 1.0 - eta - eta ) ;

            ad2NdXi2( 0,   6 ) = 1.0 - eta2 ;
            ad2NdXi2( 1,   6 ) = - xi * ( xi + 1.0 ) ;
            ad2NdXi2( 2,   6 ) = 0.75 ;
            ad2NdXi2( 4,   6 ) = -0.75 ;
            ad2NdXi2( 5,   6 ) = - eta * ( xi + xi + 1.0 ) ;

            ad2NdXi2( 0,   7 ) = - eta * ( eta + 1.0 ) ;
            ad2NdXi2( 1,   7 ) = 1.0 - xi2 ;
            ad2NdXi2( 2,   7 ) = 0.75 ;
            ad2NdXi2( 3,   7 ) = -0.75 ;
            ad2NdXi2( 5,   7 ) = - xi * ( eta + eta + 1.0 ) ;

            ad2NdXi2( 0,   8 ) = 1.0 - eta2 ;
            ad2NdXi2( 1,   8 ) = xi * ( 1.0 - xi ) ;
            ad2NdXi2( 2,   8 ) = 0.75 ;
            ad2NdXi2( 4,   8 ) = 0.75 ;
            ad2NdXi2( 5,   8 ) = eta * ( 1.0 - xi - xi ) ;

            ad2NdXi2( 2,   9 ) = -2.0 ;
            ad2NdXi2( 3,   9 ) = xi + xi - 1.0 ;
            ad2NdXi2( 4,   9 ) = eta + eta - 1.0 ;
            ad2NdXi2( 5,   9 ) = zeta + zeta ;

            ad2NdXi2( 2,  10 ) = -2.0 ;
            ad2NdXi2( 3,  10 ) = - xi - xi - 1.0 ;
            ad2NdXi2( 4,  10 ) = 1.0 - eta - eta ;
            ad2NdXi2( 5,  10 ) = - zeta - zeta ;

            ad2NdXi2( 2,  11 ) = -2.0 ;
            ad2NdXi2( 3,  11 ) = xi + xi + 1.0 ;
            ad2NdXi2( 4,  11 ) = eta + eta + 1.0 ;
            ad2NdXi2( 5,  11 ) = zeta + zeta ;

            ad2NdXi2( 2,  12 ) = -2.0 ;
            ad2NdXi2( 3,  12 ) = 1.0 - xi - xi ;
            ad2NdXi2( 4,  12 ) = -1.0 - eta - eta ;
            ad2NdXi2( 5,  12 ) = - zeta - zeta ;

            ad2NdXi2( 0,  13 ) = eta2 + eta2 - 2.0 ;
            ad2NdXi2( 1,  13 ) = xi2 + xi2 - 2.0 ;
            ad2NdXi2( 2,  13 ) = 0.5 ;
            ad2NdXi2( 5,  13 ) = 4.0 * chi ;
        }

//------------------------------------------------------------------------------
    }
}


#endif //BELFEM_CL_IF_PYRA14_HPP
