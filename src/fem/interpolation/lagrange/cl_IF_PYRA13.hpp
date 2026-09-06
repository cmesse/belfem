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

#ifndef BELFEM_CL_IF_PYRA13_HPP
#define BELFEM_CL_IF_PYRA13_HPP
#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >
        ::interpolation_order() const
        {
            return InterpolationOrder::SERENDIPITY ;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >
        ::element_type() const
        {
            return ElementType::PYRA13;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 13 );

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

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >
        ::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real xi = aXi( 0 );
            const real eta = aXi( 1 );
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

            aN.set_size( 1, 13 );

            aN( 0,  0 ) = 0.25 * ( zeta + chi - alpha - beta + xi2 + eta2 ) + 0.125 * ( phi + psi ) - 0.5*gamma - 0.25 ;
            aN( 0,  1 ) = 0.25 * ( zeta - chi - alpha + beta + xi2 + eta2 ) + 0.125 * ( phi - psi ) + 0.5*gamma - 0.25 ;
            aN( 0,  2 ) = 0.25 * ( zeta + chi + alpha + beta + xi2 + eta2 ) - 0.125 * ( phi + psi ) - 0.5*gamma - 0.25 ;
            aN( 0,  3 ) = 0.25 * ( zeta - chi + alpha - beta + xi2 + eta2 ) - 0.125 * ( phi - psi ) + 0.5*gamma - 0.25 ;

            aN( 0,  4 ) = zeta2 + zeta2 - zeta ;

            aN( 0,  5 ) =  0.5*( zeta2 - xi2 + alpha - eta ) + 0.75 * phi - zeta + 0.5 ;
            aN( 0,  6 ) = 0.5*( zeta2 - eta2 - beta + xi  ) - 0.75 * psi - zeta + 0.5 ;
            aN( 0,  7 ) = 0.5*( zeta2 - xi2 - alpha + eta ) - 0.75 * phi - zeta + 0.5 ;
            aN( 0,  8 ) = 0.5*( zeta2 - eta2 + beta - xi )  + 0.75 * psi - zeta + 0.5 ;

            aN( 0,  9 ) = zeta - phi - psi - zeta2 + gamma2 ;
            aN( 0, 10 ) = zeta - phi + psi - zeta2 - gamma2 ;
            aN( 0, 11 ) = zeta + phi + psi - zeta2 + gamma2 ;
            aN( 0, 12 ) = zeta + phi - psi - zeta2 - gamma2 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >
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
            real zeta8 = 0.125 * zeta ;

            adNdXi.set_size( 3, 13 );

            adNdXi( 0,  0 ) = 0.5 * ( xi - chi - phi ) + 0.25 * ( eta - eta2 ) + zeta8 ;
            adNdXi( 1,  0 ) = 0.5 * ( eta - chi - psi ) + 0.25 * ( xi - xi2 )  + zeta8 ;
            adNdXi( 2,  0 ) = 0.25 + 0.125 * ( xi + eta ) - 0.5 * chi ;

            adNdXi( 0,  1 ) = 0.5 * ( xi - chi + phi ) + 0.25 * ( eta2 - eta ) - zeta8 ;
            adNdXi( 1,  1 ) = 0.5 * ( eta + chi + psi ) - 0.25 * ( xi + xi2 )  + zeta8 ;
            adNdXi( 2,  1 ) = 0.25 + 0.125 * ( eta - xi ) + 0.5 * chi ;

            adNdXi( 0,  2 ) = 0.5 * ( xi + chi - phi )  + 0.25 * ( eta2 + eta ) - zeta8 ;
            adNdXi( 1,  2 ) = 0.5 * ( eta + chi - psi ) + 0.25 * ( xi + xi2 )   - zeta8 ;
            adNdXi( 2,  2 ) = 0.25 - 0.125 * ( xi + eta ) - 0.5 * chi ;

            adNdXi( 0,  3 ) = 0.5 * ( xi + chi + phi )  - 0.25 * ( eta2 + eta ) + zeta8 ;
            adNdXi( 1,  3 ) = 0.5 * ( eta - chi + psi ) + 0.25 * ( xi2 - xi )   - zeta8 ;
            adNdXi( 2,  3 ) = 0.25 + 0.125 * ( xi - eta ) + 0.5 * chi ;

            adNdXi( 0,  4 ) = 0.0 ;
            adNdXi( 1,  4 ) = 0.0 ;
            adNdXi( 2,  4 ) = 4.0 * zeta - 1.0 ;

            adNdXi( 0,  5 ) = chi - xi ;
            adNdXi( 1,  5 ) = 0.5 * xi2 + 0.75 *  zeta - 0.5 ;
            adNdXi( 2,  5 ) = zeta + 0.75 * eta - 1.0 ;

            adNdXi( 0,  6 ) = 0.5 - 0.5  * eta2 - 0.75 * zeta ;
            adNdXi( 1,  6 ) = -eta - chi ;
            adNdXi( 2,  6 ) = zeta - 0.75 * xi - 1.0 ;

            adNdXi( 0,  7 ) = -chi - xi ;
            adNdXi( 1,  7 ) = 0.5 - 0.5 * xi2 - 0.75 * zeta ;
            adNdXi( 2,  7 ) = zeta - 0.75 * eta - 1.0 ;

            adNdXi( 0,  8 ) = 0.5 * eta2 + 0.75 * zeta - 0.5 ;
            adNdXi( 1,  8 ) = chi - eta ;
            adNdXi( 2,  8 ) = zeta + 0.75 * xi - 1.0 ;

            adNdXi( 0,  9 ) = phi2 - zeta ;
            adNdXi( 1,  9 ) = psi2 - zeta ;
            adNdXi( 2,  9 ) = 1.0 - xi - eta + 2.0 * (chi - zeta);

            adNdXi( 0, 10 ) =  zeta - phi2 ;
            adNdXi( 1, 10 ) = -psi2 - zeta ;
            adNdXi( 2, 10 ) = 1.0 + xi - eta - 2.0 * (chi + zeta);

            adNdXi( 0, 11 ) = zeta + phi2 ;
            adNdXi( 1, 11 ) = zeta + psi2 ;
            adNdXi( 2, 11 ) = 1.0 + xi + eta + 2.0 * (chi - zeta);

            adNdXi( 0, 12 ) = -phi2 - zeta ;
            adNdXi( 1, 12 ) = zeta - psi2 ;
            adNdXi( 2, 12 ) = 1.0 - xi + eta - 2.0 * (chi + zeta);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {

            const real xi   = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            const real hxi    = 0.5 * xi ;
            const real heta   = 0.5 * eta ;
            const real hzeta  = 0.5 * zeta ;

            const real xi2    = xi + xi ;
            const real eta2   = eta + eta ;
            const real zeta2  = zeta + zeta ;

            ad2NdXi2.set_size( 6, 13, 0.0 );

            ad2NdXi2( 0,   0 ) = 0.5 - heta ;
            ad2NdXi2( 1,   0 ) = 0.5 - hxi ;
            ad2NdXi2( 2,   0 ) = 0.0 ;
            ad2NdXi2( 3,   0 ) = 0.125 - hxi ;
            ad2NdXi2( 4,   0 ) = 0.125 - heta ;
            ad2NdXi2( 5,   0 ) = 0.25 - hxi - hzeta - heta ;

            ad2NdXi2( 0,   1 ) = 0.5 - heta ;
            ad2NdXi2( 1,   1 ) = hxi + 0.5 ;
            ad2NdXi2( 2,   1 ) = 0.0 ;
            ad2NdXi2( 3,   1 ) = hxi + 0.125 ;
            ad2NdXi2( 4,   1 ) = heta - 0.125 ;
            ad2NdXi2( 5,   1 ) = heta - hxi + hzeta - 0.25 ;

            ad2NdXi2( 0,   2 ) = heta + 0.5 ;
            ad2NdXi2( 1,   2 ) = hxi + 0.5 ;
            ad2NdXi2( 2,   2 ) = 0.0 ;
            ad2NdXi2( 3,   2 ) =  - hxi - 0.125 ;
            ad2NdXi2( 4,   2 ) =  - heta - 0.125 ;
            ad2NdXi2( 5,   2 ) = heta + hxi - hzeta + 0.25 ;

            ad2NdXi2( 0,   3 ) = heta + 0.5 ;
            ad2NdXi2( 1,   3 ) = 0.5 - hxi ;
            ad2NdXi2( 2,   3 ) = 0.0 ;
            ad2NdXi2( 3,   3 ) = hxi - 0.125 ;
            ad2NdXi2( 4,   3 ) = heta + 0.125 ;
            ad2NdXi2( 5,   3 ) = hxi - heta + hzeta - 0.25 ;

            ad2NdXi2( 0,   4 ) = 0.0 ;
            ad2NdXi2( 1,   4 ) = 0.0 ;
            ad2NdXi2( 2,   4 ) = 4.0 ;
            ad2NdXi2( 3,   4 ) = 0.0 ;
            ad2NdXi2( 4,   4 ) = 0.0 ;
            ad2NdXi2( 5,   4 ) = 0.0 ;

            ad2NdXi2( 0,   5 ) = eta - 1.0 ;
            ad2NdXi2( 1,   5 ) = 0.0 ;
            ad2NdXi2( 2,   5 ) = 1.0 ;
            ad2NdXi2( 3,   5 ) = 0.75 ;
            ad2NdXi2( 4,   5 ) = 0.0 ;
            ad2NdXi2( 5,   5 ) = xi ;

            ad2NdXi2( 0,   6 ) = 0.0 ;
            ad2NdXi2( 1,   6 ) =  - 1.0 - xi ;
            ad2NdXi2( 2,   6 ) = 1.0 ;
            ad2NdXi2( 3,   6 ) = 0.0 ;
            ad2NdXi2( 4,   6 ) = -0.75 ;
            ad2NdXi2( 5,   6 ) = -eta ;

            ad2NdXi2( 0,   7 ) = - 1.0 - eta ;
            ad2NdXi2( 1,   7 ) = 0 ;
            ad2NdXi2( 2,   7 ) = 1.0 ;
            ad2NdXi2( 3,   7 ) = -0.75 ;
            ad2NdXi2( 4,   7 ) = 0.0 ;
            ad2NdXi2( 5,   7 ) = - xi ;

            ad2NdXi2( 0,   8 ) = 0.0 ;
            ad2NdXi2( 1,   8 ) = xi - 1.0 ;
            ad2NdXi2( 2,   8 ) = 1.0 ;
            ad2NdXi2( 3,   8 ) = 0.0 ;
            ad2NdXi2( 4,   8 ) = 0.75 ;
            ad2NdXi2( 5,   8 ) = eta ;

            ad2NdXi2( 0,   9 ) =  0.0 ;
            ad2NdXi2( 1,   9 ) =  0.0 ;
            ad2NdXi2( 2,   9 ) = -2.0 ;
            ad2NdXi2( 3,   9 ) = xi2 - 1.0 ;
            ad2NdXi2( 4,   9 ) = eta2 - 1.0 ;
            ad2NdXi2( 5,   9 ) = zeta2 ;

            ad2NdXi2( 0,  10 ) =  0.0 ;
            ad2NdXi2( 1,  10 ) =  0.0 ;
            ad2NdXi2( 2,  10 ) = -2.0 ;
            ad2NdXi2( 3,  10 ) = -xi2 - 1.0 ;
            ad2NdXi2( 4,  10 ) = 1.0 - eta2 ;
            ad2NdXi2( 5,  10 ) = - zeta2 ;

            ad2NdXi2( 0,  11 ) =  0.0 ;
            ad2NdXi2( 1,  11 ) =  0.0 ;
            ad2NdXi2( 2,  11 ) = -2.0 ;
            ad2NdXi2( 3,  11 ) = xi2 + 1.0 ;
            ad2NdXi2( 4,  11 ) = eta2 + 1.0 ;
            ad2NdXi2( 5,  11 ) = zeta2 ;

            ad2NdXi2( 0,  12 ) =  0.0 ;
            ad2NdXi2( 1,  12 ) =  0.0 ;
            ad2NdXi2( 2,  12 ) = -2.0 ;
            ad2NdXi2( 3,  12 ) = 1.0 - xi2 ;
            ad2NdXi2( 4,  12 ) = -1.0 - eta2 ;
            ad2NdXi2( 5,  12 ) = - zeta2 ;

        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_PYRA13_HPP
