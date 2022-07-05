//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_QUAD9_HPP
#define BELFEM_CL_IF_QUAD9_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >
        ::element_type() const
        {
            return ElementType::QUAD9;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 9 );

            aXiHat( 0, 0 ) = -1.000000;
            aXiHat( 1, 0 ) = -1.000000;

            aXiHat( 0, 1 ) =  1.000000;
            aXiHat( 1, 1 ) = -1.000000;

            aXiHat( 0, 2 ) =  1.000000;
            aXiHat( 1, 2 ) =  1.000000;

            aXiHat( 0, 3 ) = -1.000000;
            aXiHat( 1, 3 ) =  1.000000;

            aXiHat( 0, 4 ) =  0.000000;
            aXiHat( 1, 4 ) = -1.000000;

            aXiHat( 0, 5 ) =  1.000000;
            aXiHat( 1, 5 ) =  0.000000;

            aXiHat( 0, 6 ) =  0.000000;
            aXiHat( 1, 6 ) =  1.000000;

            aXiHat( 0, 7 ) = -1.000000;
            aXiHat( 1, 7 ) =  0.000000;

            aXiHat( 0, 8 ) =  0.000000;
            aXiHat( 1, 8 ) =  0.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real    c = xi * eta * 0.25;
            const real  xi2 = xi*xi;
            const real eta2 = eta*eta;

            aN.set_size( 1, 9 );

            aN( 0, 0 ) = ( c * ( eta - 1.0 ) * (xi - 1.0) );
            aN( 0, 1 ) = ( c * ( eta - 1.0 ) * (xi + 1.0) );
            aN( 0, 2 ) = ( c * ( eta + 1.0 ) * (xi + 1.0) );
            aN( 0, 3 ) = ( c * ( eta + 1.0 ) * (xi - 1.0) );
            aN( 0, 4 ) = ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
            aN( 0, 5 ) = ( xi * ( 1.0 - eta2)*( xi + 1.0 ) )*0.5;
            aN( 0, 6 ) = ( eta * (1.0 - xi2)*( eta + 1.0 ) )*0.5;
            aN( 0, 7 ) = ( xi*( 1.0 - eta2 )*( xi - 1.0 ) )*0.5;
            aN( 0, 8 ) = ( eta2 - 1.0 )*( xi2 - 1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            // often used constants
            const real    c = xi*eta;
            const real  xi2 = xi*xi;
            const real eta2 = eta*eta;

            adNdXi.set_size( 2, 9 );

            adNdXi( 0, 0 ) = ( eta * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 0 ) =  ( xi * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 1 ) = ( eta * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 1 ) =  ( xi * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 2 ) = ( eta * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 2 ) =  ( xi * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 3 ) = ( eta * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 3 ) =  ( xi * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 4 ) =  - c * ( eta - 1.0 );
            adNdXi( 1, 4 ) =  -( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;

            adNdXi( 0, 5 ) = -( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.5;
            adNdXi( 1, 5 ) = - c * ( xi + 1.0 );

            adNdXi( 0, 6 ) = - c * ( eta + 1.0 );
            adNdXi( 1, 6 ) =  -( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;

            adNdXi( 0, 7 ) = -( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.5;
            adNdXi( 1, 7 ) = - c * ( xi - 1.0 );

            adNdXi( 0, 8 ) = 2.0 * xi * ( eta2 - 1.0 );
            adNdXi( 1, 8 ) = 2.0 * eta * ( xi2 - 1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            // often used constants
            const real  xi2 =xi*xi;
            const real eta2 = eta*eta;

            ad2NdXi2.set_size( 3, 9 );

            ad2NdXi2( 0, 0 ) = ( eta * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 0 ) = ( xi * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 0 ) = ( ( 2.0 * eta - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 1 ) = ( eta * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 1 ) = ( xi * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 1 ) = ( ( 2.0 * eta - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 2 ) = ( eta * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 2 ) = ( xi * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 2 ) = ( ( 2.0 * eta + 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 3 ) = ( eta * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 3 ) = ( xi * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 3 ) = ( ( 2.0 * eta + 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 4 ) = -eta * ( eta - 1.0 );
            ad2NdXi2( 1, 4 ) = 1.0 - xi2;
            ad2NdXi2( 2, 4 ) = -xi * ( 2.0 * eta - 1.0 );

            ad2NdXi2( 0, 5 ) = 1.0 - eta2;
            ad2NdXi2( 1, 5 ) = -xi * ( xi + 1.0 );
            ad2NdXi2( 2, 5 ) = -eta * ( 2.0 * xi + 1.0 );

            ad2NdXi2( 0, 6 ) = -eta * ( eta + 1.0 );
            ad2NdXi2( 1, 6 ) = 1.0 - xi2;
            ad2NdXi2( 2, 6 ) = -xi * ( 2.0 * eta + 1.0 );

            ad2NdXi2( 0, 7 ) = 1.0 - eta2;
            ad2NdXi2( 1, 7 ) = -xi * ( xi - 1.0 );
            ad2NdXi2( 2, 7 ) = -eta * ( 2.0 * xi - 1.0 );

            ad2NdXi2( 0, 8 ) = 2.0 * eta2 - 2.0;
            ad2NdXi2( 1, 8 ) = 2.0 * xi2 - 2.0;
            ad2NdXi2( 2, 8 ) = 4.0 * eta * xi;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_QUAD9_HPP
