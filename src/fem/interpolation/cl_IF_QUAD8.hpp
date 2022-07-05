//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_QUAD8_HPP
#define BELFEM_CL_IF_QUAD8_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >
        ::interpolation_order() const
        {
            return InterpolationOrder::SERENDIPITY;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >
        ::element_type() const
        {
            return ElementType::QUAD8;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 8 );

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
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real  xi2 = xi*xi;
            const real eta2 = eta*eta;

            aN.set_size( 1, 8 );

            aN( 0, 0 ) = - ( eta - 1.0 ) * ( xi - 1.0 ) * ( eta + xi + 1.0 ) * 0.25;
            aN( 0, 1 ) =  ( eta - 1.0 ) * ( xi + 1.0 ) * ( eta - xi + 1.0 ) * 0.25;
            aN( 0, 2 ) =  ( eta + 1.0 ) * ( xi + 1.0 ) * ( eta + xi - 1.0 ) * 0.25;
            aN( 0, 3 ) =  ( eta + 1.0 ) * ( xi - 1.0 ) * (  - eta + xi + 1.0 ) * 0.25;
            aN( 0, 4 ) =   ( xi2 - 1.0 ) * ( eta - 1.0 ) * 0.5;
            aN( 0, 5 ) = - ( eta2 - 1.0 ) * ( xi + 1.0 ) * 0.5;
            aN( 0, 6 ) = - ( xi2 - 1.0 ) * ( eta + 1.0 ) * 0.5;
            aN( 0, 7 ) =   ( eta2 - 1.0 ) * ( xi - 1.0 ) * 0.5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real  xi2 = xi*xi;
            const real eta2 = eta*eta;

            adNdXi.set_size( 2, 8 );

            adNdXi( 0, 0 ) = -( eta+xi*2.0 )*( eta-1.0 )*0.25;
            adNdXi( 1, 0 ) = -( eta*2.0+xi )*( xi-1.0 )*0.25;

            adNdXi( 0, 1 ) = ( eta-xi*2.0 )*( eta-1.0 )*0.25;
            adNdXi( 1, 1 ) = ( xi+1.0 )*( eta*2.0-xi )*0.25;

            adNdXi( 0, 2 ) = ( eta+xi*2.0 )*( eta+1.0 )*0.25;
            adNdXi( 1, 2 ) = ( eta*2.0+xi )*( xi+1.0 )*0.25;

            adNdXi( 0, 3 ) = -( eta-xi*2.0 )*( eta+1.0 )*0.25;
            adNdXi( 1, 3 ) = -( xi-1.0 )*( eta*2.0-xi )*0.25;

            adNdXi( 0, 4 ) = xi*( eta-1.0 );
            adNdXi( 1, 4 ) = 0.5*xi2-0.5;

            adNdXi( 0, 5 ) = -0.5*eta2+0.5;
            adNdXi( 1, 5 ) = -eta*( xi+1.0 );

            adNdXi( 0, 6 ) = -xi*( eta+1.0 );
            adNdXi( 1, 6 ) = -0.5*xi2+0.5;

            adNdXi( 0, 7 ) = 0.5*eta2-0.5;
            adNdXi( 1, 7 ) = eta*( xi-1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            // unpack xi and eta from input vector
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );



            ad2NdXi2.set_size( 3, 8 );

            ad2NdXi2( 0, 0 ) = 0.5 - eta * 0.5;
            ad2NdXi2( 1, 0 ) = 0.5 - xi * 0.5;
            ad2NdXi2( 2, 0 ) = 0.25 - 0.5 * ( xi + eta );

            ad2NdXi2( 0, 1 ) = 0.5 - eta * 0.5;
            ad2NdXi2( 1, 1 ) = xi * 0.5 + 0.5;
            ad2NdXi2( 2, 1 ) = 0.5 * ( eta - xi ) - 0.25;

            ad2NdXi2( 0, 2 ) = eta * 0.5 + 0.5;
            ad2NdXi2( 1, 2 ) = xi * 0.5 + 0.5;
            ad2NdXi2( 2, 2 ) = 0.5 * ( xi + eta ) + 0.25;

            ad2NdXi2( 0, 3 ) = eta * 0.5 + 0.5;
            ad2NdXi2( 1, 3 ) = 0.5 - xi * 0.5;
            ad2NdXi2( 2, 3 ) = 0.5 * ( xi - eta ) - 0.25;

            ad2NdXi2( 0, 4 ) = eta - 1.0;
            ad2NdXi2( 1, 4 ) = 0.0;
            ad2NdXi2( 2, 4 ) = xi;

            ad2NdXi2( 0, 5 ) = 0.0;
            ad2NdXi2( 1, 5 ) =  - xi - 1.0;
            ad2NdXi2( 2, 5 ) =  - eta;

            ad2NdXi2( 0, 6 ) =  - eta - 1.0;
            ad2NdXi2( 1, 6 ) = 0.0;
            ad2NdXi2( 2, 6 ) =  - xi;

            ad2NdXi2( 0, 7 ) = 0.0;
            ad2NdXi2( 1, 7 ) = xi - 1.0;
            ad2NdXi2( 2, 7 ) = eta;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_QUAD8_HPP
