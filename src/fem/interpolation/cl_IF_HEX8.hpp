//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_HEX8_HPP
#define BELFEM_CL_IF_HEX8_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >
        ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >
        ::element_type() const
        {
            return ElementType::HEX8;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 8 );

            aXiHat( 0, 0 ) = -1.000000;
            aXiHat( 1, 0 ) = -1.000000;
            aXiHat( 2, 0 ) = -1.000000;

            aXiHat( 0, 1 ) =  1.000000;
            aXiHat( 1, 1 ) = -1.000000;
            aXiHat( 2, 1 ) = -1.000000;

            aXiHat( 0, 2 ) =  1.000000;
            aXiHat( 1, 2 ) =  1.000000;
            aXiHat( 2, 2 ) = -1.000000;

            aXiHat( 0, 3 ) = -1.000000;
            aXiHat( 1, 3 ) =  1.000000;
            aXiHat( 2, 3 ) = -1.000000;

            aXiHat( 0, 4 ) = -1.000000;
            aXiHat( 1, 4 ) = -1.000000;
            aXiHat( 2, 4 ) =  1.000000;

            aXiHat( 0, 5 ) =  1.000000;
            aXiHat( 1, 5 ) = -1.000000;
            aXiHat( 2, 5 ) =  1.000000;

            aXiHat( 0, 6 ) =  1.000000;
            aXiHat( 1, 6 ) =  1.000000;
            aXiHat( 2, 6 ) =  1.000000;

            aXiHat( 0, 7 ) = -1.000000;
            aXiHat( 1, 7 ) =  1.000000;
            aXiHat( 2, 7 ) =  1.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            aN.set_size( 1, 8 );

            aN( 0, 0 ) =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 0, 1 ) =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 0, 2 ) =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 0, 3 ) =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 0, 4 ) =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 0, 5 ) =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 0, 6 ) =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 0, 7 ) =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi  = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            adNdXi.set_size( 3, 8 );

            adNdXi( 0, 0 ) = -(  eta - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 0 ) = -(   xi - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 0 ) = -(  eta - 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 1 ) =  (  eta - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 1 ) =  (   xi + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 1 ) =  (  eta - 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 2 ) = -(  eta + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 2 ) = -(   xi + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 2 ) = -(  eta + 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 3 ) =  (  eta + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 3 ) =  (   xi - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 3 ) =  (  eta + 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 4 ) =  (  eta - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 4 ) =  (   xi - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 4 ) =  (  eta - 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 5 ) = -(  eta - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 5 ) = -(   xi + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 5 ) = -(  eta - 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 6 ) =  (  eta + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 6 ) =  (   xi + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 6 ) =  (  eta + 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 7 ) = -(  eta + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 7 ) = -(   xi - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 7 ) = -(  eta + 1 ) * (   xi - 1 ) * 0.125;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            real   xi = aXi( 0 );
            real  eta = aXi( 1 );
            real zeta = aXi( 2 );

            ad2NdXi2.set_size( 6, 8 );

            ad2NdXi2( 0, 0 ) = 0.0;
            ad2NdXi2( 1, 0 ) = 0.0;
            ad2NdXi2( 2, 0 ) = 0.0;
            ad2NdXi2( 3, 0 ) = 0.125 * ( - xi + 1.0 );
            ad2NdXi2( 4, 0 ) = 0.125 * ( - eta + 1.0 );
            ad2NdXi2( 5, 0 ) = 0.125 * ( - zeta + 1.0 );

            ad2NdXi2( 0, 1 ) = 0.0;
            ad2NdXi2( 1, 1 ) = 0.0;
            ad2NdXi2( 2, 1 ) = 0.0;
            ad2NdXi2( 3, 1 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 1 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 1 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 0, 2 ) = 0.0;
            ad2NdXi2( 1, 2 ) = 0.0;
            ad2NdXi2( 2, 2 ) = 0.0;
            ad2NdXi2( 3, 2 ) = 0.125 * ( - xi - 1.0 );
            ad2NdXi2( 4, 2 ) = 0.125 * ( - eta - 1.0 );
            ad2NdXi2( 5, 2 ) = 0.125 * ( - zeta + 1.0 );

            ad2NdXi2( 0, 3 ) = 0.0;
            ad2NdXi2( 1, 3 ) = 0.0;
            ad2NdXi2( 2, 3 ) = 0.0;
            ad2NdXi2( 3, 3 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 3 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 3 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 0, 4 ) = 0.0;
            ad2NdXi2( 1, 4 ) = 0.0;
            ad2NdXi2( 2, 4 ) = 0.0;
            ad2NdXi2( 3, 4 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 4 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 4 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 0, 5 ) = 0.0;
            ad2NdXi2( 1, 5 ) = 0.0;
            ad2NdXi2( 2, 5 ) = 0.0;
            ad2NdXi2( 3, 5 ) = 0.125 * ( - xi - 1.0 );
            ad2NdXi2( 4, 5 ) = 0.125 * ( - eta + 1.0 );
            ad2NdXi2( 5, 5 ) = 0.125 * ( - zeta - 1.0 );

            ad2NdXi2( 0, 6 ) = 0.0;
            ad2NdXi2( 1, 6 ) = 0.0;
            ad2NdXi2( 2, 6 ) = 0.0;
            ad2NdXi2( 3, 6 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 6 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 6 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 0, 7 ) = 0.0;
            ad2NdXi2( 1, 7 ) = 0.0;
            ad2NdXi2( 2, 7 ) = 0.0;
            ad2NdXi2( 3, 7 ) = 0.125 * ( - xi + 1.0 );
            ad2NdXi2( 4, 7 ) = 0.125 * ( - eta - 1.0 );
            ad2NdXi2( 5, 7 ) = 0.125 * ( - zeta - 1.0 );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_HEX8_HPP
