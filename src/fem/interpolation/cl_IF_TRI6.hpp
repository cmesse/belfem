//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TRI6_HPP
#define BELFEM_CL_IF_TRI6_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >
        ::element_type() const
        {
            return ElementType::TRI6;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 6 );

            aXiHat( 0, 0 ) = 1.0;
            aXiHat( 1, 0 ) = 0.0;

            aXiHat( 0, 1 ) = 0.0;
            aXiHat( 1, 1 ) = 1.0;

            aXiHat( 0, 2 ) = 0.0;
            aXiHat( 1, 2 ) = 0.0;

            aXiHat( 0, 3 ) = 0.5;
            aXiHat( 1, 3 ) = 0.5;

            aXiHat( 0, 4 ) = 0.0;
            aXiHat( 1, 4 ) = 0.5;

            aXiHat( 0, 5 ) = 0.5;
            aXiHat( 1, 5 ) = 0.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            aN.set_size( 1, 6 );

            aN( 0, 0 ) = xi * ( 2.0 * xi - 1.0 );
            aN( 0, 1 ) = eta * (2.0 * eta - 1.0 );
            aN( 0, 2 ) = ( xi + eta -1.0 ) * ( 2.0 * ( xi + eta ) - 1.0 );
            aN( 0, 3 ) = 4.0 * xi * eta;
            aN( 0, 4 ) = 4.0 * eta * ( 1.0 - xi - eta );
            aN( 0, 5 ) = 4.0 * xi * ( 1.0 - xi - eta );
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            adNdXi.set_size( 2, 6 );

            adNdXi( 0, 0 ) = xi * 4.0 - 1.0;
            adNdXi( 1, 0 ) = 0.0;

            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 1, 1 ) = eta * 4.0 - 1.0;

            adNdXi( 0, 2 ) = 4.0 * ( xi + eta ) - 3.0;
            adNdXi( 1, 2 ) = 4.0 * ( xi + eta ) - 3.0;

            adNdXi( 0, 3 ) = 4.0 * eta;
            adNdXi( 1, 3 ) = 4.0 * xi;

            adNdXi( 0, 4 ) =  - 4.0 * eta;
            adNdXi( 1, 4 ) = 4.0  - 4.0 * xi - 8.0 * eta;

            adNdXi( 0, 5 ) = 4.0  - 8.0 * xi - 4.0 * eta;
            adNdXi( 1, 5 ) = - 4.0 * xi;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 3, 6 );

            ad2NdXi2( 0, 0 ) =  4.0;
            ad2NdXi2( 1, 0 ) =  0.0;
            ad2NdXi2( 2, 0 ) =  0.0;

            ad2NdXi2( 0, 1 ) =  0.0;
            ad2NdXi2( 1, 1 ) =  4.0;
            ad2NdXi2( 2, 1 ) =  0.0;

            ad2NdXi2( 0, 2 ) =  4.0;
            ad2NdXi2( 1, 2 ) =  4.0;
            ad2NdXi2( 2, 2 ) =  4.0;

            ad2NdXi2( 0, 3 ) =  0.0;
            ad2NdXi2( 1, 3 ) =  0.0;
            ad2NdXi2( 2, 3 ) =  4.0;

            ad2NdXi2( 0, 4 ) =  0.0;
            ad2NdXi2( 1, 4 ) = -8.0;
            ad2NdXi2( 2, 4 ) = -4.0;

            ad2NdXi2( 0, 5 ) = -8.0;
            ad2NdXi2( 1, 5 ) =  0.0;
            ad2NdXi2( 2, 5 ) = -4.0;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TRI6_HPP
