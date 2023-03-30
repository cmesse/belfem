//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_QUAD4_HPP
#define BELFEM_CL_IF_QUAD4_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >
        ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >
        ::element_type() const
        {
            return ElementType::QUAD4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 4 );

            aXiHat( 0, 0 ) = -1.000000;
            aXiHat( 1, 0 ) = -1.000000;

            aXiHat( 0, 1 ) =  1.000000;
            aXiHat( 1, 1 ) = -1.000000;

            aXiHat( 0, 2 ) =  1.000000;
            aXiHat( 1, 2 ) =  1.000000;

            aXiHat( 0, 3 ) = -1.000000;
            aXiHat( 1, 3 ) =  1.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            aN.set_size( 1, 4 );

            std::cout << "quad4 " << xi << " " << eta << std::endl ;
            aN( 0, 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
            aN( 0, 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
            aN( 0, 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
            aN( 0, 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;

            aN.print("N");
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            adNdXi.set_size( 2, 4 );

            adNdXi( 0, 0 ) =  0.25 * ( eta - 1.0 );
            adNdXi( 1, 0 ) =  0.25 * ( xi - 1.0 );

            adNdXi( 0, 1 ) = -0.25 * ( eta - 1.0 );
            adNdXi( 1, 1 ) = -0.25 * ( xi + 1.0 );

            adNdXi( 0, 2 ) =  0.25 * ( eta + 1.0 );
            adNdXi( 1, 2 ) =  0.25 * ( xi + 1.0 );

            adNdXi( 0, 3 ) = -0.25 * ( eta + 1.0 );
            adNdXi( 1, 3 ) = -0.25 * ( xi - 1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 3, 4 );

            ad2NdXi2( 0, 0 ) =  0.0;
            ad2NdXi2( 1, 0 ) =  0.0;
            ad2NdXi2( 2, 0 ) =  0.25;

            ad2NdXi2( 0, 1 ) =  0.0;
            ad2NdXi2( 1, 1 ) =  0.0;
            ad2NdXi2( 2, 1 ) = -0.25;

            ad2NdXi2( 0, 2 ) =  0.0;
            ad2NdXi2( 1, 2 ) =  0.0;
            ad2NdXi2( 2, 2 ) =  0.25;

            ad2NdXi2( 0, 3 ) =  0.0;
            ad2NdXi2( 1, 3 ) =  0.0;
            ad2NdXi2( 2, 3 ) = -0.25;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_QUAD4_HPP
