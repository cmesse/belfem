//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TRI3_HPP
#define BELFEM_CL_IF_TRI3_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >
        ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >
        ::element_type() const
        {
            return ElementType::TRI3;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 3 );

            aXiHat( 0, 0 ) = 1.0;
            aXiHat( 1, 0 ) = 0.0;

            aXiHat( 0, 1 ) = 0.0;
            aXiHat( 1, 1 ) = 1.0;

            aXiHat( 0, 2 ) = 0.0;
            aXiHat( 1, 2 ) = 0.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            aN.set_size( 1, 3 );
            aN( 0, 0 ) = xi;
            aN( 0, 1 ) = eta;
            aN( 0, 2 ) = 1.0 - xi - eta;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            adNdXi.set_size( 2, 3 );
            adNdXi( 0, 0 ) =   1.0;
            adNdXi( 1, 0 ) =   0.0;

            adNdXi( 0, 1 ) =   0.0;
            adNdXi( 1, 1 ) =   1.0;

            adNdXi( 0, 2 ) =  -1.0;
            adNdXi( 1, 2 ) =  -1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 3, 3, 0.0 );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TRI3_HPP
