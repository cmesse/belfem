//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_LINE2_HPP
#define BELFEM_CL_IF_LINE2_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >
                ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
            GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >
        ::element_type() const
        {
            return ElementType::LINE2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 1, 2 );

            aXiHat( 0, 0 ) = -1.0;
            aXiHat( 0, 1 ) =  1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >::N(
                        const Vector< real > & aXi,
                              Matrix< real > & aN  ) const
        {
            const real xi = aXi( 0 );

            aN.set_size( 1, 2 );
            aN( 0, 0 ) =  0.5 * ( 1.0 - xi );
            aN( 0, 1 ) =  0.5 * ( 1.0 + xi );
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            adNdXi.set_size( 1, 2 );
            adNdXi( 0, 0 ) =  -0.5;
            adNdXi( 0, 1 ) =   0.5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 1, 2, 0.0 );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_LINE2_HPP
