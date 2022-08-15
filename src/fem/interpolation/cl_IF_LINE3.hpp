//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_LINE3_HPP
#define BELFEM_CL_IF_LINE3_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3>
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3>
        ::element_type() const
        {
            return ElementType::LINE3;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3>
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 1, 3 );

            aXiHat( 0, 0 ) = -1.0;
            aXiHat( 0, 1 ) = 1.0;
            aXiHat( 0, 2 ) = 0.0;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3>::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {
            const real xi = aXi( 0 );
            const real xi2 = xi * xi;

            aN.set_size( 1, 3 );
            aN( 0, 0 ) = 0.5 * ( xi2 - xi );
            aN( 0, 1 ) = 0.5 * ( xi2 + xi );
            aN( 0, 2 ) = 1.0 - xi2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
            const real xi = aXi( 0 );
            adNdXi.set_size( 1, 3 );
            adNdXi( 0, 0 ) = xi - 0.5;
            adNdXi( 0, 1 ) = xi + 0.5;
            adNdXi( 0, 2 ) = -2.0 * xi;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
            ad2NdXi2.set_size( 1, 3, 0.0 );
            ad2NdXi2( 0, 0 ) = 1.0;
            ad2NdXi2( 0, 1 ) = 1.0;
            ad2NdXi2( 0, 2 ) = -2.0;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_LINE3_HPP
