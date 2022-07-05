//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TRUSS4_HPP
#define BELFEM_CL_IF_TRUSS4_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >
        ::element_type() const
        {
            return ElementType::LINE4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 1, 4 );

            aXiHat( 0, 0 ) = -1.000000;
            aXiHat( 0, 1 ) = 1.000000;
            aXiHat( 0, 2 ) = -1.0 / 3.0;
            aXiHat( 0, 3 ) = 1.0 / 3.0;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {

            aN.set_size( 1, 4 );

            const real xi = aXi( 0 );

            aN.set_size( 1, 4, 0.0 );
            aN( 0, 0 ) = ( 3.0 * xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( 1.0 - xi );
            aN( 0, 1 ) = ( 3.0 * xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( 1.0 + xi );
            aN( 0, 2 ) = 9.0 * ( xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( xi - 1.0 );
            aN( 0, 3 ) = 9.0 * ( xi + 1.0 ) * ( 3.0 * xi + 1.0 ) * ( 1.0 - xi );

            aN /= 16.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
            const real xi = aXi( 0 );
            const real xi2 = xi * xi;

            // set adNdXi
            adNdXi.set_size( 1, 4 );
            adNdXi( 0, 0 ) = -27.0 * xi2 + 18.0 * xi + 1.0;
            adNdXi( 0, 1 ) = 27.0 * xi2 + 18.0 * xi - 1.0;
            adNdXi( 0, 2 ) = 9.0 * ( 9.0 * xi2 - 2.0 * xi - 3.0 );
            adNdXi( 0, 3 ) = 9.0 * ( -9.0 * xi2 - 2.0 * xi + 3.0 );

            adNdXi /= 16.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
            const real xi = aXi( 0 );

            ad2NdXi2.set_size( 1, 4 );

            ad2NdXi2( 0, 0 ) = 1.0 - 3.0 * xi;
            ad2NdXi2( 0, 1 ) = 1.0 + 3.0 * xi;
            ad2NdXi2( 0, 2 ) = 9.0 * xi - 1.0;
            ad2NdXi2( 0, 3 ) = -9.0 * xi - 1.0;

            ad2NdXi2 *= 9.0 / 8.0;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_TRUSS4_HPP
