//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TET4_HPP
#define BELFEM_CL_IF_TET4_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >
        ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >
        ::element_type() const
        {
            return ElementType::TET4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 4 );

            aXiHat( 0, 0 ) = 1.0;
            aXiHat( 1, 0 ) = 0.0;
            aXiHat( 2, 0 ) = 0.0;

            aXiHat( 0, 1 ) = 0.0;
            aXiHat( 1, 1 ) = 1.0;
            aXiHat( 2, 1 ) = 0.0;

            aXiHat( 0, 2 ) = 0.0;
            aXiHat( 1, 2 ) = 0.0;
            aXiHat( 2, 2 ) = 1.0;

            aXiHat( 0, 3 ) = 0.0;
            aXiHat( 1, 3 ) = 0.0;
            aXiHat( 2, 3 ) = 0.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            aN.set_size( 1, 4 );

            aN( 0, 0 ) = xi;
            aN( 0, 1 ) = eta;
            aN( 0, 2 ) = zeta;
            aN( 0, 3 ) = 1.0 - xi - eta - zeta;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            adNdXi.set_size( 3, 4 );

            adNdXi( 0, 0 ) = 1.0;
            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 2, 0 ) = 0.0;

            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 1, 1 ) = 1.0;
            adNdXi( 2, 1 ) = 0.0;

            adNdXi( 0, 2 ) = 0.0;
            adNdXi( 1, 2 ) = 0.0;
            adNdXi( 2, 2 ) = 1.0;

            adNdXi( 0, 3 ) = -1.0;
            adNdXi( 1, 3 ) = -1.0;
            adNdXi( 2, 3 ) = -1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 6, 3, 0.0 );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TET4_HPP
