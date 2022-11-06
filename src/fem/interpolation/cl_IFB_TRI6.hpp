//
// Created by Christian Messe on 11/05/22
//

#ifndef BELFEM_CL_IFB_TRI6_HPP
#define BELFEM_CL_IFB_TRI6_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >
        ::element_type() const
        {
            return ElementType::TRI6;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >
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
                GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = 1. - xi - eta ;

            aN.set_size( 1, 6 );

            aN( 0, 0 ) = xi*xi ;
            aN( 0, 1 ) = eta*eta ;
            aN( 0, 2 ) = zeta*zeta ;
            aN( 0, 3 ) = 2. * eta * xi ;
            aN( 0, 4 ) = 2. * eta * zeta ;
            aN( 0, 5 ) = 2. * xi * zeta ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );
            const real zeta = 1. - xi - eta ;

            adNdXi.set_size( 2, 6 );

            adNdXi( 0, 0 ) = xi + xi ;
            adNdXi( 1, 0 ) = 0 ;

            adNdXi( 0, 1 ) = 0 ;
            adNdXi( 1, 1 ) = eta + eta ;

            adNdXi( 0, 2 ) = - zeta - zeta ;
            adNdXi( 1, 2 ) = - zeta - zeta ;

            adNdXi( 0, 3 ) = eta + eta ;
            adNdXi( 1, 3 ) = xi + xi ;

            adNdXi( 0, 4 ) = - eta - eta ;
            adNdXi( 1, 4 ) = 2.-xi-xi-eta-eta-eta-eta;

            adNdXi( 0, 5 ) = 2.-xi-xi-xi-xi-eta-eta;
            adNdXi( 1, 5 ) = -xi - xi ;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 3, 6 );

            ad2NdXi2( 0, 0 ) =  2. ;
            ad2NdXi2( 1, 0 ) =  0. ;
            ad2NdXi2( 2, 0 ) =  0. ;

            ad2NdXi2( 0, 1 ) =  0. ;
            ad2NdXi2( 1, 1 ) =  2. ;
            ad2NdXi2( 2, 1 ) =  0. ;

            ad2NdXi2( 0, 2 ) =  2. ;
            ad2NdXi2( 1, 2 ) =  2. ;
            ad2NdXi2( 2, 2 ) =  2. ;

            ad2NdXi2( 0, 3 ) =  0. ;
            ad2NdXi2( 1, 3 ) =  0. ;
            ad2NdXi2( 2, 3 ) =  2. ;

            ad2NdXi2( 0, 4 ) =  0. ;
            ad2NdXi2( 1, 4 ) = -4. ;
            ad2NdXi2( 2, 4 ) = -2. ;

            ad2NdXi2( 0, 5 ) = -4. ;
            ad2NdXi2( 1, 5 ) = 0. ;
            ad2NdXi2( 2, 5 ) = -2. ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IFB_TRI6_HPP
