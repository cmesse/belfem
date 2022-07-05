//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_PENTA6_HPP
#define BELFEM_CL_IF_PENTA6_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >
        ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >
        ::element_type() const
        {
            return ElementType::PENTA6;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 6 );

            aXiHat( 0, 0 ) =  1.0;
            aXiHat( 1, 0 ) =  0.0;
            aXiHat( 2, 0 ) = -1.0;

            aXiHat( 0, 1 ) =  0.0;
            aXiHat( 1, 1 ) =  1.0;
            aXiHat( 2, 1 ) = -1.0;

            aXiHat( 0, 2 ) =  0.0;
            aXiHat( 1, 2 ) =  0.0;
            aXiHat( 2, 2 ) = -1.0;

            aXiHat( 0, 3 ) =  1.0;
            aXiHat( 1, 3 ) =  0.0;
            aXiHat( 2, 3 ) =  1.0;

            aXiHat( 0, 4 ) =  0.0;
            aXiHat( 1, 4 ) =  1.0;
            aXiHat( 2, 4 ) =  1.0;

            aXiHat( 0, 5 ) =  0.0;
            aXiHat( 1, 5 ) =  0.0;
            aXiHat( 2, 5 ) =  1.0;

            aXiHat( 0, 6 ) =  0.5;
            aXiHat( 1, 6 ) =  0.5;
            aXiHat( 2, 6 ) = -1.0;
        }
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            aN.set_size( 1, 6 );
            aN( 0, 0 ) = -(xi*(zeta-1.0))*0.5;
            aN( 0, 1 ) = -(eta*(zeta-1.0))*0.5;
            aN( 0, 2 ) = ((zeta-1.0)*(eta+xi-1.0))*0.5;
            aN( 0, 3 ) = (xi*(zeta + 1.0))*0.5;
            aN( 0, 4 ) = (eta*(zeta + 1.0))*0.5;
            aN( 0, 5 ) = ((zeta+1.0)*(1.0-xi-eta))*0.5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            adNdXi.set_size( 3, 6 );

            adNdXi( 0, 0 ) = 0.5*(1.0-zeta);
            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 2, 0 ) = -0.5*xi;

            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 1, 1 ) = 0.5*(1.0-zeta);
            adNdXi( 2, 1 ) = -0.5*eta;

            adNdXi( 0, 2 ) = 0.5*(zeta-1.0);
            adNdXi( 1, 2 ) = 0.5*(zeta-1.0);
            adNdXi( 2, 2 ) = 0.5*(xi+eta-1.0);

            adNdXi( 0, 3 ) = 0.5*(zeta+1);
            adNdXi( 1, 3 ) = 0.0;
            adNdXi( 2, 3 ) = 0.5*xi;

            adNdXi( 0, 4 ) = 0.0;
            adNdXi( 1, 4 ) = 0.5*(zeta+1.0);
            adNdXi( 2, 4 ) = 0.5*eta;

            adNdXi( 0, 5 ) = -0.5*(zeta+1.0);
            adNdXi( 1, 5 ) = -0.5*(zeta+1.0);
            adNdXi( 2, 5 ) = 0.5*(1.0-xi-eta);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 6, 6 );

            ad2NdXi2( 0, 0 ) =  0.0;
            ad2NdXi2( 1, 0 ) =  0.0;
            ad2NdXi2( 2, 0 ) =  0.0;
            ad2NdXi2( 3, 0 ) =  0.0;
            ad2NdXi2( 4, 0 ) = -0.5;
            ad2NdXi2( 5, 0 ) =  0.0;

            ad2NdXi2( 0, 1 ) =  0.0;
            ad2NdXi2( 1, 1 ) =  0.0;
            ad2NdXi2( 2, 1 ) =  0.0;
            ad2NdXi2( 3, 1 ) = -0.5;
            ad2NdXi2( 4, 1 ) =  0.0;
            ad2NdXi2( 5, 1 ) =  0.0;

            ad2NdXi2( 0, 2 ) =  0.0;
            ad2NdXi2( 1, 2 ) =  0.0;
            ad2NdXi2( 2, 2 ) =  0.0;
            ad2NdXi2( 3, 2 ) =  0.5;
            ad2NdXi2( 4, 2 ) =  0.5;
            ad2NdXi2( 5, 2 ) =  0.0;

            ad2NdXi2( 0, 3 ) =  0.0;
            ad2NdXi2( 1, 3 ) =  0.0;
            ad2NdXi2( 2, 3 ) =  0.0;
            ad2NdXi2( 3, 3 ) =  0.0;
            ad2NdXi2( 4, 3 ) =  0.5;
            ad2NdXi2( 5, 3 ) =  0.0;

            ad2NdXi2( 0, 4 ) =  0.0;
            ad2NdXi2( 1, 4 ) =  0.0;
            ad2NdXi2( 2, 4 ) =  0.0;
            ad2NdXi2( 3, 4 ) =  0.5;
            ad2NdXi2( 4, 4 ) =  0.0;
            ad2NdXi2( 5, 4 ) =  0.0;

            ad2NdXi2( 0, 5 ) =  0.0;
            ad2NdXi2( 1, 5 ) =  0.0;
            ad2NdXi2( 2, 5 ) =  0.0;
            ad2NdXi2( 3, 5 ) = -0.5;
            ad2NdXi2( 4, 5 ) = -0.5;
            ad2NdXi2( 5, 5 ) =  0.0;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_PENTA6_HPP
