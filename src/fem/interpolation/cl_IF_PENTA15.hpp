//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_PENTA15_HPP
#define BELFEM_CL_IF_PENTA15_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >
        ::interpolation_order() const
        {
            return InterpolationOrder::SERENDIPITY;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >
        ::element_type() const
        {
            return ElementType::PENTA15;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 15 );

            aXiHat( 0,  0 ) =  1.0;
            aXiHat( 1,  0 ) =  0.0;
            aXiHat( 2,  0 ) = -1.0;

            aXiHat( 0,  1 ) =  0.0;
            aXiHat( 1,  1 ) =  1.0;
            aXiHat( 2,  1 ) = -1.0;

            aXiHat( 0,  2 ) =  0.0;
            aXiHat( 1,  2 ) =  0.0;
            aXiHat( 2,  2 ) = -1.0;

            aXiHat( 0,  3 ) =  1.0;
            aXiHat( 1,  3 ) =  0.0;
            aXiHat( 2,  3 ) =  1.0;

            aXiHat( 0,  4 ) =  0.0;
            aXiHat( 1,  4 ) =  1.0;
            aXiHat( 2,  4 ) =  1.0;

            aXiHat( 0,  5 ) =  0.0;
            aXiHat( 1,  5 ) =  0.0;
            aXiHat( 2,  5 ) =  1.0;

            aXiHat( 0,  6 ) =  0.5;
            aXiHat( 1,  6 ) =  0.5;
            aXiHat( 2,  6 ) = -1.0;

            aXiHat( 0,  7 ) =  0.0;
            aXiHat( 1,  7 ) =  0.5;
            aXiHat( 2,  7 ) = -1.0;

            aXiHat( 0,  8 ) =  0.5;
            aXiHat( 1,  8 ) =  0.0;
            aXiHat( 2,  8 ) = -1.0;

            aXiHat( 0,  9 ) =  1.0;
            aXiHat( 1,  9 ) =  0.0;
            aXiHat( 2,  9 ) =  0.0;

            aXiHat( 0, 10 ) =  0.0;
            aXiHat( 1, 10 ) =  1.0;
            aXiHat( 2, 10 ) =  0.0;

            aXiHat( 0, 11 ) =  0.0;
            aXiHat( 1, 11 ) =  0.0;
            aXiHat( 2, 11 ) =  0.0;

            aXiHat( 0, 12 ) =  0.5;
            aXiHat( 1, 12 ) =  0.5;
            aXiHat( 2, 12 ) =  1.0;

            aXiHat( 0, 13 ) =  0.0;
            aXiHat( 1, 13 ) =  0.5;
            aXiHat( 2, 13 ) =  1.0;

            aXiHat( 0, 14 ) =  0.5;
            aXiHat( 1, 14 ) =  0.0;
            aXiHat( 2, 14 ) =  1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real zeta2 = zeta*zeta;

            aN.set_size( 1, 15 );

            aN( 0,  0 ) = (xi*(zeta-1.0)*(xi*(-2.0)+zeta+2.0))*0.5;
            aN( 0,  1 ) = (eta*(zeta-1.0)*(eta*(-2.0)+zeta+2.0))*0.5;
            aN( 0,  2 ) = (1.0-zeta)*(eta+xi-1.0)*(eta*2.0+xi*2.0+zeta)*0.5;
            aN( 0,  3 ) = (xi*(zeta+1.0)*(xi*2.0+zeta-2.0))*0.5;
            aN( 0,  4 ) = (eta*(zeta+1.0)*(eta*2.0+zeta-2.0))*0.5;
            aN( 0,  5 ) = ((zeta+1.0)*(eta+xi-1.0)*(eta*2.0+xi*2.0-zeta))*0.5;
            aN( 0,  6 ) = eta*xi*(1.0-zeta)*2.0;
            aN( 0,  7 ) = eta*(zeta-1.0)*(eta+xi-1.0)*2.0;
            aN( 0,  8 ) = xi*(zeta-1.0)*(eta+xi-1.0)*2.0;
            aN( 0,  9 ) = xi*(1.0-zeta2);
            aN( 0, 10 ) = eta*(1.0-zeta2);
            aN( 0, 11 ) = (zeta2-1.0)*(eta+xi-1.0);
            aN( 0, 12 ) = eta*xi*(zeta+1.0)*2.0;
            aN( 0, 13 ) = eta*(zeta+1.0)*(1.0-eta-xi)*2.0;
            aN( 0, 14 ) = xi*(zeta+1.0)*(1.0-eta-xi)*2.0;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            adNdXi.set_size( 3, 15 );

            adNdXi( 0,  5 ) = 0.5*((zeta - 1.0)*(zeta - 4.0*xi + 2.0));
            adNdXi( 1,  0 ) = 0.0;
            adNdXi( 2,  0 ) = 0.5*xi*(1.0+2.0*(zeta-xi));

            adNdXi( 0,  1 ) = 0.0;
            adNdXi( 1,  1 ) = 0.5*(zeta-1.0)*(zeta-4.0*eta+2.0);
            adNdXi( 2,  1 ) = 0.5*eta*(1.0+2.0*(zeta-eta));

            adNdXi( 0,  2 ) = 0.5*(1.0-zeta)*(zeta+4.0*(xi+eta)-2.0);
            adNdXi( 1,  2 ) = 0.5*(1.0-zeta)*(zeta+4.0*(xi+eta)-2.0);
            adNdXi( 2,  2 ) = 0.5*(1-xi-eta)*(2.0*(xi+eta+zeta)-1.0);

            adNdXi( 0,  3 ) = 0.5*(1.0+zeta)*(4.0*xi+zeta-2.0);
            adNdXi( 1,  3 ) = 0.0;
            adNdXi( 2,  3 ) = xi*(xi+zeta-0.5);

            adNdXi( 0,  4 ) = 0.0;
            adNdXi( 1,  4 ) = 0.5*(1.0+zeta)*(4.0*eta+zeta-2.0);
            adNdXi( 2,  4 ) = eta*(eta+zeta-0.5);

            adNdXi( 0,  5 ) = 0.5*(1.0+zeta)*(4.0*(xi+eta)-zeta-2.0);
            adNdXi( 1,  5 ) = 0.5*(1.0+zeta)*(4.0*(xi+eta)-zeta-2.0);
            adNdXi( 2,  5 ) = 0.5*(xi+eta-1.0)*(2.0*(xi+eta-zeta)-1.0);

            adNdXi( 0,  6 ) = 2.0*eta*(1.0-zeta);
            adNdXi( 1,  6 ) = 2.0*xi*(1.0-zeta);
            adNdXi( 2,  6 ) = -2.0*xi*eta;

            adNdXi( 0,  7 ) = 2.0*eta*(zeta-1.0);
            adNdXi( 1,  7 ) = 2.0*(zeta-1.0)*(xi+2.0*eta-1.0);
            adNdXi( 2,  7 ) = 2.0*eta*(xi+eta-1.0);

            adNdXi( 0,  8 ) = 2.0*(zeta-1.0)*(2.0*xi+eta-1.0);
            adNdXi( 1,  8 ) = 2.0*xi*(zeta-1.0);
            adNdXi( 2,  8 ) = 2.0*xi*(xi+eta-1.0);

            adNdXi( 0,  9 ) = 1.0-zeta*zeta;
            adNdXi( 1,  9 ) = 0.0;
            adNdXi( 2,  9 ) = -2.0*xi*zeta;

            adNdXi( 0, 10 ) = 0.0;
            adNdXi( 1, 10 ) = 1.0-zeta*zeta;
            adNdXi( 2, 10 ) = -2.0*eta*zeta;

            adNdXi( 0, 11 ) = zeta*zeta-1.0;
            adNdXi( 1, 11 ) = zeta*zeta-1.0;
            adNdXi( 2, 11 ) = 2.0*zeta*(xi+eta-1.0);

            adNdXi( 0, 12 ) = 2.0*eta*(1.0+zeta);
            adNdXi( 1, 12 ) = 2.0*xi*(1.0+zeta);
            adNdXi( 2, 12 ) = 2.0*xi*eta;

            adNdXi( 0, 13 ) = -2.0*eta*(1.0+zeta);
            adNdXi( 1, 13 ) = 2.0*(1.0-xi-2.0*eta)*(1.0+zeta);
            adNdXi( 2, 13 ) = 2.0*eta*(1.0-xi-eta);

            adNdXi( 0, 14 ) = 2.0*(1.0-2.0*xi-eta)*(1.0+zeta);
            adNdXi( 1, 14 ) = -2.0*xi*(1.0+zeta);
            adNdXi( 2, 14 ) = 2.0*xi*(1.0-xi-eta);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            ad2NdXi2.set_size( 6, 15 );

            ad2NdXi2( 0,  0 ) = 2.0-2.0*zeta;
            ad2NdXi2( 1,  0 ) = 0.0;
            ad2NdXi2( 2,  0 ) = xi;
            ad2NdXi2( 3,  0 ) = 0.0;
            ad2NdXi2( 4,  0 ) = zeta-2.0*xi+0.5;
            ad2NdXi2( 5,  0 ) = 0.0;

            ad2NdXi2( 0,  1 ) = 0.0;
            ad2NdXi2( 1,  1 ) = 2.0-2.0*zeta;
            ad2NdXi2( 2,  1 ) = eta;
            ad2NdXi2( 3,  1 ) = zeta-2.0*eta+0.5;
            ad2NdXi2( 4,  1 ) = 0.0;
            ad2NdXi2( 5,  1 ) = 0.0;

            ad2NdXi2( 0,  2 ) = 2.0-2.0*zeta;
            ad2NdXi2( 1,  2 ) = 2.0-2.0*zeta;
            ad2NdXi2( 2,  2 ) = 1.0-xi-eta;
            ad2NdXi2( 3,  2 ) = 1.5-2.0*(xi+eta)-zeta;
            ad2NdXi2( 4,  2 ) = 1.5-2.0*(xi+eta)-zeta;
            ad2NdXi2( 5,  2 ) = 2.0-2.0*zeta;

            ad2NdXi2( 0,  3 ) = 2.0*zeta+2.0;
            ad2NdXi2( 1,  3 ) = 0.0;
            ad2NdXi2( 2,  3 ) = xi;
            ad2NdXi2( 3,  3 ) = 0.0;
            ad2NdXi2( 4,  3 ) = 2.0*xi+zeta-0.5;
            ad2NdXi2( 5,  3 ) = 0.0;

            ad2NdXi2( 0,  4 ) = 0.0;
            ad2NdXi2( 1,  4 ) = 2.0*zeta+2.0;
            ad2NdXi2( 2,  4 ) = eta;
            ad2NdXi2( 3,  4 ) = 2.0*eta+zeta-0.5;
            ad2NdXi2( 4,  4 ) = 0.0;
            ad2NdXi2( 5,  4 ) = 0.0;

            ad2NdXi2( 0,  5 ) = 2.0*zeta+2.0;
            ad2NdXi2( 1,  5 ) = 2.0*zeta+2.0;
            ad2NdXi2( 2,  5 ) = 1.0-eta-xi;
            ad2NdXi2( 3,  5 ) = 2.0*eta+2.0*xi-zeta-1.5;
            ad2NdXi2( 4,  5 ) = 2.0*eta+2.0*xi-zeta-1.5;
            ad2NdXi2( 5,  5 ) = 2.0*zeta+2.0;

            ad2NdXi2( 0,  6 ) = 0.0;
            ad2NdXi2( 1,  6 ) = 0.0;
            ad2NdXi2( 2,  6 ) = 0.0;
            ad2NdXi2( 3,  6 ) = -2.0*xi;
            ad2NdXi2( 4,  6 ) = -2.0*eta;
            ad2NdXi2( 5,  6 ) =  2.0-2.0*zeta;

            ad2NdXi2( 0,  7 ) = 0.0;
            ad2NdXi2( 1,  7 ) = zeta*4.0-4.0;
            ad2NdXi2( 2,  7 ) = 0.0;
            ad2NdXi2( 3,  7 ) = eta*4.0+2.0*xi-2.0;
            ad2NdXi2( 4,  7 ) = 2.0*eta;
            ad2NdXi2( 5,  7 ) = 2.0*zeta-2.0;

            ad2NdXi2( 0,  8 ) = zeta*4.0-4.0;
            ad2NdXi2( 1,  8 ) = 0.0;
            ad2NdXi2( 2,  8 ) = 0.0;
            ad2NdXi2( 3,  8 ) = 2.0*xi;
            ad2NdXi2( 4,  8 ) = 4.0*xi+2.0*eta-2.0;
            ad2NdXi2( 5,  8 ) = 2.0*zeta-2.0;

            ad2NdXi2( 0,  9 ) = 0.0;
            ad2NdXi2( 1,  9 ) = 0.0;
            ad2NdXi2( 2,  9 ) = -2.0*xi;
            ad2NdXi2( 3,  9 ) = 0.0;
            ad2NdXi2( 4,  9 ) = -2.0*zeta;
            ad2NdXi2( 5,  9 ) = 0.0;

            ad2NdXi2( 0, 10 ) = 0.0;
            ad2NdXi2( 1, 10 ) = 0.0;
            ad2NdXi2( 2, 10 ) = -2.0*eta;
            ad2NdXi2( 3, 10 ) = -2.0*zeta;
            ad2NdXi2( 4, 10 ) = 0.0;
            ad2NdXi2( 5, 10 ) = 0.0;

            ad2NdXi2( 0, 11 ) = 0.0;
            ad2NdXi2( 1, 11 ) = 0.0;
            ad2NdXi2( 2, 11 ) = 2.0*xi+2.0*eta-2.0;
            ad2NdXi2( 3, 11 ) = 2.0*zeta;
            ad2NdXi2( 4, 11 ) = 2.0*zeta;
            ad2NdXi2( 5, 11 ) = 0.0;

            ad2NdXi2( 0, 12 ) = 0.0;
            ad2NdXi2( 1, 12 ) = 0.0;
            ad2NdXi2( 2, 12 ) = 0.0;
            ad2NdXi2( 3, 12 ) = 2.0*xi;
            ad2NdXi2( 4, 12 ) = 2.0*eta;
            ad2NdXi2( 5, 12 ) = 2.0*zeta+2.0;

            ad2NdXi2( 0, 13 ) = 0.0;
            ad2NdXi2( 1, 13 ) = -4.0*zeta-4.0;
            ad2NdXi2( 2, 13 ) = 0.0;
            ad2NdXi2( 3, 13 ) = 2.0-2.0*xi-4.0*eta;
            ad2NdXi2( 4, 13 ) = -2.0*eta;
            ad2NdXi2( 5, 13 ) = -2.0*zeta-2.0;

            ad2NdXi2( 0, 14 ) = -4.0*zeta-4.0;
            ad2NdXi2( 1, 14 ) = 0.0;
            ad2NdXi2( 2, 14 ) = 0.0;
            ad2NdXi2( 3, 14 ) =  -2.0*xi;
            ad2NdXi2( 4, 14 ) =  2.0-4.0*xi-2.0*eta;
            ad2NdXi2( 5, 14 ) = -2.0*zeta-2.0;

        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_PENTA15_HPP
