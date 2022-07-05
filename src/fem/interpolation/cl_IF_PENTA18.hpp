//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_PENTA18_HPP
#define BELFEM_CL_IF_PENTA18_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >
        ::element_type() const
        {
            return ElementType::PENTA18;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 18 );

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

            aXiHat( 0, 15 ) =  0.5;
            aXiHat( 1, 15 ) =  0.5;
            aXiHat( 2, 15 ) =  0.0;

            aXiHat( 0, 16 ) =  0.0;
            aXiHat( 1, 16 ) =  0.5;
            aXiHat( 2, 16 ) =  0.0;

            aXiHat( 0, 17 ) =  0.5;
            aXiHat( 1, 17 ) =  0.0;
            aXiHat( 2, 17 ) =  0.0;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real zeta2 = zeta*zeta;

            aN.set_size( 1, 18 );

            aN( 0,  0 ) = 0.5*(xi*zeta*(xi*2.0-1.0)*(zeta-1.0));
            aN( 0,  1 ) = 0.5*(eta*zeta*(eta*2.0-1.0)*(zeta-1.0));
            aN( 0,  2 ) = 0.5*((1.0-xi-eta)*(2*(xi+eta)-1.0)*(1.0-zeta)*zeta);
            aN( 0,  3 ) = 0.5*(xi*zeta*(xi*2.0-1.0)*(zeta+1.0));
            aN( 0,  4 ) = 0.5*(eta*zeta*(eta*2.0-1.0)*(zeta+1.0));
            aN( 0,  5 ) = 0.5*((1.0-xi-eta)*(1.0-2.0*(xi+eta))*zeta*(1.0+zeta));
            aN( 0,  6 ) = xi*eta*zeta*(zeta-1.0)*2.0;
            aN( 0,  7 ) = 2.0*eta*zeta*(1.0-zeta)*(eta+xi-1.0);
            aN( 0,  8 ) = 2.0*xi*zeta*(1.0-zeta)*(eta+xi-1.0);
            aN( 0,  9 ) = xi*(1.0-2.0*xi)*(zeta2-1.0);
            aN( 0, 10 ) = eta*(1.0-2.0*eta)*(zeta2-1.0);
            aN( 0, 11 ) = (1-xi-eta)*(2.0*(xi+eta)-1)*(zeta2-1.0);
            aN( 0, 12 ) = xi*eta*zeta*(1.0+zeta)*2.0;
            aN( 0, 13 ) = 2.0*eta*zeta*(1.0+zeta)*(1.0-xi-eta);
            aN( 0, 14 ) = 2.0*xi*zeta*(1.0+zeta)*(1.0-xi-eta);
            aN( 0, 15 ) = 4.0*eta*xi*(1.0-zeta2);
            aN( 0, 16 ) = 4.0*eta*(1.0-zeta2)*(1.0-xi-eta);
            aN( 0, 17 ) = 4.0*xi*(1.0-zeta2)*(1.0-xi-eta);

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real zeta2 = aXi( 2 );

            adNdXi.set_size( 3, 18 );

            adNdXi( 0,  0 ) = zeta*(0.5-0.5*zeta+xi*2.0*(zeta-1.0));
            adNdXi( 1,  0 ) = 0.0;
            adNdXi( 2,  0 ) = xi*(0.5-zeta+xi*(2.0*zeta-1.0));

            adNdXi( 0,  1 ) = 0.0;
            adNdXi( 1,  1 ) = zeta*(0.5-0.5*zeta+2.0*eta*(zeta-1.0));
            adNdXi( 2,  1 ) = eta*(0.5-zeta+eta*(2.0*zeta-1.0));

            adNdXi( 0,  2 ) = zeta*(1.5-1.5*zeta+2.0*(xi+eta)*(zeta-1.0));
            adNdXi( 1,  2 ) = zeta*(1.5-1.5*zeta+2.0*(xi+eta)*(zeta-1.0));
            adNdXi( 2,  2 ) = 0.5*(1.0-xi-eta)*(1.0-2.0*(xi+eta))*(2.0*zeta-1.0);

            adNdXi( 0,  3 ) = zeta*(xi*(2.0+2.0*zeta)-0.5-0.5*zeta);
            adNdXi( 1,  3 ) = 0.0;
            adNdXi( 2,  3 ) = xi*(xi*(1.+2.0*zeta)-0.5-zeta);

            adNdXi( 0,  4 ) = 0.0;
            adNdXi( 1,  4 ) = zeta*(eta*(2.0+2.0*zeta)-0.5-0.5*zeta);
            adNdXi( 2,  4 ) = eta*(eta*(1.+2.0*zeta)-0.5-zeta);

            adNdXi( 0,  5 ) = zeta*(2.0*(xi+eta)*(1.0+zeta)-1.5-1.5*zeta);
            adNdXi( 1,  5 ) = 0.5*(zeta*(zeta + 1.0)*(4.0*(xi+eta)-3.0));
            adNdXi( 2,  5 ) = 2.0*(1.0-xi-eta)*(0.5-xi-eta)*(0.5+zeta);

            adNdXi( 0,  6 ) = 2.0*eta*zeta*(zeta-1.0);
            adNdXi( 1,  6 ) = 2.0*xi*zeta*(zeta-1.0);
            adNdXi( 2,  6 ) = 2.0*eta*xi*(2.0*zeta-1.0);

            adNdXi( 0,  7 ) = 2.0*eta*zeta*(1.0-zeta);
            adNdXi( 1,  7 ) = 4.0*zeta*(eta+0.5*xi-0.5)*(1.0-zeta);
            adNdXi( 2,  7 ) = 4.0*eta*(1.0-xi-eta)*(zeta-0.5);

            adNdXi( 0,  8 ) = 2.0*(1.0-2.0*xi-eta)*(zeta-1.0)*zeta;
            adNdXi( 1,  8 ) = 2.0*xi*zeta*(1.0-zeta);
            adNdXi( 2,  8 ) = 2.0*xi*(1.0-xi-eta)*(2.0*zeta-1.0);

            adNdXi( 0,  9 ) = xi*(4.-4.0*zeta2)+zeta2-1.0;
            adNdXi( 1,  9 ) = 0.0;
            adNdXi( 2,  9 ) = 2.0*xi*zeta*(1.0-2.0*xi);

            adNdXi( 0, 10 ) = 0.0;
            adNdXi( 1, 10 ) = zeta2-1.0+eta*(4.-4.0*zeta2);
            adNdXi( 2, 10 ) = 2.0*eta*zeta*(1.0-2.0*eta);

            adNdXi( 0, 11 ) = (1.0-zeta2)*(4.0*(xi+eta)-3.0);
            adNdXi( 1, 11 ) = (1.0-zeta2)*(4.0*(xi+eta)-3.0);
            adNdXi( 2, 11 ) = 4.0*(1.0-xi-eta)*(xi+eta-0.5)*zeta;

            adNdXi( 0, 12 ) = 2.0*eta*zeta*(1.0+zeta);
            adNdXi( 1, 12 ) = 2.0*xi*zeta*(1.0+zeta);
            adNdXi( 2, 12 ) = 2.0*xi*eta*(1.0+2.0*zeta);

            adNdXi( 0, 13 ) = -2.0*eta*zeta*(1.0+zeta);
            adNdXi( 1, 13 ) =  2.0*zeta*(1.0-2.0*eta-xi)*(1.0+zeta);
            adNdXi( 2, 13 ) =  2.0*eta*(1.0+2.0*zeta-(xi+eta)*(1.0+2.0*zeta));

            adNdXi( 0, 14 ) =  2.0*(1.0-eta-2.0*xi)*zeta*(1.0+zeta);
            adNdXi( 1, 14 ) = -2.0*xi*zeta*(1.0+zeta);
            adNdXi( 2, 14 ) =  2.0*xi*(1.0-(xi+eta)*(1.0+2.0*zeta)+2.0*zeta);

            adNdXi( 0, 15 ) =  4.0*eta*(1.0-zeta2);
            adNdXi( 1, 15 ) =  4.0*xi*(1.0-zeta2);
            adNdXi( 2, 15 ) = -8.0*xi*eta*zeta;

            adNdXi( 0, 16 ) = 4.0*eta*(zeta2-1.0);
            adNdXi( 1, 16 ) = 4.0*(1.0-zeta2+(xi+2.0*eta)*(zeta2-1.0));
            adNdXi( 2, 16 ) = 8.0*eta*zeta*(eta+xi-1.0);

            adNdXi( 0, 17 ) = 4.0*(1.0-2.0*xi-eta)*(1.0-zeta2);
            adNdXi( 1, 17 ) = 4.0*xi*(zeta2-1.0);
            adNdXi( 2, 17 ) = 8.0*xi*zeta*(xi+eta-1.0);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real zeta2 = aXi( 2 );

            ad2NdXi2.set_size( 6, 18 );

            ad2NdXi2( 0,  0 ) = 2.0*(zeta-1.0)*zeta;
            ad2NdXi2( 1,  0 ) = 0.0;
            ad2NdXi2( 2,  0 ) = xi*(2.0*xi-1.0);
            ad2NdXi2( 3,  0 ) = 0.0;
            ad2NdXi2( 4,  0 ) = 0.5-zeta+2.0*xi*(2.0*zeta-1.0);
            ad2NdXi2( 5,  0 ) = 0.0;

            ad2NdXi2( 0,  1 ) = 0.0;
            ad2NdXi2( 1,  1 ) = 2.0*(zeta-1.0)*zeta;
            ad2NdXi2( 2,  1 ) = eta*(2.0*eta-1.0);
            ad2NdXi2( 3,  1 ) = 0.5-zeta+2.0*eta*(2.0*zeta-1.0);
            ad2NdXi2( 4,  1 ) = 0.0;
            ad2NdXi2( 5,  1 ) = 0.0;

            ad2NdXi2( 0,  2 ) = 2.0*zeta*(zeta-1.0);
            ad2NdXi2( 1,  2 ) = 2.0*zeta*(zeta-1.0);
            ad2NdXi2( 2,  2 ) = 2.0*(1.0-xi-eta)*(0.5-eta-xi);
            ad2NdXi2( 3,  2 ) = 1.5-3.0*zeta+2.0*(xi+eta)*(2.0*zeta-1.0);
            ad2NdXi2( 4,  2 ) = 1.5-3.0*zeta+2.0*(xi+eta)*(2.0*zeta-1.0);
            ad2NdXi2( 5,  2 ) = 2.0*zeta*(zeta-1.0);

            ad2NdXi2( 0,  3 ) = 2.0*zeta*(1.0+zeta);
            ad2NdXi2( 1,  3 ) = 0.0;
            ad2NdXi2( 2,  3 ) = xi*(2.0*xi-1.0);
            ad2NdXi2( 3,  3 ) = 0.0;
            ad2NdXi2( 4,  3 ) = 2.0*xi*(1.0+2.0*zeta)-zeta-0.5;
            ad2NdXi2( 5,  3 ) = 0.0;

            ad2NdXi2( 0,  4 ) = 0.0;
            ad2NdXi2( 1,  4 ) = 2.0*zeta*(1.0+zeta);
            ad2NdXi2( 2,  4 ) = eta*(2.0*eta-1.0);
            ad2NdXi2( 3,  4 ) = 2.0*eta*(1.0+2.0*zeta)-0.5-zeta;
            ad2NdXi2( 4,  4 ) = 0.0;
            ad2NdXi2( 5,  4 ) = 0.0;

            ad2NdXi2( 0,  5 ) = 2.0*zeta*(1.0+zeta);
            ad2NdXi2( 1,  5 ) = 2.0*zeta*(1.0+zeta);
            ad2NdXi2( 2,  5 ) = 2.0*(1.0-xi-eta)*(0.5-xi-eta);
            ad2NdXi2( 3,  5 ) = 2.0*(xi+eta)*(1.0+2.0*zeta)-3.0*zeta-1.5;
            ad2NdXi2( 4,  5 ) = 2.0*(xi+eta)*(1.0+2.0*zeta)-3.0*zeta-1.5;
            ad2NdXi2( 5,  5 ) = 2.0*zeta*(1.0+zeta);

            ad2NdXi2( 0,  6 ) = 0.0;
            ad2NdXi2( 1,  6 ) = 0.0;
            ad2NdXi2( 2,  6 ) = 4.0*xi*eta;
            ad2NdXi2( 3,  6 ) = 2.0*xi*(2.0*zeta-1.0);
            ad2NdXi2( 4,  6 ) = 2.0*eta*(2.0*zeta-1.0);
            ad2NdXi2( 5,  6 ) = 2.0*zeta*(zeta-1.0);

            ad2NdXi2( 0,  7 ) = 0.0;
            ad2NdXi2( 1,  7 ) = 4.0*(1.0-zeta)*zeta;
            ad2NdXi2( 2,  7 ) = 4.0*eta*(1.0-xi-eta);
            ad2NdXi2( 3,  7 ) = 2.0*(xi+2.0*eta)*(1.0-2.0*zeta)+4.0*zeta-2.0;
            ad2NdXi2( 4,  7 ) = eta*(2.0-4.0*zeta);
            ad2NdXi2( 5,  7 ) = 2.0*(1.0-zeta)*zeta;

            ad2NdXi2( 0,  8 ) = 4.0*(1.0-zeta)*zeta;
            ad2NdXi2( 1,  8 ) = 0.0;
            ad2NdXi2( 2,  8 ) = 4.0*xi*(1.0-xi-eta);
            ad2NdXi2( 3,  8 ) = 2.0*xi*(1.0-2.0*zeta);
            ad2NdXi2( 4,  8 ) = 4.0*zeta+2.0*(1.0-2.0*zeta)*(2.0*xi+eta)-2.0;
            ad2NdXi2( 5,  8 ) = 2.0*(1.0-zeta)*zeta;

            ad2NdXi2( 0,  9 ) = 4.0*(1.0-zeta2);
            ad2NdXi2( 1,  9 ) = 0.0;
            ad2NdXi2( 2,  9 ) = (2.0-4.0*xi)*xi;
            ad2NdXi2( 3,  9 ) = 0.0;
            ad2NdXi2( 4,  9 ) = (2.0-8.0*xi)*zeta;
            ad2NdXi2( 5,  9 ) = 0.0;

            ad2NdXi2( 0, 10 ) = 0.0;
            ad2NdXi2( 1, 10 ) = 4.0*(1.0-zeta2);
            ad2NdXi2( 2, 10 ) = (2.0-4.0*eta)*eta;
            ad2NdXi2( 3, 10 ) = (2.0-8.0*eta)*zeta;
            ad2NdXi2( 4, 10 ) = 0.0;
            ad2NdXi2( 5, 10 ) = 0.0;

            ad2NdXi2( 0, 11 ) = 4.0*(1.0-zeta2);
            ad2NdXi2( 1, 11 ) = 4.0*(1.0-zeta2);
            ad2NdXi2( 2, 11 ) = 4.0*(1.0-xi-eta)*(xi+eta-0.5);
            ad2NdXi2( 3, 11 ) = (6.0-8.0*(xi+eta))*zeta;
            ad2NdXi2( 4, 11 ) = (6.0-8.0*(xi+eta))*zeta;
            ad2NdXi2( 5, 11 ) = 4.0*(1.0-zeta2);

            ad2NdXi2( 0, 12 ) = 0.0;
            ad2NdXi2( 1, 12 ) = 0.0;
            ad2NdXi2( 2, 12 ) = 4.0*xi*eta;
            ad2NdXi2( 3, 12 ) = xi*(2.0+4.0*zeta);
            ad2NdXi2( 4, 12 ) = eta*(2.0+4.0*zeta);
            ad2NdXi2( 5, 12 ) = 2.0*zeta*(1.0+zeta);

            ad2NdXi2( 0, 13 ) = 0.0;
            ad2NdXi2( 1, 13 ) = -4.0*zeta*(1.0+zeta);
            ad2NdXi2( 2, 13 ) = 4.0*eta*(1.0-xi-eta);
            ad2NdXi2( 3, 13 ) =2.0*(1.0+2.0*zeta)*(1.0-xi-2.0*eta);
            ad2NdXi2( 4, 13 ) = -2.0*eta*(1.0+2.0*zeta);
            ad2NdXi2( 5, 13 ) = -2.0*zeta*(1.0+zeta);

            ad2NdXi2( 0, 14 ) = -4.0*zeta*(1.0+zeta);
            ad2NdXi2( 1, 14 ) = 0.0;
            ad2NdXi2( 2, 14 ) = 4.0*xi*(1.0-xi-eta);
            ad2NdXi2( 3, 14 ) = -xi*(2.0+4.0*zeta);
            ad2NdXi2( 4, 14 ) = 2.0*(1.0+2.0*zeta)*(1.0-2.0*xi-eta);
            ad2NdXi2( 5, 14 ) = -2.0*zeta*(1.0+zeta);

            ad2NdXi2( 0, 15 ) = 0.0;
            ad2NdXi2( 1, 15 ) = 0.0;
            ad2NdXi2( 2, 15 ) = -8.0*xi*eta;
            ad2NdXi2( 3, 15 ) = -8.0*xi*zeta;
            ad2NdXi2( 4, 15 ) = -8.0*eta*zeta;
            ad2NdXi2( 5, 15 ) = 4.0-4.0*zeta2;

            ad2NdXi2( 0, 16 ) = 0.0;
            ad2NdXi2( 1, 16 ) = 8.0*zeta2-8.0;
            ad2NdXi2( 2, 16 ) = 8.0*eta*(xi+eta-1.0);
            ad2NdXi2( 3, 16 ) = (8.0*xi+16.0*eta-8.0)*zeta;
            ad2NdXi2( 4, 16 ) = 8.0*eta*zeta;
            ad2NdXi2( 5, 16 ) = 4.0*zeta2-4.0;

            ad2NdXi2( 0, 17 ) = 8.0*zeta2-8.0;
            ad2NdXi2( 1, 17 ) = 0.0;
            ad2NdXi2( 2, 17 ) = 8.0*xi*(xi+eta-1.0);
            ad2NdXi2( 3, 17 ) = 8.0*xi*zeta;
            ad2NdXi2( 4, 17 ) = (16.0*xi+8.0*eta-8.0)*zeta;
            ad2NdXi2( 5, 17 ) = 4.0*zeta2-4.0;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_PENTA18_HPP
