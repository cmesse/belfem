//
// Created by Christian Messe on 26.10.19.
//

#ifndef BELFEM_CL_IF_PLATE_HPP
#define BELFEM_CL_IF_PLATE_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >
        ::element_type() const
        {
            return ElementType::QUAD4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >
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
                GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real  a = xi-1.0;
            const real  b = xi+1.0;
            const real  c = eta-1.0;
            const real  d = eta+1.0;
            const real  e = xi-2.0;
            const real  f = xi+2.0;
            const real  g = eta-2.0;
            const real  h = eta+2.0;
            const real  a2 = a*a;
            const real  b2 = b*b;
            const real  c2 = c*c;
            const real  d2 = d*d;
            aN.set_size( 1, 16 );

            aN( 0,  0 ) =  c2*h*a2*f;
            aN( 0,  1 ) =  c2*h*a2*b;
            aN( 0,  2 ) =  c2*d*a2*f;
            aN( 0,  3 ) =  c2*d*a2*b;
            aN( 0,  4 ) = -c2*h*b2*e;
            aN( 0,  5 ) =  c2*h*a*b2;
            aN( 0,  6 ) = -c2*d*b2*e;
            aN( 0,  7 ) =  c2*d*a*b2;
            aN( 0,  8 ) =  d2*g*b2*e;
            aN( 0,  9 ) = -d2*g*a*b2;
            aN( 0, 10 ) = -c*d2*b2*e;
            aN( 0, 11 ) =  c*d2*a*b2;
            aN( 0, 12 ) = -d2*g*a2*f;
            aN( 0, 13 ) = -d2*g*a2*b;
            aN( 0, 14 ) =  c*d2*a2*f;
            aN( 0, 15 ) =  c*d2*a2*b;

            aN *= 0.0625;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real  a = xi-1.0;
            const real  b = xi+1.0;
            const real  c = eta-1.0;
            const real  d = eta+1.0;
            const real  e = xi-2.0;
            const real  f = xi+2.0;
            const real  g = eta-2.0;
            const real  h = eta+2.0;
            const real  p = a*b;
            const real  q = c*d;
            const real  r = 1.0+xi*(2.0-3.0*xi);
            const real  s = (3.0*xi+2.0)*xi-1.0;
            const real  u = 1.0 + eta*( 2.0-3.0*eta);
            const real  v = (3.0*eta+2.0)*eta - 1.0;
            const real  a2 = a*a;
            const real  b2 = b*b;
            const real  c2 = c*c;
            const real  d2 = d*d;

            adNdXi.set_size( 2, 16 );

            adNdXi( 0,  0 ) =  3.0*p*c2*h;
            adNdXi( 1,  0 ) =  3.0*q*a2*f;

            adNdXi( 0,  1 ) = -c2*h*r;
            adNdXi( 1,  1 ) =  3.0*q*a2*b;

            adNdXi( 0,  2 ) =  3.0*p*c2*d;
            adNdXi( 1,  2 ) = -a2*f*u;

            adNdXi( 0,  3 ) = -c2*d*r;
            adNdXi( 1,  3 ) = -a2*b*u;

            adNdXi( 0,  4 ) = -3.0*p*c2*h;
            adNdXi( 1,  4 ) = -3.0*q*b2*e;

            adNdXi( 0,  5 ) =  c2*h*s;
            adNdXi( 1,  5 ) =  3.0*q*a*b2;

            adNdXi( 0,  6 ) = -3.0*p*c2*d;
            adNdXi( 1,  6 ) =  b2*e*u;

            adNdXi( 0,  7 ) =  c2*d*s;
            adNdXi( 1,  7 ) = -a*b2*u;

            adNdXi( 0,  8 ) =  3.0*p*d2*g;
            adNdXi( 1,  8 ) =  3.0*q*b2*e;

            adNdXi( 0,  9 ) = -d2*g*s;
            adNdXi( 1,  9 ) = -3.0*q*a*b2;

            adNdXi( 0, 10 ) = -3.0*p*c*d2;
            adNdXi( 1, 10 ) = -b2*e*v;

            adNdXi( 0, 11 ) =  c*d2*s;
            adNdXi( 1, 11 ) =  a*b2*v;

            adNdXi( 0, 12 ) = -3.0*p*d2*g;
            adNdXi( 1, 12 ) = -3.0*q*a2*f;

            adNdXi( 0, 13 ) =  d2*g*r;
            adNdXi( 1, 13 ) = -3.0*q*a2*b;

            adNdXi( 0, 14 ) =  3.0*p*c*d2;
            adNdXi( 1, 14 ) =  a2*f*v;

            adNdXi( 0, 15 ) = -c*d2*r;
            adNdXi( 1, 15 ) =  a2*b*v;

            adNdXi *= 0.0625;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real  a = xi-1.0;
            const real  b = xi+1.0;
            const real  c = eta-1.0;
            const real  d = eta+1.0;
            const real  e = xi-2.0;
            const real  f = xi+2.0;
            const real  g = eta-2.0;
            const real  h = eta+2.0;
            const real  p = a*b;
            const real  q = c*d;
            const real  r = 1.0+xi*(2.0-3.0*xi);
            const real  s = (3.0*xi+2.0)*xi-1.0;
            const real  u = 1.0 + eta*( 2.0-3.0*eta);
            const real  v = (3.0*eta+2.0)*eta - 1.0;

            const real  w = (3.0*xi-1.0);
            const real  x = (3.0*xi+1.0);
            const real  y = (3.0*eta-1.0);
            const real  z = (3.0*eta+1.0);

            const real  a2 = a*a;
            const real  b2 = b*b;
            const real  c2 = c*c;
            const real  d2 = d*d;

            ad2NdXi2.set_size( 3, 16 );

            ad2NdXi2( 0,  0 ) =  6.0*xi*c2*h;
            ad2NdXi2( 1,  0 ) =  6.0*eta*a2*f;
            ad2NdXi2( 2,  0 ) =  9.0*q*p;

            ad2NdXi2( 0,  1 ) =  2.0*w*c2*h;
            ad2NdXi2( 1,  1 ) =  6.0*eta*a2*b;
            ad2NdXi2( 2,  1 ) = -3.0*q*r;

            ad2NdXi2( 0,  2 ) =  6.0*xi*c2*d;
            ad2NdXi2( 1,  2 ) =  2.0*y*a2*f;
            ad2NdXi2( 2,  2 ) = -3.0*p*u;

            ad2NdXi2( 0,  3 ) =  2.0*w*c2*d;
            ad2NdXi2( 1,  3 ) =  2.0*y*a2*b;
            ad2NdXi2( 2,  3 ) =  u*r;

            ad2NdXi2( 0,  4 ) = -6.0*xi*c2*h;
            ad2NdXi2( 1,  4 ) = -6.0*eta*b2*e;
            ad2NdXi2( 2,  4 ) = -9.0*q*p;

            ad2NdXi2( 0,  5 ) =  2.0*x*c2*h;
            ad2NdXi2( 1,  5 ) =  6.0*eta*a*b2;
            ad2NdXi2( 2,  5 ) =  3.0*q*s;

            ad2NdXi2( 0,  6 ) = -6.0*xi*c2*d;
            ad2NdXi2( 1,  6 ) = -2.0*y*b2*e;
            ad2NdXi2( 2,  6 ) =  3.0*p*u;

            ad2NdXi2( 0,  7 ) =  2.0*x*c2*d;
            ad2NdXi2( 1,  7 ) =  2.0*y*a*b2;
            ad2NdXi2( 2,  7 ) = -u*s;

            ad2NdXi2( 0,  8 ) =  6.0*xi*d2*g;
            ad2NdXi2( 1,  8 ) =  6.0*eta*b2*e;
            ad2NdXi2( 2,  8 ) =  9.0*q*p;

            ad2NdXi2( 0,  9 ) = -2.0*x*d2*g;
            ad2NdXi2( 1,  9 ) = -6.0*eta*a*b2;
            ad2NdXi2( 2,  9 ) = -3.0*q*s;

            ad2NdXi2( 0, 10 ) = -6.0*xi*c*d2;
            ad2NdXi2( 1, 10 ) = -2.0*z*b2*e;
            ad2NdXi2( 2, 10 ) = -3.0*p*v;

            ad2NdXi2( 0, 11 ) =  2.0*x*c*d2;
            ad2NdXi2( 1, 11 ) =  2.0*z*a*b2;
            ad2NdXi2( 2, 11 ) =  v*s;

            ad2NdXi2( 0, 12 ) = -6.0*xi*d2*g;
            ad2NdXi2( 1, 12 ) = -6.0*eta*a2*f;
            ad2NdXi2( 2, 12 ) = -9.0*q*p;

            ad2NdXi2( 0, 13 ) = -2.0*w*d2*g;
            ad2NdXi2( 1, 13 ) = -6.0*eta*a2*b;
            ad2NdXi2( 2, 13 ) =  3.0*q*r;

            ad2NdXi2( 0, 14 ) =  6.0*xi*c*d2;
            ad2NdXi2( 1, 14 ) =  2.0*z*a2*f;
            ad2NdXi2( 2, 14 ) =  3.0*p*v;

            ad2NdXi2( 0, 15 ) =  2.0*w*c*d2;
            ad2NdXi2( 1, 15 ) =  2.0*z*a2*b;
            ad2NdXi2( 2, 15 ) = -v*r;

            ad2NdXi2*= 0.0625;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_PLATE_HPP
