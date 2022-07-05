//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TRI10_HPP
#define BELFEM_CL_IF_TRI10_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >
        ::element_type() const
        {
            return ElementType::TRI10;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            const real a = 1.0/3.0;
            const real b = 2.0/3.0;

            aXiHat.set_size( 2, 10 );

            aXiHat( 0, 0 ) = 1.0;
            aXiHat( 1, 0 ) = 0.0;

            aXiHat( 0, 1 ) = 0.0;
            aXiHat( 1, 1 ) = 1.0;

            aXiHat( 0, 2 ) = 0.0;
            aXiHat( 1, 2 ) = 0.0;

            aXiHat( 0, 3 ) = a;
            aXiHat( 1, 3 ) = b;

            aXiHat( 0, 4 ) = b;
            aXiHat( 1, 4 ) = a;

            aXiHat( 0, 5 ) = 0;
            aXiHat( 1, 5 ) = a;

            aXiHat( 0, 6 ) = 0;
            aXiHat( 1, 6 ) = b;

            aXiHat( 0, 7 ) = a;
            aXiHat( 1, 7 ) = 0;

            aXiHat( 0, 8 ) = b;
            aXiHat( 1, 8 ) = 0;

            aXiHat( 0, 9 ) = a;
            aXiHat( 1, 9 ) = a;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >::N(
                const Vector< real > & aXi,
                      Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = 1.0 - xi - eta;

            aN.set_size( 1, 10 );

            aN( 0, 0 ) = (9.0*(xi-1)*xi+2.0)*xi*0.5;

            aN( 0, 1 ) = 0.5*eta*(9.0*eta*(eta-1.0)+2.0);

            aN( 0, 2 ) = 0.5*zeta*(3.0*(xi+eta)-2.0)*(3.0*(xi+eta)-1.0);

            aN( 0, 3 ) = 4.5*eta*xi*(3.0*eta-1.0);

            aN( 0, 4 ) = 4.5*eta*xi*(3.0*xi-1.0);

            aN( 0, 5 ) = 4.5*eta*((3.0*(xi+eta)-5.0)*(xi+eta)+2.0);

            aN( 0, 6 ) = 4.5*zeta*eta*(3.0*eta-1.0);

            aN( 0, 7 ) = 4.5*xi*((3.0*(xi+eta)-5.0)*(xi+eta)+2.0);

            aN( 0, 8 ) = 4.5*xi*zeta*(3.0*xi-1.0);

            aN( 0, 9 ) = 27.0*xi*eta*zeta;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );

            adNdXi.set_size( 2, 10 );

            adNdXi(0,0)=1.0+xi*(13.5*xi-9.0);
            adNdXi(1,0)=0.0;

            adNdXi(0,1)=0.0;
            adNdXi(1,1)=1.0+eta*(13.5*eta-9.0);

            adNdXi(0,2)=eta*((18.0-27.0*xi)-13.5*eta)+(18.0-13.5*xi)*xi-5.5;
            adNdXi(1,2)=eta*((18.0-27.0*xi)-13.5*eta)+(18.0-13.5*xi)*xi-5.5;

            adNdXi(0,3)=eta*(13.5*eta-4.5);
            adNdXi(1,3)=xi*(27.0*eta-4.5);

            adNdXi(0,4)=eta*(27.*xi-4.5);
            adNdXi(1,4)=xi*(13.5*xi-4.5);

            adNdXi(0,5)=eta*(27.0*(eta+xi)-22.5);
            adNdXi(1,5)=9.0+eta*(40.5*eta+54.0*xi-45.0)+xi*(13.5*xi-22.5);

            adNdXi(0,6)=(4.5-13.5*eta)*eta;
            adNdXi(1,6)=eta*(36.0-40.5*eta-27.0*xi)+4.5*xi-4.5;

            adNdXi(0,7)=xi*(40.5*xi-45.0)+eta*(13.5*eta+54.0*xi-22.5)+9.0;
            adNdXi(1,7)=xi*(27.0*(xi+eta)-22.5);

            adNdXi(0,8)=(36.0-40.5*xi)*xi+eta*(4.5-27.0*xi)-4.5;
            adNdXi(1,8)=(4.5-13.5*xi)*xi;

            adNdXi(0,9)=27.0*eta*(1.0-eta-2.0*xi);
            adNdXi(1,9)=27.0*xi*(1.0-2.0*eta-xi);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );

            ad2NdXi2.set_size( 3, 10 );

            ad2NdXi2(0,0)=27.0*xi-9.0;
            ad2NdXi2(1,0)=0.0;
            ad2NdXi2(2,0)=0.0;

            ad2NdXi2(0,1)=0.0;
            ad2NdXi2(1,1)=27.0*eta-9.0;
            ad2NdXi2(2,1)=0.0;

            ad2NdXi2(0,2)=18.0-27.0*(xi+eta);
            ad2NdXi2(1,2)=18.0-27.0*(xi+eta);
            ad2NdXi2(2,2)=18.0-27.0*(xi+eta);

            ad2NdXi2(0,3)=0.0;
            ad2NdXi2(1,3)=27.0*xi;
            ad2NdXi2(2,3)=27.0*eta-4.5;

            ad2NdXi2(0,4)=27.0*eta;
            ad2NdXi2(1,4)=0.0;
            ad2NdXi2(2,4)=27.0*xi-4.5;

            ad2NdXi2(0,5)=27.0*eta;
            ad2NdXi2(1,5)=54.0*xi+81.0*eta-45.0;
            ad2NdXi2(2,5)=27.0*xi+54.0*eta-22.5;

            ad2NdXi2(0,6)=0.0;
            ad2NdXi2(1,6)=36.0-81.0*eta-27.0*xi;
            ad2NdXi2(2,6)=4.5-27.0*eta;

            ad2NdXi2(0,7)=81.0*xi+54.0*eta-45.0;
            ad2NdXi2(1,7)=27.0*xi;
            ad2NdXi2(2,7)=54.0*xi+27.0*eta-22.5;

            ad2NdXi2(0,8)=36.-27.0*eta-81.0*xi;
            ad2NdXi2(1,8)=0.0;
            ad2NdXi2(2,8)=4.5-27.0*xi;

            ad2NdXi2(0,9)=-54.0*eta;
            ad2NdXi2(1,9)=-54.0*xi;
            ad2NdXi2(2,9)=27.0*(1.0-2.0*(xi+eta));
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TRI10_HPP
