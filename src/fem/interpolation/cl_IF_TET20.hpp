//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TET20_HPP
#define BELFEM_CL_IF_TET20_HPP

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >
        ::element_type() const
        {
            return ElementType::TET20;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 20 );

            const real a = 1.0/3.0;
            const real b = 2.0/3.0;

            aXiHat( 0,  0 ) = 1.0;
            aXiHat( 1,  0 ) = 0.0;
            aXiHat( 2,  0 ) = 0.0;

			aXiHat( 0,  1 ) = 0.0;
            aXiHat( 1,  1 ) = 0.0;
            aXiHat( 2,  1 ) = 1.0;
            
            aXiHat( 0,  2 ) = 0.0;
            aXiHat( 1,  2 ) = 1.0;
            aXiHat( 2,  2 ) = 0.0;

            aXiHat( 0,  3 ) = 0.0;
            aXiHat( 1,  3 ) = 0.0;
            aXiHat( 2,  3 ) = 0.0;

			aXiHat( 0,  4 ) = b;
            aXiHat( 1,  4 ) = 0.0;
            aXiHat( 2,  4 ) = a;
            
            aXiHat( 0,  5 ) = a;
            aXiHat( 1,  5 ) = 0.0;
            aXiHat( 2,  5 ) = b;
 
			aXiHat( 0,  6 ) = 0.0;
            aXiHat( 1,  6 ) = a;
            aXiHat( 2,  6 ) = b;
            
            aXiHat( 0,  7 ) = 0.0;
            aXiHat( 1,  7 ) = b;
            aXiHat( 2,  7 ) = a;

            aXiHat( 0,  8 ) = a;
            aXiHat( 1,  8 ) = b;
			aXiHat( 2,  8 ) = 0.0;
			
			aXiHat( 0,  9 ) = b;
            aXiHat( 1,  9 ) = a;
            aXiHat( 2,  9 ) = 0.0;

            aXiHat( 0, 10 ) = a;
            aXiHat( 1, 10 ) = 0.0;
            aXiHat( 2, 10 ) = 0.0;

            aXiHat( 0, 11 ) = b;
            aXiHat( 1, 11 ) = 0.0;
            aXiHat( 2, 11 ) = 0.0;

			aXiHat( 0, 12 ) = 0.0;
            aXiHat( 1, 12 ) = a;
            aXiHat( 2, 12 ) = 0.0;

            aXiHat( 0, 13 ) = 0.0;
            aXiHat( 2, 13 ) = b;
            aXiHat( 3, 13 ) = 0.0;
            
            aXiHat( 0, 14 ) = 0.0;
            aXiHat( 1, 14 ) = 0.0;
            aXiHat( 2, 14 ) = a;

            aXiHat( 0, 15 ) = 0.0;
            aXiHat( 1, 15 ) = 0.0;
            aXiHat( 2, 15 ) = b;

            aXiHat( 0, 16 ) = a;
            aXiHat( 1, 16 ) = a;
            aXiHat( 2, 16 ) = a;

			aXiHat( 0, 17 ) = a;
            aXiHat( 1, 17 ) = 0.0;
            aXiHat( 2, 17 ) = a;
            
            aXiHat( 0, 18 ) = a;
            aXiHat( 1, 18 ) = a;
            aXiHat( 2, 18 ) = 0.0;

            aXiHat( 0, 19 ) = 0.0;
            aXiHat( 1, 19 ) = a;
            aXiHat( 2, 19 ) = a;
        }
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 ) ;

            const real tau = 1.0 - xi - eta - zeta;
            const real psi = 3.0*(xi+eta+zeta);

            aN.set_size( 1, 20 );

            aN( 0,  0 ) = (xi*(9.0*xi*(xi-1.0)+ 2.0))* 0.5;
            aN( 0,  1 ) = (zeta*(9.0*zeta*(zeta-1.0) + 2.0)) * 0.5;
            aN( 0,  2 ) = (eta*(9.0*eta*(eta-1.0)+2.0)) * 0.5;
            aN( 0,  3 ) = 0.5 * tau*( psi-2.0 )*( psi-1.0 );
            aN( 0,  4 ) = 4.5*xi*zeta*(3.0*xi-1.0);
            aN( 0,  5 ) = 4.5*xi*zeta*(3.0*zeta-1.0);
            aN( 0,  6 ) = 4.5*eta*zeta*(3.0*zeta-1.0);
            aN( 0,  7 ) = 4.5*eta*zeta*(3.0*eta-1.0);
            aN( 0,  8 ) = 4.5*eta*xi*(3.0*eta-1.0);
            aN( 0,  9 ) = 4.5*eta*xi*(3.0*xi-1.0);
            aN( 0, 10 ) = 4.5*xi*tau*( 2.0-psi );
            aN( 0, 11 ) = 4.5*xi*(3.0*xi-1.0)*tau;
            aN( 0, 12 ) = 4.5*eta*tau*(2.0-psi);
            aN( 0, 13 ) = 4.5*eta*(3.0*eta-1.0)*tau;
            aN( 0, 14 ) = 4.5*zeta*tau*(2.0-psi);
            aN( 0, 15 ) = 4.5*zeta*(3.0*zeta-1.0)*tau;
			aN( 0, 16 ) = 27.0*xi*eta*zeta;
            aN( 0, 17 ) = 27.0*xi*zeta*tau;
            aN( 0, 18 ) = 27.0*xi*eta*tau;
            aN( 0, 19 ) = 27.0*eta*zeta*tau;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real   xi2 = xi*xi;
            const real  eta2 = eta*eta;

            adNdXi.set_size( 3, 20 );

            adNdXi( 0,  0 ) = 1.0+xi*(13.5*xi-9.0);
            adNdXi( 1,  0 ) = 0.0;
            adNdXi( 2,  0 ) = 0.0;

			adNdXi( 0,  1 ) = 0.0;
            adNdXi( 1,  1 ) = 0.0;
            adNdXi( 2,  1 ) = 1.0+zeta*(13.5*zeta-9.0);
            
            adNdXi( 0,  2 ) = 0.0;
            adNdXi( 1,  2 ) = 1.0+eta*(13.5*eta-9.0);
            adNdXi( 2,  2 ) = 0.0;

            adNdXi( 0,  3 ) = xi*(18.0-27.0*zeta)-13.5*(xi2+eta2)-5.5+eta*(18.0-27.0*(xi+zeta))+zeta*(18.0-13.5*zeta);
            adNdXi( 1,  3 ) = xi*(18.0-27.0*zeta)-13.5*(xi2+eta2)-5.5+eta*(18.0-27.0*(xi+zeta))+zeta*(18.0-13.5*zeta);
            adNdXi( 2,  3 ) = xi*(18.0-27.0*zeta)-13.5*(xi2+eta2)-5.5+eta*(18.0-27.0*(xi+zeta))+zeta*(18.0-13.5*zeta);

			adNdXi( 0,  4 ) = zeta*(27.0*xi-4.5);
            adNdXi( 1,  4 ) = 0.0;
            adNdXi( 2,  4 ) = xi*(13.5*xi-4.5);
            
            adNdXi( 0,  5 ) = zeta*(13.5*zeta-4.5);
            adNdXi( 1,  5 ) = 0.0;
            adNdXi( 2,  5 ) = xi*(27.0*zeta-4.5);

			adNdXi( 0,  6 ) = 0.0;
            adNdXi( 1,  6 ) = zeta*(13.5*zeta-4.5);
            adNdXi( 2,  6 ) = eta*(27.0*zeta-4.5);
            
			adNdXi( 0,  7 ) = 0.0;
            adNdXi( 1,  7 ) = zeta*(27.0*eta-4.5);
            adNdXi( 2,  7 ) = eta*(13.5*eta-4.5);
            
            adNdXi( 0,  8 ) = eta*(13.5*eta-4.5);
            adNdXi( 1,  8 ) = xi*(27.0*eta-4.5);
            adNdXi( 2,  8 ) = 0.0;
            
            adNdXi( 0,  9 ) = eta*(27.0*xi-4.5);
            adNdXi( 1,  9 ) = xi*(13.5*xi-4.5);
            adNdXi( 2,  9 ) = 0.0;

            adNdXi( 0, 10 ) = 9.0+13.5*eta2+40.5*xi2+zeta*(13.5*zeta-22.5)+eta*(54.0*xi+27.0*zeta-22.5)+xi*(54.0*zeta-45.0);
            adNdXi( 1, 10 ) = xi*(27.0*(xi+eta+zeta)-22.5);
            adNdXi( 2, 10 ) = xi*(27.0*(xi+eta+zeta)-22.5);

            adNdXi( 0, 11 ) = xi*(36.0-40.5*xi-27.0*zeta)+eta*(4.5-27.0*xi)+4.5*zeta-4.5;
            adNdXi( 1, 11 ) = (4.5-13.5*xi)*xi;
            adNdXi( 2, 11 ) = (4.5-13.5*xi)*xi;

			adNdXi( 0, 12 ) = eta*(27.0*(xi+eta+zeta)-22.5);
            adNdXi( 1, 12 ) = 9.0+40.5*eta2+13.5*xi2+xi*(27.0*zeta-22.5)+eta*(54.0*xi+54.0*zeta-45.0)+zeta*(13.5*zeta-22.5);
            adNdXi( 2, 12 ) = eta*(27.0*(xi+eta+zeta)-22.5);

			adNdXi( 0, 13 ) = eta*(4.5-13.5*eta);
            adNdXi( 1, 13 ) = 4.5*(xi+zeta)+eta*(36.0-40.5*eta-27.0*(xi+zeta))-4.5;
            adNdXi( 2, 13 ) = eta*(4.5-13.5*eta);

            adNdXi( 0, 14 ) = zeta*(27.0*(xi+eta+zeta)-22.5);
            adNdXi( 1, 14 ) = zeta*(27.0*(xi+eta+zeta)-22.5);
            adNdXi( 2, 14 ) = 9.0+13.5*(xi2+eta2)+zeta*(40.5*zeta-45.0)+xi*(54.0*zeta-22.5)+eta*(27.0*xi+54.0*zeta-22.5);

            adNdXi( 0, 15 ) = zeta*(4.5-13.5*zeta);
            adNdXi( 1, 15 ) = zeta*(4.5-13.5*zeta);
            adNdXi( 2, 15 ) = xi*(4.5-27.0*zeta)+eta*(4.5-27.0*zeta)+zeta*(36.0-40.5*zeta)-4.5;

            adNdXi( 0, 16 ) = 27.0*eta*zeta;
            adNdXi( 1, 16 ) = 27.0*xi*zeta;
            adNdXi( 2, 16 ) = 27.0*xi*eta;

			adNdXi( 0, 17 ) = 27.0*zeta*(1.0-(2.0*xi+eta+zeta));
            adNdXi( 1, 17 ) = -27.0*xi*zeta;
            adNdXi( 2, 17 ) = 27.0*xi*(1.0-(xi+eta+2.0*zeta));
            
            adNdXi( 0, 18 ) = 27.0*eta*(1.0 - (eta+zeta+2.0*xi));
            adNdXi( 1, 18 ) = 27.0*xi*(1.0-(xi+2.0*eta+zeta));
            adNdXi( 2, 18 ) = -27.0*xi*eta;

            adNdXi( 0, 19 ) = -27.0*eta*zeta;
            adNdXi( 1, 19 ) = 27.0*zeta*(1.0-(xi+2.0*eta+zeta));
            adNdXi( 2, 19 ) = 27.0*eta*(1.0-(xi+eta+2.0*zeta));
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            ad2NdXi2.set_size( 6, 20 );

            ad2NdXi2( 0,  0 ) = 27.0*xi-9.0;
            ad2NdXi2( 1,  0 ) = 0.0;
            ad2NdXi2( 2,  0 ) = 0.0;
            ad2NdXi2( 3,  0 ) = 0.0;
            ad2NdXi2( 4,  0 ) = 0.0;
            ad2NdXi2( 5,  0 ) = 0.0;

			ad2NdXi2( 0,  1 ) = 0.0;
            ad2NdXi2( 1,  1 ) = 0.0;
            ad2NdXi2( 2,  1 ) = 27.0*zeta-9.0;
            ad2NdXi2( 3,  1 ) = 0.0;
            ad2NdXi2( 4,  1 ) = 0.0;
            ad2NdXi2( 5,  1 ) = 0.0;
            
            ad2NdXi2( 0,  2 ) = 0.0;
            ad2NdXi2( 1,  2 ) = 27.0*eta-9.0;
            ad2NdXi2( 2,  2 ) = 0.0;
            ad2NdXi2( 3,  2 ) = 0.0;
            ad2NdXi2( 4,  2 ) = 0.0;
            ad2NdXi2( 5,  2 ) = 0.0;

            ad2NdXi2( 0,  3 ) = 18.0-27.0*(xi+eta+zeta);
            ad2NdXi2( 1,  3 ) = 18.0-27.0*(xi+eta+zeta);
            ad2NdXi2( 2,  3 ) = 18.0-27.0*(xi+eta+zeta);
            ad2NdXi2( 3,  3 ) = 18.0-27.0*(xi+eta+zeta);
            ad2NdXi2( 4,  3 ) = 18.0-27.0*(xi+eta+zeta);
            ad2NdXi2( 5,  3 ) = 18.0-27.0*(xi+eta+zeta);

			ad2NdXi2( 0,  4 ) = 27.0*zeta;
            ad2NdXi2( 1,  4 ) = 0.0;
            ad2NdXi2( 2,  4 ) = 0.0;
            ad2NdXi2( 3,  4 ) = 0.0;
            ad2NdXi2( 4,  4 ) = 27.0*xi-4.5;
            ad2NdXi2( 5,  4 ) = 0.0;
            
            ad2NdXi2( 0,  5 ) = 0.0;
            ad2NdXi2( 1,  5 ) = 0.0;
            ad2NdXi2( 2,  5 ) = 27.0*xi;
            ad2NdXi2( 3,  5 ) = 0.0;
            ad2NdXi2( 4,  5 ) = 27.0*zeta-4.5;
            ad2NdXi2( 5,  5 ) = 0.0;
            
            ad2NdXi2( 0,  6 ) = 0.0;
            ad2NdXi2( 1,  6 ) = 0.0;
            ad2NdXi2( 2,  6 ) = 27.0*eta;
            ad2NdXi2( 3,  6 ) = 27.0*zeta-4.5;
            ad2NdXi2( 4,  6 ) = 0.0;
            ad2NdXi2( 5,  6 ) = 0.0;

            ad2NdXi2( 0,  7 ) = 0.0;
            ad2NdXi2( 1,  7 ) = 27.0*zeta;
            ad2NdXi2( 2,  7 ) = 0.0;
            ad2NdXi2( 3,  7 ) = 27.0*eta-4.5;
            ad2NdXi2( 4,  7 ) = 0.0;
            ad2NdXi2( 5,  7 ) = 0.0;

            ad2NdXi2( 0,  8 ) = 0.0;
            ad2NdXi2( 1,  8 ) = 27.0*xi;
            ad2NdXi2( 2,  8 ) = 0.0;
            ad2NdXi2( 3,  8 ) = 0.0;
            ad2NdXi2( 4,  8 ) = 0.0;
            ad2NdXi2( 5,  8 ) = 27.0*eta-4.5;

            ad2NdXi2( 0,  9 ) = 27.0*eta;
            ad2NdXi2( 1,  9 ) = 0.0;
            ad2NdXi2( 2,  9 ) = 0.0;
            ad2NdXi2( 3,  9 ) = 0.0;
            ad2NdXi2( 4,  9 ) = 0.0;
            ad2NdXi2( 5,  9 ) = 27.0*xi-4.5;


            ad2NdXi2( 0, 10 ) = 81.0*xi+54.0*(eta+zeta)-45.0;
            ad2NdXi2( 1, 10 ) = 27.0*xi;
            ad2NdXi2( 2, 10 ) = 27.0*xi;
            ad2NdXi2( 3, 10 ) = 27.0*xi;
            ad2NdXi2( 4, 10 ) = 54.0*xi+27.0*(eta+zeta)-22.5;
            ad2NdXi2( 5, 10 ) = 54.0*xi+27.0*(eta+zeta)-22.5;

            ad2NdXi2( 0, 11 ) = 36.0-81.0*xi-27.0*(eta+zeta);
            ad2NdXi2( 1, 11 ) = 0.0;
            ad2NdXi2( 2, 11 ) = 0.0;
            ad2NdXi2( 3, 11 ) = 0.0;
            ad2NdXi2( 4, 11 ) = 4.5-27.0*xi;
            ad2NdXi2( 5, 11 ) = 4.5-27.0*xi;

			ad2NdXi2( 0, 12 ) = 27.0*eta;
            ad2NdXi2( 1, 12 ) = 54.0*(xi+zeta)+81.0*eta-45.0;
            ad2NdXi2( 2, 12 ) = 27.0*eta;
            ad2NdXi2( 3, 12 ) = 27.0*(xi+zeta)+54.0*eta-22.5;
            ad2NdXi2( 4, 12 ) = 27.0*eta;
            ad2NdXi2( 5, 12 ) = 27.0*(xi+zeta)+54.0*eta-22.5;

            ad2NdXi2( 0, 13 ) = 0.0;
            ad2NdXi2( 1, 13 ) = 36.0-27.0*(xi+zeta)-81.0*eta;
            ad2NdXi2( 2, 13 ) = 0.0;
            ad2NdXi2( 3, 13 ) = 4.5-27.0*eta;
            ad2NdXi2( 4, 13 ) = 0.0;
            ad2NdXi2( 5, 13 ) = 4.5-27.0*eta;
            
            ad2NdXi2( 0, 14 ) = 27.0*zeta;
            ad2NdXi2( 1, 14 ) = 27.0*zeta;
            ad2NdXi2( 2, 14 ) = 54.0*(xi+eta)+81.0*zeta-45.0;
            ad2NdXi2( 3, 14 ) = 27.0*(xi+eta)+54.0*zeta-22.5;
            ad2NdXi2( 4, 14 ) = 27.0*(xi+eta)+54.0*zeta-22.5;
            ad2NdXi2( 5, 14 ) = 27.0*zeta;

            ad2NdXi2( 0, 15 ) = 0.0;
            ad2NdXi2( 1, 15 ) = 0.0;
            ad2NdXi2( 2, 15 ) = 36.0-27.0*(xi+eta)-81.0*zeta;
            ad2NdXi2( 3, 15 ) = 4.5-27.0*zeta;
            ad2NdXi2( 4, 15 ) = 4.5-27.0*zeta;
            ad2NdXi2( 5, 15 ) = 0.0;

            ad2NdXi2( 0, 16 ) = 0.0;
            ad2NdXi2( 1, 16 ) = 0.0;
            ad2NdXi2( 2, 16 ) = 0.0;
            ad2NdXi2( 3, 16 ) = 27.0*xi;
            ad2NdXi2( 4, 16 ) = 27.0*eta;
            ad2NdXi2( 5, 16 ) = 27.0*zeta;

			ad2NdXi2( 0, 17 ) = -54.0*zeta;
            ad2NdXi2( 1, 17 ) = 0.0;
            ad2NdXi2( 2, 17 ) = -54.0*xi;
            ad2NdXi2( 3, 17 ) = -27.0*xi;
            ad2NdXi2( 4, 17 ) = 27.0-54.0*(xi+zeta)-27.0*eta;
            ad2NdXi2( 5, 17 ) = -27.0*zeta;

            ad2NdXi2( 0, 18 ) = -54.0*eta;
            ad2NdXi2( 1, 18 ) = -54.0*xi;
            ad2NdXi2( 2, 18 ) = 0.0;
            ad2NdXi2( 3, 18 ) = -27.0*xi;
            ad2NdXi2( 4, 18 ) = -27.0*eta;
            ad2NdXi2( 5, 18 ) = 27.0-54.0*(xi+eta)-27.0*zeta;

            ad2NdXi2( 0, 19 ) = 0.0;
            ad2NdXi2( 1, 19 ) = -54.0*zeta;
            ad2NdXi2( 2, 19 ) = -54.0*eta;
            ad2NdXi2( 3, 19 ) = 27.0-27.0*xi-54.0*(eta+zeta);
            ad2NdXi2( 4, 19 ) = -27.0*eta;
            ad2NdXi2( 5, 19 ) = -27.0*zeta;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_TET20_HPP
