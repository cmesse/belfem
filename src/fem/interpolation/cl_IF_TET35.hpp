//
// Created by Christian Messe on 25.07.20.
//

#ifndef BELFEM_CL_IF_TET35_HPP
#define BELFEM_CL_IF_TET35_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUARTIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >
        ::element_type() const
        {
            return ElementType::TET35;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 35 );

            aXiHat( 0 ,  0 ) =  1.00 ;
            aXiHat( 1 ,  0 ) =  0.00 ;
            aXiHat( 2 ,  0 ) =  0.00 ;

            aXiHat( 0 ,  2 ) =  0.00 ;
            aXiHat( 1 ,  2 ) =  1.00 ;
            aXiHat( 2 ,  2 ) =  0.00 ;

            aXiHat( 0 ,  1 ) =  0.00 ;
            aXiHat( 1 ,  1 ) =  0.00 ;
            aXiHat( 2 ,  1 ) =  1.00 ;

            aXiHat( 0 ,  3 ) =  0.00 ;
            aXiHat( 1 ,  3 ) =  0.00 ;
            aXiHat( 2 ,  3 ) =  0.00 ;

			aXiHat( 0 ,  4 ) =  0.75 ;
            aXiHat( 1 ,  4 ) =  0.00 ;
            aXiHat( 2 ,  4 ) =  0.25 ;
            
            aXiHat( 0 ,  5 ) =  0.50 ;
            aXiHat( 1 ,  5 ) =  0.00 ;
            aXiHat( 2 ,  5 ) =  0.50 ;
            
            aXiHat( 0 ,  6 ) =  0.25 ;
            aXiHat( 1 ,  6 ) =  0.00 ;
            aXiHat( 2 ,  6 ) =  0.75 ;
            
            aXiHat( 0 ,  7 ) =  0.00 ;
            aXiHat( 1 ,  7 ) =  0.25 ;
            aXiHat( 2 ,  7 ) =  0.75 ;
            
            aXiHat( 0 ,  8 ) =  0.00 ;
            aXiHat( 1 ,  8 ) =  0.50 ;
            aXiHat( 2 ,  8 ) =  0.50 ;
            
			aXiHat( 0 ,  9 ) =  0.00 ;
            aXiHat( 1 ,  9 ) =  0.75 ;
            aXiHat( 2 ,  9 ) =  0.25 ;

            aXiHat( 0 , 10 ) =  0.25 ;
            aXiHat( 1 , 10 ) =  0.75 ;
            aXiHat( 2 , 10 ) =  0.00 ;

            aXiHat( 0 , 11 ) =  0.50 ;
            aXiHat( 1 , 11 ) =  0.50 ;
            aXiHat( 2 , 11 ) =  0.00 ;

            aXiHat( 0 , 12 ) =  0.75 ;
            aXiHat( 1 , 12 ) =  0.25 ;
            aXiHat( 2 , 12 ) =  0.00 ;

            aXiHat( 0 , 13 ) =  0.25 ;
            aXiHat( 1 , 13 ) =  0.00 ;
            aXiHat( 2 , 13 ) =  0.00 ;

            aXiHat( 0 , 14 ) =  0.50 ;
            aXiHat( 1 , 14 ) =  0.00 ;
            aXiHat( 2 , 14 ) =  0.00 ;

            aXiHat( 0 , 15 ) =  0.75 ;
            aXiHat( 1 , 15 ) =  0.00 ;
            aXiHat( 2 , 15 ) =  0.00 ;

			aXiHat( 0 , 16 ) =  0.00 ;
            aXiHat( 1 , 16 ) =  0.25 ;
            aXiHat( 2 , 16 ) =  0.00 ;

            aXiHat( 0 , 17 ) =  0.00 ;
            aXiHat( 1 , 17 ) =  0.50 ;
            aXiHat( 2 , 17 ) =  0.00 ;

            aXiHat( 0 , 18 ) =  0.00 ;
            aXiHat( 1 , 18 ) =  0.75 ;
            aXiHat( 2 , 18 ) =  0.00 ;
            
            aXiHat( 0 , 19 ) =  0.00 ;
            aXiHat( 1 , 19 ) =  0.00 ;
            aXiHat( 2 , 19 ) =  0.25 ;

            aXiHat( 0 , 20 ) =  0.00 ;
            aXiHat( 1 , 20 ) =  0.00 ;
            aXiHat( 2 , 20 ) =  0.50 ;

            aXiHat( 0 , 21 ) =  0.00 ;
            aXiHat( 1 , 21 ) =  0.00 ;
            aXiHat( 2 , 21 ) =  0.75 ;

            aXiHat( 0 , 22 ) =  0.50 ;
            aXiHat( 1 , 22 ) =  0.25 ;
            aXiHat( 2 , 22 ) =  0.25 ;

			aXiHat( 0 , 23 ) =  0.25 ;
            aXiHat( 1 , 23 ) =  0.50 ;
            aXiHat( 2 , 23 ) =  0.25 ;
            
            aXiHat( 0 , 24 ) =  0.25 ;
            aXiHat( 1 , 24 ) =  0.25 ;
            aXiHat( 2 , 24 ) =  0.50 ;

            
			aXiHat( 0 , 25 ) =  0.50 ;
            aXiHat( 1 , 25 ) =  0.00 ;
            aXiHat( 2 , 25 ) =  0.25 ;

			aXiHat( 0 , 26 ) =  0.25 ;
            aXiHat( 1 , 26 ) =  0.00 ;
            aXiHat( 2 , 26 ) =  0.50 ;

            aXiHat( 0 , 27 ) =  0.25 ;
            aXiHat( 1 , 27 ) =  0.00 ;
            aXiHat( 2 , 27 ) =  0.25 ;
            
            aXiHat( 0 , 28 ) =  0.50 ;
            aXiHat( 1 , 28 ) =  0.25 ;
            aXiHat( 2 , 28 ) =  0.00 ;

			aXiHat( 0 , 29 ) =  0.25 ;
            aXiHat( 1 , 29 ) =  0.25 ;
            aXiHat( 2 , 29 ) =  0.00 ;
            
            aXiHat( 0 , 30 ) =  0.25 ;
            aXiHat( 1 , 30 ) =  0.50 ;
            aXiHat( 2 , 30 ) =  0.00 ;

            aXiHat( 0 , 31 ) =  0.00 ;
            aXiHat( 1 , 31 ) =  0.25 ;
            aXiHat( 2 , 31 ) =  0.25 ;

			aXiHat( 0 , 32 ) =  0.00 ;
            aXiHat( 1 , 32 ) =  0.25 ;
            aXiHat( 2 , 32 ) =  0.50 ;
            
            aXiHat( 0 , 33 ) =  0.00 ;
            aXiHat( 1 , 33 ) =  0.50 ;
            aXiHat( 2 , 33 ) =  0.25 ;

            aXiHat( 0 , 34 ) =  0.25 ;
            aXiHat( 1 , 34 ) =  0.25 ;
            aXiHat( 2 , 34 ) =  0.25 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real a = 1./3.;
            const real b = 4.0;
            const real c = 16.0/3.0;
            const real d = 32.0;

            const real tau = 1.-xi-eta-zeta;
            const real alpha = (3.-4.*(xi+eta+zeta));
            const real beta  = (1.-4.*(xi+eta+zeta));
            const real gamma = (1.-2.*(xi+eta+zeta));
            const real delta = (eta*(8.*eta-6.)+1.);
            const real phi = (xi*(8.*xi-6.)+1.);
            const real psi = (zeta*(8.*zeta-6.)+1.);


            aN.set_size( 1, 35 );

            aN( 0,  0 ) = a*xi*(2.*xi-1.)*(4.*xi-3.)*(4.*xi-1.) ;
            aN( 0,  1 ) = a*zeta*(2.*zeta-1.)*(4.*zeta-3.)*(4.*zeta-1.) ;
            aN( 0,  2 ) = a*eta*(2.*eta-1.)*(4.*eta-3.)*(4.*eta-1.) ;
            aN( 0,  3 ) = a*tau*(2.*tau-1.)*alpha*beta ;
            aN( 0,  4 ) = c*xi*zeta*phi ;
            aN( 0,  5 ) = b*xi*zeta*(4.*xi-1.)*(4.*zeta-1.) ;
            aN( 0,  6 ) = c*xi*zeta*psi ;
            aN( 0,  7 ) = c*eta*zeta*psi ;
            aN( 0,  8 ) = b*eta*zeta*(4.*eta-1.)*(4.*zeta-1.) ;
			aN( 0,  9 ) = c*eta*zeta*delta ;
            aN( 0, 10 ) = c*eta*xi*delta ;
            aN( 0, 11 ) = b*eta*xi*(4.*eta-1.)*(4.*xi-1.) ;
            aN( 0, 12 ) = c*eta*xi*(1.+xi*(8.*xi-6.)) ;
            aN( 0, 13 ) = c*xi*tau*gamma*alpha ;
            aN( 0, 14 ) = b*xi*(4.*xi-1.)*tau*alpha ;
            aN( 0, 15 ) = c*xi*tau*phi ;
            aN( 0, 16 ) = c*eta*tau*gamma*alpha ;
            aN( 0, 17 ) = b*eta*tau*(4.*eta-1.)*alpha ;
            aN( 0, 18 ) = c*eta*tau*delta ;
            aN( 0, 19 ) = c*zeta*tau*gamma*alpha ;
            aN( 0, 20 ) = b*zeta*tau*(4.*zeta-1.)*alpha ;
            aN( 0, 21 ) = c*zeta*tau*psi ;
            aN( 0, 22 ) = d*eta*xi*zeta*(4.*xi-1.) ;
            aN( 0, 23 ) = d*eta*xi*zeta*(4.*eta-1.) ;
            aN( 0, 24 ) = d*eta*xi*zeta*(4.*zeta-1.) ;
            aN( 0, 25 ) = d*xi*zeta*tau*(4.*xi-1.) ;
            aN( 0, 26 ) = d*xi*zeta*(4.*zeta-1.)*tau ;
            aN( 0, 27 ) = d*xi*zeta*tau*alpha ;
            aN( 0, 28 ) = d*eta*xi*tau*(4.*xi-1.) ;
            aN( 0, 29 ) = d*eta*xi*tau*alpha ;
            aN( 0, 30 ) = d*eta*xi*tau*(4.*eta-1.) ;
            aN( 0, 31 ) = d*eta*zeta*tau*alpha ;
            aN( 0, 33 ) = d*eta*zeta*tau*(4.*eta-1.) ;
            aN( 0, 32 ) = d*eta*zeta*tau*(4.*zeta-1.) ;
            aN( 0, 34 ) = 8.*d*xi*eta*zeta*tau ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real a = 1./3.;
            const real b = 4.0;
            const real c = 16.0/3.0;
            const real d = 32.0;

            adNdXi.set_size( 3, 35 );

            adNdXi( 0,  0 ) = a*(8.*xi-3.)*(1.+4.*xi*(4.*xi-3.)) ;
            adNdXi( 1,  0 ) = 0.0 ;
            adNdXi( 2,  0 ) = 0.0 ;

			adNdXi( 0,  1 ) = 0.0 ;
            adNdXi( 1,  1 ) = 0.0 ;
            adNdXi( 2,  1 ) = a*(8.*zeta-3.)*(1.+4.*zeta*(4.*zeta-3.)) ;
            
            adNdXi( 0,  2 ) = 0.0 ;
            adNdXi( 1,  2 ) = a*(8.*eta-3.)*(1.+4.*eta*(4*eta-3.)) ;
            adNdXi( 2,  2 ) = 0.0 ;

            adNdXi( 0,  3 ) = a*(8.*(xi+eta+zeta)-5.)*(5.-20.*(xi+zeta)+4.*(eta*(4.*eta+(8.*(xi+zeta)-5.))+4.*(xi+zeta)*(xi+zeta))) ;
            adNdXi( 1,  3 ) = a*(8.*(xi+eta+zeta)-5.)*(5.-20.*(xi+zeta)+4.*(eta*(4.*eta+8.*(xi+zeta)-5.)+4.*(xi+zeta)*(xi+zeta))) ;
            adNdXi( 2,  3 ) = a*(8.*(xi+eta+zeta)-5.)*((5.-20.*(xi+zeta))+eta*(16.*eta+32.*(xi+zeta)-20.)+16.*(xi+zeta)*(xi+zeta)) ;

            adNdXi( 0,  4 ) = c*zeta*(12.*xi*(2.*xi-1.)+1.) ;
            adNdXi( 1,  4 ) = 0.0 ;
            adNdXi( 2 , 4 ) = c*xi*(xi*(8.*xi-6.)+1.) ;

			adNdXi( 0,  5 ) = b*zeta*(8.*xi-1.)*(4.*zeta-1.) ;
            adNdXi( 1,  5 ) = 0.0 ;
            adNdXi( 2,  5 ) = b*xi*(4.*xi-1.)*(8.*zeta-1.) ;

			adNdXi( 0,  6 ) = c*zeta*(zeta*(8.*zeta-6.)+1.) ;
            adNdXi( 1,  6 ) = 0.0 ;
            adNdXi( 2,  6 ) = c*xi*(12.*zeta*(2.*zeta-1.)+1.) ;

			adNdXi( 0,  7 ) = 0.0 ;
            adNdXi( 1,  7 ) = c*zeta*(zeta*(8.*zeta-6.)+1.) ;
            adNdXi( 2,  7 ) = c*eta*(12.*zeta*(2.*zeta-1.)+1.) ;
            
            adNdXi( 0,  8 ) = 0.0 ;
            adNdXi( 1,  8 ) = b*zeta*(8.*eta-1.)*(4.*zeta-1.) ;
            adNdXi( 2,  8 ) = b*eta*(4.*eta-1.)*(8.*zeta-1.) ;

            adNdXi( 0,  9 ) = 0.0 ;
            adNdXi( 1,  9 ) = c*zeta*(12.*eta*(2.*eta-1.)+1.) ;
            adNdXi( 2,  9 ) = c*eta*(eta*(8.*eta-6.)+1.) ;

			adNdXi( 0, 10 ) = c*eta*(eta*(8.*eta-6.)+1.) ;
            adNdXi( 1, 10 ) = c*xi*(12.*eta*(2.*eta-1.)+1.) ;
            adNdXi( 2, 10 ) = 0.0 ;

			adNdXi( 0, 11 ) = b*eta*(4.*eta-1.)*(8.*xi-1.) ;
            adNdXi( 1, 11 ) = b*xi*(8.*eta-1.)*(4.*xi-1.) ;
            adNdXi( 2, 11 ) = 0.0 ;

            adNdXi( 0, 12 ) = c*eta*(12.*xi*(2.*xi-1.)+1.) ;
            adNdXi( 1, 12 ) = c*xi*(xi*(8.*xi-6.)+1.) ;
            adNdXi( 2, 12 ) = 0.0 ;
            
            adNdXi( 0, 13 ) = c*(3.-26.*xi-13.*zeta-2.*(xi+zeta)*(xi*(16.*xi+20.*zeta-27.)+zeta*(4.*zeta-9.))-eta*((13.-36.*zeta+24.*(3.*(xi-1.)*xi+zeta*(4.*xi+zeta)))+eta*(48.*xi+24.*zeta-18.+8.*eta))) ;
            adNdXi( 1, 13 ) = c*xi*(36.*(xi+zeta)-12.*(eta*(2.*eta+4.*(xi+zeta)-3.)+2.*(xi+zeta)*(xi+zeta))-13.) ;
            adNdXi( 2, 13 ) = c*xi*(36*(xi+zeta)-12.*(eta*(2.*eta+4.*(xi+zeta)-3.)+2*(xi+zeta)*(xi+zeta))-13.) ;

            adNdXi( 0, 14 ) = b*(eta+2.*xi+zeta-1.)*(3.+4.*eta*(8.*xi-1.)-4.*zeta+32.*xi*(xi+zeta-1.)) ;
            adNdXi( 1, 14 ) = b*xi*(4.*xi-1.)*(8.*(xi+eta+zeta)-7.) ;
            adNdXi( 2, 14 ) = b*xi*(4.*xi-1.)*(8.*(xi+eta+zeta)-7.) ;

            adNdXi( 0, 15 ) = c*(1.+eta*(12.*xi*(1.-2.*xi)-1.)-zeta-2.*xi*(7.-6.*zeta+xi*(16.*xi+12.*zeta-21.))) ;
            adNdXi( 1, 15 ) = c*xi*(xi*(6.-8.*xi)-1.) ;
            adNdXi( 2, 15 ) = c*xi*(xi*(6.-8.*xi)-1.) ;

			adNdXi( 0, 16 ) = c*(eta*(36.*(xi+zeta)-12.*(eta*(2.*eta+4.*(xi+zeta)-3.)+2*(xi+zeta)*(xi+zeta))-13.)) ;
            adNdXi( 1, 16 ) = c*((1.0-xi-zeta)*(2.*(xi+zeta)-1.)*(4.*(xi+zeta)-3.)-eta*(eta*(32.*eta+(72.*(xi+zeta)-54.))+(26.+xi*(48.*xi+96.*zeta-72.)+zeta*(48.*zeta-72.)))) ;
            adNdXi( 2, 16 ) = c*eta*(36.*(xi+zeta)-eta*(24.*eta+48.*(xi+zeta)-36.)-24.*(xi+zeta)*(xi+zeta)-13.) ;

            adNdXi( 0, 17 ) = b*eta*(4.*eta-1.)*(8.*(xi+eta+zeta)-7.) ;
            adNdXi( 1, 17 ) = b*(2.*eta+xi+zeta-1.)*(3.-4.*(xi+zeta)+32.*eta*(xi+eta+zeta-1.)) ;
            adNdXi( 2, 17 ) = b*eta*(4.*eta-1.)*(8.*(xi+eta+zeta)-7.) ;

            adNdXi( 0, 18 ) = c*eta*(eta*(6.-8*eta)-1.) ;
            adNdXi( 1, 18 ) = c*(1.0-xi-zeta-2.*eta*(7.+eta*(16*eta+12.*(xi+zeta)-21.)-6.*(xi+zeta))) ;
            adNdXi( 2, 18 ) = c*eta*(eta*(6.-8.*eta)-1.) ;
            
            adNdXi( 0, 19 ) = c*zeta*(36.*(xi+zeta)-12.*(eta*(2.*eta+4.*(xi+zeta)-3.)+2.*(xi+zeta)*(xi+zeta))-13.) ;
            adNdXi( 1, 19 ) = c*zeta*(36.*(xi+zeta)-12*(eta*(2.*eta+4*(xi+zeta)-3.)+2*(xi+zeta)*(xi+zeta))-13.) ;
            adNdXi( 2, 19 ) = c*(3.0-13.*xi-26*zeta-2.*(xi+zeta)*(xi*(4.*xi+20.*zeta-9.)+zeta*(16.*zeta-27.))-eta*(eta*(8.*eta+24.*xi+48.*zeta-18.)+(13.+xi*(24.*xi+12.*(8.*zeta-3.))+72.*(zeta-1.)*zeta))) ;

            adNdXi( 0, 20 ) = b*zeta*(4.*zeta-1.)*(8.*(xi+eta+zeta)-7.) ;
            adNdXi( 1, 20 ) = b*zeta*(4.*zeta-1.)*(8.*(xi+eta+zeta)-7.) ;
            adNdXi( 2, 20 ) = b*(xi+eta+2.*zeta-1.)*(3.-4.*(xi+eta)+32.*zeta*(xi+eta+zeta-1.)) ;

            adNdXi( 0, 21 ) = c*zeta*(zeta*(6.-8.*zeta)-1.) ;
            adNdXi( 1, 21 ) = c*zeta*(zeta*(6.-8.*zeta)-1.) ;
            adNdXi( 2, 21 ) = c*(1-eta-xi+zeta*(2.*(6.*(xi+eta)-7.)-zeta*(24.*(xi+eta)+32.*zeta-42.))) ;

            adNdXi( 0, 22 ) = d*eta*zeta*(8.*xi-1.) ;
            adNdXi( 1, 22 ) = d*xi*zeta*(4.*xi-1.) ;
            adNdXi( 2, 22 ) = d*eta*xi*(4.*xi-1.) ;
            
            adNdXi( 0, 23 ) = d*eta*zeta*(4.*eta-1.) ;
            adNdXi( 1, 23 ) = d*xi*zeta*(8.*eta-1.) ;
            adNdXi( 2, 23 ) = d*eta*xi*(4.*eta-1.) ;

			adNdXi( 0, 24 ) = d*eta*zeta*(4.*zeta-1.) ;
            adNdXi( 1, 24 ) = d*xi*zeta*(4.*zeta-1.) ;
            adNdXi( 2, 24 ) = d*eta*xi*(8.*zeta-1.) ;
            
            adNdXi( 0, 25 ) = d*zeta*(eta+zeta-xi*(8.*eta+12.*xi+8.*zeta-10.)-1.) ;
            adNdXi( 1, 25 ) = d*xi*zeta*(1.-4.*xi) ;
            adNdXi( 2, 25 ) = d*xi*(1.-4.*xi)*(xi+eta+2.*zeta-1.) ;
            
            adNdXi( 0, 26 ) = d*zeta*(1.-4.*zeta)*(eta+2.*xi+zeta-1.) ;
            adNdXi( 1, 26 ) = d*xi*zeta*(1.-4.*zeta) ;
            adNdXi( 2, 26 ) = d*xi*(xi+eta+zeta*(10.-8.*(eta+xi)-12.*zeta)-1.) ;

            adNdXi( 0, 27 ) = d*zeta*(3+eta*(4*eta+16.*xi+8.*zeta-7.)-14.*xi-7.*zeta+4.*(xi+zeta)*(3.*xi+zeta)) ;
            adNdXi( 1, 27 ) = d*xi*zeta*(8.*(xi+eta+zeta)-7.) ;
            adNdXi( 2, 27 ) = d*xi*(3.+eta*(4.*eta+8.*xi+16.*zeta-7.)-7.*xi-14.*zeta+4.*(xi+zeta)*(xi+3.*zeta)) ;

            adNdXi( 0, 28 ) = d*eta*(eta-8.*eta*xi+zeta-2.*xi*(6.*xi+4.*zeta-5.)-1.) ;
            adNdXi( 1, 28 ) = d*xi*(1.-4*xi)*(xi+2.*eta+zeta-1.) ;
            adNdXi( 2, 28 ) = d*eta*xi*(1.-4*xi) ;

			adNdXi( 0, 29 ) = d*eta*(3.+eta*(4.*eta+16.*xi+8.*zeta-7.)-14.*xi-7.*zeta+4.*(xi+zeta)*(3.*xi+zeta)) ;
            adNdXi( 1, 29 ) = d*xi*(eta*(12.*eta+16.*(xi+zeta)-14.)+(xi+zeta-1.)*(4.*(xi+zeta)-3.)) ;
            adNdXi( 2, 29 ) = d*eta*xi*(8.*(xi+eta+zeta)-7.) ;
            
            adNdXi( 0, 30 ) = d*eta*(4.*eta-1.)*(1.-eta-2.*xi-zeta) ;
            adNdXi( 1, 30 ) = d*xi*(xi+zeta-2.*eta*(6.*eta+4.*(xi+zeta)-5.)-1) ;
            adNdXi( 2, 30 ) = d*eta*xi*(1.-4.*eta) ;

            adNdXi( 0, 31 ) = d*eta*zeta*(8.*(xi+eta+zeta)-7.) ;
            adNdXi( 1, 31 ) = d*zeta*(eta*(12.*eta+16.*(xi+zeta)-14.)+(xi+zeta-1.)*(4.*(xi+zeta)-3.)) ;
            adNdXi( 2, 31 ) = d*eta*(3.+eta*(4.*eta+8.*xi+16.*zeta-7.)-7.*xi-14.*zeta+4.*(xi+zeta)*(xi+3.*zeta)) ;

			adNdXi( 0, 32 ) = d*eta*zeta*(1.-4.*zeta) ;
            adNdXi( 1, 32 ) = d*zeta*(1.-4.*zeta)*(xi+2.*eta+zeta-1.) ;
            adNdXi( 2, 32 ) = d*eta*(xi+eta+zeta*(10.-8.*(eta+xi)-12.*zeta)-1.) ;


            adNdXi( 0, 33 ) = d*eta*zeta*(1.-4.*eta) ;
            adNdXi( 1, 33 ) = d*zeta*(xi+zeta-2.*eta*(6.*eta+4.*(xi+zeta)-5.)-1.) ;
            adNdXi( 2, 33 ) = d*eta*(1.-4.*eta)*(xi+eta+2.*zeta-1.) ;

            adNdXi( 0, 34 ) = 8.*d*eta*zeta*(1.-2.*xi-eta-zeta) ;
            adNdXi( 1, 34 ) = 8.*d*xi*zeta*(1.-xi-2.*eta-zeta) ;
            adNdXi( 2, 34 ) = 8.*d*eta*xi*(1.-xi-eta-2.*zeta) ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            const real a = 1./3.;
            const real b = 4.0;
            const real c = 16.0/3.0;
            const real d = 32.0;

            ad2NdXi2.set_size( 6, 35, 0.0 );

            ad2NdXi2( 0,  0 ) = 4.*a*(xi*(96.*xi-72.)+11.) ;
            ad2NdXi2( 1,  0 ) = 0.0 ;
            ad2NdXi2( 2,  0 ) = 0.0 ;
            ad2NdXi2( 3,  0 ) = 0.0 ;
            ad2NdXi2( 4,  0 ) = 0.0 ;
            ad2NdXi2( 5,  0 ) = 0.0 ;

			ad2NdXi2( 0,  1 ) = 0.0 ;
            ad2NdXi2( 1,  1 ) = 0.0 ;
            ad2NdXi2( 2,  1 ) = 4.*a*(zeta*(96.*zeta-72.)+11.) ;
            ad2NdXi2( 3,  1 ) = 0.0 ;
            ad2NdXi2( 4,  1 ) = 0.0 ;
            ad2NdXi2( 5,  1 ) = 0.0 ;
            
            ad2NdXi2( 0,  2 ) = 0.0 ;
            ad2NdXi2( 1,  2 ) = 4.*a*(eta*(96.*eta-72.)+11.) ;
            ad2NdXi2( 2,  2 ) = 0.0 ;
            ad2NdXi2( 3,  2 ) = 0.0 ;
            ad2NdXi2( 4,  2 ) = 0.0 ;
            ad2NdXi2( 5,  2 ) = 0.0 ;

            ad2NdXi2( 0,  3 ) = 4.*a*(35.-120.*(xi+zeta)+24.*(eta*(4.*eta+8.*(xi+zeta)-5.)+4.*(xi+zeta)*(xi+zeta))) ;
            ad2NdXi2( 1,  3 ) = 4.*a*(35.-120.*(xi+zeta)+eta*(96.*eta+192.*(xi+zeta)-120.)+96.*(xi+zeta)*(xi+zeta)) ;
            ad2NdXi2( 2,  3 ) = 4.*a*((35.-120.*(xi+zeta))+24.*eta*(4.*eta+8.*(xi+zeta)-5.)+96.*(xi+zeta)*(xi+zeta)) ;
            ad2NdXi2( 3,  3 ) = 4.*a*(35.-120.*(xi+zeta)+eta*(96.*eta+192.*(xi+zeta)-120.)+96.*(xi+zeta)*(xi+zeta)) ;
            ad2NdXi2( 4,  3 ) = 4.*a*(35.-120.*(xi+zeta)+24.*eta*(4*eta+8.*(xi+zeta)-5.)+96.*(xi+zeta)*(xi+zeta)) ;
            ad2NdXi2( 5,  3 ) = 4.*a*(35.-120.*(xi+zeta)+eta*(96.*eta+192.*(xi+zeta)-120.)+96*(xi+zeta)*(xi+zeta)) ;
            
            ad2NdXi2( 0,  4 ) = c*zeta*(48.*xi-12.) ;
            ad2NdXi2( 1,  4 ) = 0.0 ;
            ad2NdXi2( 2,  4 ) = 0.0 ;
            ad2NdXi2( 3,  4 ) = 0.0 ;
            ad2NdXi2( 4,  4 ) = c*(12.*xi*(2.*xi-1.)+1.) ;
            ad2NdXi2( 5,  4 ) = 0.0 ;

			ad2NdXi2( 0,  5 ) = 8.*b*zeta*(4.*zeta-1.) ;
            ad2NdXi2( 1,  5 ) = 0.0 ;
            ad2NdXi2( 2,  5 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 3,  5 ) = 0.0 ;
            ad2NdXi2( 4,  5 ) = b*(8.*xi-1.)*(8.*zeta-1.) ;
            ad2NdXi2( 5,  5 ) = 0.0 ;
            
            ad2NdXi2( 0,  6 ) = 0.0 ;
            ad2NdXi2( 1,  6 ) = 0.0 ;
            ad2NdXi2( 2,  6 ) = c*xi*(48.*zeta-12.) ;
            ad2NdXi2( 3,  6 ) = 0.0 ;
            ad2NdXi2( 4,  6 ) = c*(12.*zeta*(2.*zeta-1.)+1.) ;
            ad2NdXi2( 5,  6 ) = 0.0 ;
            
            ad2NdXi2( 0,  7 ) = 0.0 ;
            ad2NdXi2( 1,  7 ) = 0.0 ;
            ad2NdXi2( 2,  7 ) = c*eta*(48.*zeta-12.) ;
            ad2NdXi2( 3,  7 ) = c*(12.*zeta*(2.*zeta-1.)+1.) ;
            ad2NdXi2( 4,  7 ) = 0.0 ;
            ad2NdXi2( 5,  7 ) = 0.0 ;
            
            ad2NdXi2( 0,  8 ) = 0.0 ;
            ad2NdXi2( 1,  8 ) = 8.*b*zeta*(4.*zeta-1.) ;
            ad2NdXi2( 2,  8 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 3,  8 ) = b*(8.*eta-1.)*(8.*zeta-1.) ;
            ad2NdXi2( 4,  8 ) = 0.0 ;
            ad2NdXi2( 5,  8 ) = 0.0 ;
            
            ad2NdXi2( 0,  9 ) = 0.0 ;
            ad2NdXi2( 1,  9 ) = c*zeta*(48.*eta-12.) ;
            ad2NdXi2( 2,  9 ) = 0.0 ;
            ad2NdXi2( 3,  9 ) = c*(12.*eta*(2.*eta-1.)+1.) ;
            ad2NdXi2( 4,  9 ) = 0.0 ;
            ad2NdXi2( 5,  9 ) = 0.0 ;

			ad2NdXi2( 0, 10 ) = 0.0 ;
            ad2NdXi2( 1, 10 ) = c*xi*(48.*eta-12.) ;
            ad2NdXi2( 2, 10 ) = 0.0 ;
            ad2NdXi2( 3, 10 ) = 0.0 ;
            ad2NdXi2( 4, 10 ) = 0.0 ;
            ad2NdXi2( 5, 10 ) = c*(12.*eta*(2.*eta-1.)+1.) ;
            
            ad2NdXi2( 0, 11 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 1, 11 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 2, 11 ) = 0.0 ;
            ad2NdXi2( 3, 11 ) = 0.0 ;
            ad2NdXi2( 4, 11 ) = 0.0 ;
            ad2NdXi2( 5, 11 ) = b*(8.*eta-1.)*(8.*xi-1.) ;
            
            ad2NdXi2( 0, 12 ) = c*eta*(48.*xi-12.) ;
            ad2NdXi2( 1, 12 ) = 0.0 ;
            ad2NdXi2( 2, 12 ) = 0.0 ;
            ad2NdXi2( 3, 12 ) = 0.0 ;
            ad2NdXi2( 4, 12 ) = 0.0 ;
            ad2NdXi2( 5, 12 ) = c*(12.*xi*(2.*xi-1.)+1.) ;
            
            ad2NdXi2( 0, 13 ) = c*(108.*xi+72.*zeta-24.*eta*(2.*eta+6.*xi+4.*zeta-3.)-48.*(xi+zeta)*(2.*xi+zeta)-26.) ;
            ad2NdXi2( 1, 13 ) = c*xi*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 2, 13 ) = c*xi*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 3, 13 ) = c*xi*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 4, 13 ) = c*(72.*xi+36.*zeta-12.*eta*(2.*eta+8.*xi+4.*zeta-3.)-24.*(xi+zeta)*(3.*xi+zeta)-13.) ;
            ad2NdXi2( 5, 13 ) = c*(72.*xi+36.*zeta-12.*eta*(2*eta+8.*xi+4.*zeta-3.)-24.*(xi+zeta)*(3.*xi+zeta)-13.) ;

            ad2NdXi2( 0, 14 ) = 2.*b*(19.+eta*(16.*eta+96.*xi+32.*zeta-36.)+96.*xi*(xi+zeta-1.)+4*zeta*(4.*zeta-9.)) ;
            ad2NdXi2( 1, 14 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 2, 14 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 3, 14 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 4, 14 ) = b*(7.+8.*(eta*(8.*xi-1.)-zeta+xi*(12.*xi+8.*zeta-9.))) ;
            ad2NdXi2( 5, 14 ) = b*(7.+8.*(eta*(8.*xi-1.)-zeta+xi*(12.*xi+8.*zeta-9.))) ;

            ad2NdXi2( 0, 15 ) = c*(12.*(zeta-eta*(4.*xi-1.)-xi*(8.*xi+4.*zeta-7.))-14.) ;
            ad2NdXi2( 1, 15 ) = 0.0 ;
            ad2NdXi2( 2, 15 ) = 0.0 ;
            ad2NdXi2( 3, 15 ) = 0.0 ;
            ad2NdXi2( 4, 15 ) = c*(12.*xi*(1.-2.*xi)-1.) ;
            ad2NdXi2( 5, 15 ) = c*(12.*xi*(1.-2.*xi)-1.) ;

			ad2NdXi2( 0, 16 ) = c*eta*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 1, 16 ) = c*(72*(xi+zeta)-12.*(8.*eta*eta+4.*(xi+zeta)*(xi+zeta)+3.*eta*(4.*(xi+zeta)-3.))-26.) ;
            ad2NdXi2( 2, 16 ) = c*eta*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 3, 16 ) = c*(36.*(xi+zeta)-24.*(eta*(3.*eta+4.*(xi+zeta)-3.)+(xi+zeta)*(xi+zeta))-13.) ;
            ad2NdXi2( 4, 16 ) = c*eta*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 5, 16 ) = c*(36.*(xi+zeta)-eta*(72*eta+96.*(xi+zeta)-72.)-24.*(xi+zeta)*(xi+zeta)-13.) ;

            ad2NdXi2( 0, 17 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 1, 17 ) = b*(6.-8.*(xi+zeta)+64*eta*(eta+xi+zeta-1.)+32.*(xi+2.*eta+zeta-1.)*(xi+2.*eta+zeta-1.)) ;
            ad2NdXi2( 2, 17 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 3, 17 ) = b*(7.-8.*(xi+zeta)+8.*eta*(8.*(xi+zeta)+12*eta-9.)) ;
            ad2NdXi2( 4, 17 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 5, 17 ) = b*(7.+8.*(eta*(12.*eta+8.*(xi+zeta)-9.)-xi-zeta)) ;

            ad2NdXi2( 0, 18 ) = 0.0 ;
            ad2NdXi2( 1, 18 ) = c*(12.*(xi+zeta-eta*(8.*eta+4.*(xi+zeta)-7.))-14.) ;
            ad2NdXi2( 2, 18 ) = 0.0 ;
            ad2NdXi2( 3, 18 ) = c*(2.*eta*(6.-12.*eta)-1) ;
            ad2NdXi2( 4, 18 ) = 0.0 ;
            ad2NdXi2( 5, 18 ) = c*(12.*eta*(1.-2.*eta)-1.) ;
            
            ad2NdXi2( 0, 19 ) = c*zeta*(36-48*(xi+eta+zeta)) ;
            ad2NdXi2( 1, 19 ) = c*zeta*(36.-48.*(xi+eta+zeta)) ;
            ad2NdXi2( 2, 19 ) = c*(72.*xi+108.*zeta-eta*(48.*eta+96.*xi+144.*zeta-72.)-48.*(xi+zeta)*(xi+2.*zeta)-26.) ;
            ad2NdXi2( 3, 19 ) = c*(36.*xi+72*zeta-eta*(24.*eta+48.*xi+96.*zeta-36.)-24.*(xi+zeta)*(xi+3.*zeta)-13.) ;
            ad2NdXi2( 4, 19 ) = c*(36.*xi+72.*zeta-12*eta*(2*eta+4.*xi+8.*zeta-3.)-24.*(xi+zeta)*(xi+3.*zeta)-13.) ;
            ad2NdXi2( 5, 19 ) = c*zeta*(36.-48*(xi+eta+zeta)) ;

            ad2NdXi2( 0, 20 ) = 8.*b*zeta*(4.*zeta-1.) ;
            ad2NdXi2( 1, 20 ) = 8.*b*zeta*(4.*zeta-1.) ;
            ad2NdXi2( 2, 20 ) = b*(6.-8.*(xi+eta)+64*zeta*(eta+xi+zeta-1.)+32.*(eta+xi+2.*zeta-1.)*(eta+xi+2.*zeta-1.)) ;
            ad2NdXi2( 3, 20 ) = b*(7.-8.*(xi+eta)+zeta*(64.*(xi+eta)+96*zeta-72.)) ;
            ad2NdXi2( 4, 20 ) = b*(8.*(zeta*(8.*(xi+eta)+12.*zeta-9.)-xi-eta)+7.) ;
            ad2NdXi2( 5, 20 ) = 8.*b*zeta*(4.*zeta-1.) ;

            ad2NdXi2( 0, 21 ) = 0.0 ;
            ad2NdXi2( 1, 21 ) = 0.0 ;
            ad2NdXi2( 2, 21 ) = c*(12*(xi+eta-zeta*(4.*(xi+eta)-7.+8.*zeta))-14.) ;
            ad2NdXi2( 3, 21 ) = c*(12.*zeta*(1.-2.*zeta)-1.) ;
            ad2NdXi2( 4, 21 ) = c*(12.*zeta*(1-2.*zeta)-1.) ;
            ad2NdXi2( 5, 21 ) = 0.0 ;

            ad2NdXi2( 0, 22 ) = 8.*d*eta*zeta ;
            ad2NdXi2( 1, 22 ) = 0.0 ;
            ad2NdXi2( 2, 22 ) = 0.0 ;
            ad2NdXi2( 3, 22 ) = d*xi*(4.*xi-1.) ;
            ad2NdXi2( 4, 22 ) = d*eta*(8.*xi-1.) ;
            ad2NdXi2( 5, 22 ) = d*zeta*(8.*xi-1.) ;

			ad2NdXi2( 0, 23 ) = 0.0 ;
            ad2NdXi2( 1, 23 ) = 8.*d*xi*zeta ;
            ad2NdXi2( 2, 23 ) = 0.0 ;
            ad2NdXi2( 3, 23 ) = d*xi*(8.*eta-1.) ;
            ad2NdXi2( 4, 23 ) = d*eta*(4.*eta-1.) ;
            ad2NdXi2( 5, 23 ) = d*zeta*(8.*eta-1.) ;

            ad2NdXi2( 0, 24 ) = 0.0 ;
            ad2NdXi2( 1, 24 ) = 0.0 ;
            ad2NdXi2( 2, 24 ) = 8.*d*xi*eta ;
            ad2NdXi2( 3, 24 ) = d*xi*(8.*zeta-1.) ;
            ad2NdXi2( 4, 24 ) = d*eta*(8.*zeta-1.) ;
            ad2NdXi2( 5, 24 ) = d*zeta*(4.*zeta-1.) ;

            ad2NdXi2( 0, 25 ) = d*zeta*(10-24.*xi-8.*(eta+zeta)) ;
            ad2NdXi2( 1, 25 ) = 0.0 ;
            ad2NdXi2( 2, 25 ) = d*xi*(2.-8*xi) ;
            ad2NdXi2( 3, 25 ) = d*xi*(1.-4.*xi) ;
            ad2NdXi2( 4, 25 ) = d*(eta-8.*eta*xi+2.*zeta-2.*xi*(6.*xi+8.*zeta-5.)-1.) ;
            ad2NdXi2( 5, 25 ) = d*zeta*(1.-8.*xi) ;

			ad2NdXi2( 0, 26 ) = d*zeta*(2.-8*zeta) ;
            ad2NdXi2( 1, 26 ) = 0.0 ;
            ad2NdXi2( 2, 26 ) = d*xi*(10-8.*(xi+eta)-24.*zeta) ;
            ad2NdXi2( 3, 26 ) = d*xi*(1.-8.*zeta) ;
            ad2NdXi2( 4, 26 ) = d*(eta*(1.-8.*zeta)+xi*(2.-16.*zeta)+zeta*(10.-12.*zeta)-1.) ;
            ad2NdXi2( 5, 26 ) = d*zeta*(1.-4.*zeta) ;
            
            ad2NdXi2( 0, 27 ) = d*zeta*(24.*xi+16*(eta+zeta)-14.) ;
            ad2NdXi2( 1, 27 ) = 8.*d*xi*zeta ;
            ad2NdXi2( 2, 27 ) = d*xi*(16.*(xi+eta)+24.*zeta-14.) ;
            ad2NdXi2( 3, 27 ) = d*xi*(8*(xi+eta)+16.*zeta-7.) ;
            ad2NdXi2( 4, 27 ) = d*(3.+eta*(4.*eta+16*(xi+zeta)-7.)+xi*(12.*xi+32.*zeta-14.)+2.*zeta*(6.*zeta-7.)) ;
            ad2NdXi2( 5, 27 ) = d*zeta*(16.*xi+8*(eta+zeta)-7.) ;
            
            ad2NdXi2( 0, 28 ) = d*eta*(10.-24.*xi-8.*(eta+zeta)) ;
            ad2NdXi2( 1, 28 ) = d*xi*(2.-8.*xi) ;
            ad2NdXi2( 2, 28 ) = 0.0 ;
            ad2NdXi2( 3, 28 ) = d*xi*(1.-4.*xi) ;
            ad2NdXi2( 4, 28 ) = d*eta*(1.-8.*xi) ;
            ad2NdXi2( 5, 28 ) = d*(eta*(2.-16.*xi)+zeta-2.*xi*(-5+6*xi+4*zeta)-1.) ;

			ad2NdXi2( 0, 29 ) = d*eta*(24.*xi+16*(eta+zeta)-14.) ;
            ad2NdXi2( 1, 29 ) = d*xi*(24.*eta+16.*(xi+zeta)-14.) ;
            ad2NdXi2( 2, 29 ) = 8.*d*xi*eta ;
            ad2NdXi2( 3, 29 ) = d*xi*(16.*eta+8.*(xi+zeta)-7.) ;
            ad2NdXi2( 4, 29 ) = d*eta*(16.*xi+8.*(zeta+eta)-7.) ;
            ad2NdXi2( 5, 29 ) = d*(3.+eta*(12.*eta+32.*xi+16.*zeta-14.)-14.*xi-7.*zeta+4.*(xi+zeta)*(3.*xi+zeta)) ;
            
            ad2NdXi2( 0, 30 ) = d*eta*(2.-8.*eta) ;
            ad2NdXi2( 1, 30 ) = d*xi*(10.-24.*eta-8.*(xi+zeta)) ;
            ad2NdXi2( 2, 30 ) = 0.0 ;
            ad2NdXi2( 3, 30 ) = d*xi*(1.-8.*eta) ;
            ad2NdXi2( 4, 30 ) = d*eta*(1.-4.*eta) ;
            ad2NdXi2( 5, 30 ) = d*(2.*xi+zeta-2.*eta*(6.*eta+8.*xi+4.*zeta-5.)-1.) ;

            ad2NdXi2( 0, 31 ) = 8.*d*eta*zeta ;
            ad2NdXi2( 1, 31 ) = d*zeta*(24.*eta+16.*(xi+zeta)-14.) ;
            ad2NdXi2( 2, 31 ) = d*eta*(16.*(xi+eta)+24.*zeta-14.) ;
            ad2NdXi2( 3, 31 ) = d*(3.+2.*eta*(6.*eta+8.*xi+16.*zeta-7.)-7.*xi-14.*zeta+4.*(xi+zeta)*(xi+3.*zeta)) ;
            ad2NdXi2( 4, 31 ) = d*eta*(8.*(xi+eta)+16.*zeta-7.) ;
            ad2NdXi2( 5, 31 ) = d*zeta*(16*eta+8*(xi+zeta)-7.) ;

			ad2NdXi2( 0, 32 ) = 0.0 ;
            ad2NdXi2( 1, 32 ) = d*zeta*(2.-8.*zeta) ;
            ad2NdXi2( 2, 32 ) = d*eta*(10-8.*(xi+eta)-24.*zeta) ;
            ad2NdXi2( 3, 32 ) = d*(xi+eta*(2.-16.*zeta)-8.*xi*zeta+(10.-12.*zeta)*zeta-1.) ;
            ad2NdXi2( 4, 32 ) = d*eta*(1.-8.*zeta) ;
            ad2NdXi2( 5, 32 ) = d*zeta*(1.-4.*zeta) ;

            ad2NdXi2( 0, 33 ) = 0.0 ;
            ad2NdXi2( 1, 33 ) = d*zeta*(10.-24.*eta-8.*(xi+zeta)) ;
            ad2NdXi2( 2, 33 ) = d*eta*(2.-8*eta) ;
            ad2NdXi2( 3, 33 ) = d*(xi+2.*zeta-2.*eta*(6.*eta+4.*xi+8.*zeta-5.)-1.) ;
            ad2NdXi2( 4, 33 ) = d*eta*(1.-4.*eta) ;
            ad2NdXi2( 5, 33 ) = d*zeta*(1.-8.*eta) ;

            
            ad2NdXi2( 0, 34 ) = -16*d*eta*zeta ;
            ad2NdXi2( 1, 34 ) = -16*d*xi*zeta ;
            ad2NdXi2( 2, 34 ) = -16.*d*xi*eta ;
            ad2NdXi2( 3, 34 ) = 8.*d*xi*(1.-xi-2.*(eta+zeta)) ;
            ad2NdXi2( 4, 34 ) = 8.*d*eta*(1.-2.*(xi+zeta)-eta) ;
            ad2NdXi2( 5, 34 ) = 8*d*zeta*(1.-2.*(xi+eta)-zeta) ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TET35_HPP
