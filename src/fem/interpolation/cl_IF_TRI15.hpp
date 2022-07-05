//
// Created by Christian Messe on 25.07.20.
//

#ifndef BELFEM_CL_IF_TRI15_HPP
#define BELFEM_CL_IF_TRI15_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUARTIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >
        ::element_type() const
        {
            return ElementType::TRI15;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 15 );

            aXiHat( 0,  0 ) = 1.00 ;
            aXiHat( 0,  1 ) = 0.00 ;
            aXiHat( 0,  2 ) = 0.00 ;
            aXiHat( 0,  3 ) = 0.75 ;
            aXiHat( 0,  4 ) = 0.50 ;
            aXiHat( 0,  5 ) = 0.25 ;
            aXiHat( 0,  6 ) = 0.00 ;
            aXiHat( 0,  7 ) = 0.00 ;
            aXiHat( 0,  8 ) = 0.00 ;
            aXiHat( 0,  9 ) = 0.25 ;
            aXiHat( 0, 10 ) = 0.50 ;
            aXiHat( 0, 11 ) = 0.75 ;
            aXiHat( 0, 12 ) = 0.50 ;
            aXiHat( 0, 13 ) = 0.25 ;
            aXiHat( 0, 14 ) = 0.25 ;

            aXiHat( 1,  0 ) = 0.00 ;
            aXiHat( 1,  1 ) = 1.00 ;
            aXiHat( 1,  2 ) = 0.00 ;
            aXiHat( 1,  3 ) = 0.25 ;
            aXiHat( 1,  4 ) = 0.50 ;
            aXiHat( 1,  5 ) = 0.75 ;
            aXiHat( 1,  6 ) = 0.75 ;
            aXiHat( 1,  7 ) = 0.50 ;
            aXiHat( 1,  8 ) = 0.25 ;
            aXiHat( 1,  9 ) = 0.00 ;
            aXiHat( 1, 10 ) = 0.00 ;
            aXiHat( 1, 11 ) = 0.00 ;
            aXiHat( 1, 12 ) = 0.25 ;
            aXiHat( 1, 13 ) = 0.50 ;
            aXiHat( 1, 14 ) = 0.25 ;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real a = 1.0/3.0 ;
            const real b = 4.0 ;
            const real c = 16.0 / 3.0 ;
            const real d = 32.0 ;

            aN.set_size( 1, 15 );

            aN( 0,  0 ) =  a*xi*(2.*xi-1.)*(4.*xi-3.)*(4.*xi-1.);
            aN( 0,  1 ) =  a*eta*(2.*eta-1.)*(4.*eta-3.)*(4.*eta-1.);
            aN( 0,  2 ) =  a*(xi+eta-1.)*(2.*(xi+eta)-1.)*(4.*(xi+eta)-3.)*(4.*(xi+eta)-1.);
            aN( 0,  3 ) =  c*eta*xi*(8.*xi*xi-6.*xi+1.);
            aN( 0,  4 ) =  b*eta*xi*(4.*eta-1.)*(4.*xi-1.);
            aN( 0,  5 ) =  c*eta*xi*(8.*eta*eta-6.*eta+1.);
            aN( 0,  6 ) = -c*eta*(xi+eta-1.)*(8.*eta*eta-6.*eta+1.);
            aN( 0,  7 ) =  b*eta*(4.*eta-1.)*(xi+eta-1.)*(4.*(xi+eta)-3.);
            aN( 0,  8 ) = -c*eta*(xi+eta-1.)*(2.*(xi+eta)-1.)*(4.*(xi+eta)-3.);
            aN( 0,  9 ) = -c*xi*(xi+eta-1.)*(2.*(xi+eta)-1.)*(4.*(xi+eta)-3.);
            aN( 0, 10 ) =  b*xi*(xi+eta-1.)*(4.*xi-1.)*(4.*(xi+eta)-3.);
            aN( 0, 11 ) = -c*xi*(xi+eta-1.)*(8.*xi*xi-6.*xi+1.);
            aN( 0, 12 ) = -d*eta*xi*(4.*xi-1.)*(xi+eta-1.);
            aN( 0, 13 ) = -d*eta*xi*(4.*eta-1.)*(xi+eta-1.);
            aN( 0, 14 ) =  d*eta*xi*(xi+eta-1.)*(4.*(xi+eta)-3.);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            adNdXi.set_size( 2, 15 );

            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real a = 1.0/3.0 ;
            const real b = 4.0 ;
            const real c = 16.0 / 3.0 ;
            const real d = 32.0 ;

            adNdXi( 0,  0 ) = a*(8.*xi-3.)*(1.+4.*xi*(4.*xi-3.)) ;
            adNdXi( 1,  0 ) = 0.0 ;

            adNdXi( 0,  1 ) = 0.0 ;
            adNdXi( 1,  1 ) = a*(8.*eta-3.)*(1+4.*eta*(4.*eta-3.)) ;

            adNdXi( 0,  2 ) = a*(8.*(xi+eta)-5.)*(5.+16.*eta*eta+4.*xi*(4.*xi-5.)+4.*eta*(8.*xi-5.)) ;
            adNdXi( 1,  2 ) = a*(8.*(xi+eta)-5.)*(5.+4*xi*(4.*xi-5.)+4*eta*((8.*xi-5.)+4.*eta)) ;

            adNdXi( 0,  3 ) = c*eta*(1.+12.*xi*(2.*xi-1.)) ;
            adNdXi( 1,  3 ) = c*xi*(8.*xi*xi-6.*xi+1.) ;

            adNdXi( 0,  4 ) = b*eta*(4.*eta-1.)*(8.*xi-1.) ;
            adNdXi( 1,  4 ) = b*xi*(8.*eta-1.)*(4.*xi-1.) ;

            adNdXi( 0,  5 ) = c*eta*(eta*(8.*eta-6.)+1.) ;
            adNdXi( 1,  5 ) = c*xi*(24.*eta*eta-12*eta+1.) ;

            adNdXi( 0,  6 ) = -c*eta*(eta*(8.*eta-6.)+1.) ;
            adNdXi( 1,  6 ) = -(c*(xi+2.*eta*(7.-6.*xi+eta*(16*eta+12.*xi-21.))-1.)) ;

            adNdXi( 0,  7 ) = b*eta*(4.*eta-1.)*(8.*(eta+xi)-7.) ;
            adNdXi( 1,  7 ) = b*(2.*eta+xi-1.)*(3.-4.*xi+32*eta*(xi+eta-1.)) ;

            adNdXi( 0,  8 ) = -(c*eta*(13.+24.*eta*eta+12.*xi*(2.*xi-3.)+12.*eta*(4.*xi-3.))) ;
            adNdXi( 1,  8 ) = -(c*(-3.+2*eta*(13.+eta*(16*eta-27.))+xi*(13.+72.*(eta-1.)*eta+xi*(6.*(8.*eta-3.)+8.*xi)))) ;

            adNdXi( 0,  9 ) = -(c*(eta*(eta*(8.*eta+6.*(8.*xi-3.))+(13.+72.*(xi-1.)*xi))+2.*xi*(13.+xi*(16.*xi-27.))-3.)) ;
            adNdXi( 1,  9 ) = -(c*xi*(13+24.*eta*eta+12.*xi*(2*xi-3.)+12.*eta*(4.*xi-3.))) ;

            adNdXi( 0, 10 ) = b*(2.*xi+eta-1.)*(3.+32.*(xi-1.)*xi+4.*eta*(8.*xi-1.)) ;
            adNdXi( 1, 10 ) = b*xi*(4.*xi-1.)*(8.*eta+8.*xi-7.) ;

            adNdXi( 0, 11 ) = -(c*(eta+12.*eta*xi*(2.*xi-1.)+2.*xi*(7.+xi*(16.*xi-21.))-1.)) ;
            adNdXi( 1, 11 ) = -c*xi*(8.*xi*xi-6.*xi+1.) ;

            adNdXi( 0, 12 ) = d*eta*(eta-8.*eta*xi+2.*(5.-6.*xi)*xi-1.) ;
            adNdXi( 1, 12 ) = -d*xi*(4.*xi-1.)*(2*eta+xi-1.) ;

            adNdXi( 0, 13 ) = -(d*eta*(4.*eta-1.)*(2.*xi+eta-1.)) ;
            adNdXi( 1, 13 ) = d*xi*(xi-2*eta*(6*eta+4.*xi-5.)-1.) ;

            adNdXi( 0, 14 ) = d*eta*(3.+4.*eta*eta+2.*xi*(6.*xi-7.)+eta*(16.*xi-7.)) ;
            adNdXi( 1, 14 ) = d*xi*(3.+2*eta*(6*eta-7.)+xi*(4.*xi-7.+16.*eta)) ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 3, 15, 0.0 );

            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real a = 1.0/3.0 ;
            const real b = 4.0 ;
            const real c = 16.0 / 3.0 ;
            const real d = 32.0 ;

            ad2NdXi2( 0,  0 ) = 4.*a*(11.+24.*xi*(4.*xi-3.)) ;
            ad2NdXi2( 1,  0 ) = 0.0 ;
            ad2NdXi2( 2,  0 ) = 0.0  ;

            ad2NdXi2( 0,  1 ) = 0.0 ;
            ad2NdXi2( 1,  1 ) = 4.*a*(11.+24.*eta*(4.*eta-3.)) ;
            ad2NdXi2( 2,  1 ) = 0.0 ;

            ad2NdXi2( 0,  2 ) = 4.*a*(35.+24.*xi*(4.*xi-5.)+eta*(96.*eta+24*(8.*xi-5.))) ;
            ad2NdXi2( 1,  2 ) = 4.*a*(35.+24.*xi*(4.*xi-5.)+eta*(96.*eta+24*(8.*xi-5.))) ;
            ad2NdXi2( 2,  2 ) = 4.*a*(35.+24*xi*(4.*xi-5.)+24.*eta*(4.*eta+8.*xi-5.)) ;

            ad2NdXi2( 0,  3 ) = c*eta*(48.*xi-12.) ;
            ad2NdXi2( 1,  3 ) = 0.0 ;
            ad2NdXi2( 2,  3 ) = c*(12.*xi*(2.*xi-1)+1.) ;

            ad2NdXi2( 0,  4 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 1,  4 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 2,  4 ) = b*(8.*eta-1.)*(8.*xi-1.) ;

            ad2NdXi2( 0,  5 ) = 0.0;
            ad2NdXi2( 1,  5 ) = c*xi*(48.*eta-12.) ;
            ad2NdXi2( 2,  5 ) = c+12.*c*eta*(2.*eta-1.) ;

            ad2NdXi2( 0,  6 ) = 0.0;
            ad2NdXi2( 1,  6 ) = 2.*c*(6*xi-6.*eta*(8.*eta+4.*xi-7.)-7) ;
            ad2NdXi2( 2,  6 ) = c*(12.*(1.-2.*eta)*eta-1.) ;

            ad2NdXi2( 0,  7 ) = 8.*b*eta*(4.*eta-1.) ;
            ad2NdXi2( 1,  7 ) = 2.*b*(19.+96.*eta*(xi+eta-1.)+4.*xi*(4.*xi-9.)) ;
            ad2NdXi2( 2,  7 ) = b*(7.-8.*xi+8.*eta*(12.*eta+8.*xi-9.)) ;

            ad2NdXi2( 0,  8 ) = c*eta*(36.-48.*(xi+eta)) ;
            ad2NdXi2( 1,  8 ) = 2.*c*(36.*xi-13.-6.*eta*(8.*eta-9.)-24.*xi*(3.*eta+xi)) ;
            ad2NdXi2( 2,  8 ) = c*(36.*xi-24.*(3.*eta*(eta-1.)+xi*(4.*eta+xi))-13.) ;

            ad2NdXi2( 0,  9 ) = 2.*c*(6.*xi*(9.-8.*xi)-eta*(24.*eta+36.*(2.*xi-1.))-13.) ;
            ad2NdXi2( 1,  9 ) = c*xi*(36-48*(xi+eta)) ;
            ad2NdXi2( 2,  9 ) = c*(72.*(1.-xi)*xi-12.*eta*(2.*eta+8.*xi-3.)-13.) ;

            ad2NdXi2( 0, 10 ) = 2.*b*(19.+96.*(xi-1.)*xi+eta*(16.*eta+12.*(8.*xi-3.))) ;
            ad2NdXi2( 1, 10 ) = 8.*b*xi*(4.*xi-1.) ;
            ad2NdXi2( 2, 10 ) = b*(7.+24.*xi*(4.*xi-3.)+8.*eta*(8.*xi-1.)) ;

            ad2NdXi2( 0, 11 ) = 2.*c*(6.*eta*(1.0-4.*xi)-6.*xi*(8.*xi-7.)-7.) ;
            ad2NdXi2( 1, 11 ) = 0.0 ;
            ad2NdXi2( 2, 11 ) = c*(12.*xi*(1.-2.*xi)-1.) ;

            ad2NdXi2( 0, 12 ) = d*eta*(10-24.*xi-8.*eta) ;
            ad2NdXi2( 1, 12 ) = 2.*d*xi*(1.0-4.*xi) ;
            ad2NdXi2( 2, 12 ) = d*(eta*(2.-16.*xi)+2.*(5. - 6.*xi)*xi -1.) ;

            ad2NdXi2( 0, 13 ) = 2.*d*eta*(1.-4.*eta) ;
            ad2NdXi2( 1, 13 ) =  d*xi*(10.-24.*eta-8.*xi) ;
            ad2NdXi2( 2, 13 ) = d*(2.*xi - 2.*eta*(6.*eta+8.*xi-5.)-1.) ;

            ad2NdXi2( 0, 14 ) = d*eta*(16.*eta+24.*xi-14.) ;
            ad2NdXi2( 1, 14 ) = d*xi*(24.*eta+16.*xi-14.) ;
            ad2NdXi2( 2, 14 ) = d*(3.+2*xi*(6.*xi-7.) + 2*eta*(6.*eta+16*xi-7.)) ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TRI15_HPP
