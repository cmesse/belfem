//
// Created by Christian Messe on 26.07.20.
//

#ifndef BELFEM_CL_IF_TRUSS5_HPP
#define BELFEM_CL_IF_TRUSS5_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUARTIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >
        ::element_type() const
        {
            return ElementType::LINE5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 1, 5 );

            aXiHat( 0, 0 ) =  1.0 ;
            aXiHat( 0, 1 ) = -1.0 ;
            aXiHat( 0, 2 ) = -0.5 ;
            aXiHat( 0, 3 ) =  0.0 ;
            aXiHat( 0, 4 ) =  0.5 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real xi = aXi( 0 );

            aN.set_size( 1, 5 );

            aN( 0, 0 ) = xi*(1.-xi*(1.+xi*(4.-4.*xi)));
            aN( 0, 1 ) = xi*((4.*(xi+1.)*xi-1.)*xi-1.);
            aN( 0, 2 ) = xi*(((8.-16.*xi)*xi+16.)*xi-8.);
            aN( 0, 3 ) = (24.*xi*xi - 30.)*xi*xi + 6.;
            aN( 0, 4 ) = xi*(8.+xi*(16.-xi*(8.+16.*xi)));

            aN /= 6.0 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real xi = aXi( 0 );

            adNdXi.set_size( 1, 5 );

            adNdXi( 0, 0 ) = xi*((16.*xi-12.)*xi-2.)+1.;
            adNdXi( 0, 1 ) = xi*((16.*xi+12.)*xi-2.)-1.;
            adNdXi( 0, 2 ) = 8.*xi*(4.+(3.-8.*xi)*xi)-8.;
            adNdXi( 0, 3 ) = xi*(96.*xi*xi-60.);
            adNdXi( 0, 4 ) = 8.+ xi*(32.-xi*(64.*xi+24.));

            adNdXi /= 6.0 ;

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 1, 5 );

            const real xi = aXi( 0 );

            ad2NdXi2( 0, 0 ) = xi*(48.*xi-24.)-2.;
            ad2NdXi2( 0, 1 ) = xi*(48.*xi+24.)-2.;
            ad2NdXi2( 0, 2 ) = 32.+xi*(48.-192.*xi);
            ad2NdXi2( 0, 3 ) = 288.*xi*xi-60;
            ad2NdXi2( 0, 4 ) = 32.- xi*(192.*xi+48);

            ad2NdXi2 /= 6.0 ;

        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TRUSS5_HPP
