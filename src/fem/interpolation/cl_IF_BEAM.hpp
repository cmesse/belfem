//
// Created by Christian Messe on 26.10.19.
//

#ifndef BELFEM_CL_IF_BEAM_HPP
#define BELFEM_CL_IF_BEAM_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >
        ::element_type() const
        {
            return ElementType::LINE2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 1, 2 );

            aXiHat( 0, 0 ) = -1.0;
            aXiHat( 0, 1 ) =  1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real xi = aXi( 0 );

            const real a = (xi-1.0);
            const real b = (xi+1-0);

            aN.set_size( 1, 4 );

            aN( 0, 0 ) = a*a*(xi+2.0);
            aN( 0, 1 ) = a*a*b;
            aN( 0, 2 ) = b*b*(2.0-xi);
            aN( 0, 3 ) = b*b*a;

            aN *= 0.25;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real xi = aXi( 0 );
            const real xi2 = xi*xi;

            adNdXi.set_size( 1, 4 );
            adNdXi( 0, 0 ) =  3.0*(xi2-1.0);
            adNdXi( 0, 1 ) = (3.0*xi-2.0)*xi-1.0;
            adNdXi( 0, 2 ) = 3.0*(1.0-xi2);
            adNdXi( 0, 3 ) = (3.0*xi+2.0)*xi-1.0;
            adNdXi *= 0.25;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            const real xi = aXi( 0 );

            ad2NdXi2.set_size( 1, 4 );
            ad2NdXi2(0, 0 ) =  3.0*xi;
            ad2NdXi2(0, 1 ) =  3.0*xi-1.0;
            ad2NdXi2(0, 2 ) = -3.0*xi;
            ad2NdXi2(0, 3 ) =  3.0*xi+1.0;

            ad2NdXi2*=0.5;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_BEAM_HPP
