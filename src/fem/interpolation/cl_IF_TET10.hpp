//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_TET10_HPP
#define BELFEM_CL_IF_TET10_HPP

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >
        ::element_type() const
        {
            return ElementType::TET10;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 10 );

            aXiHat(0, 0) = 1.0;
            aXiHat(1, 0) = 0.0;
            aXiHat(2, 0) = 0.0;

            aXiHat(0, 1) = 0.0;
            aXiHat(1, 1) = 1.0;
            aXiHat(2, 1) = 0.0;

            aXiHat(0, 2) = 0.0;
            aXiHat(1, 2) = 0.0;
            aXiHat(2, 2) = 1.0;

            aXiHat(0, 3) = 0.0;
            aXiHat(1, 3) = 0.0;
            aXiHat(2, 3) = 0.0;

            aXiHat(0, 4) = 0.5;
            aXiHat(1, 4) = 0.5;
            aXiHat(2, 0) = 0.0;

            aXiHat(0, 5) = 0.0;
            aXiHat(1, 5) = 0.5;
            aXiHat(2, 5) = 0.5;

            aXiHat(0, 6) = 0.5;
            aXiHat(1, 6) = 0.0;
            aXiHat(2, 6) = 0.5;

            aXiHat(0, 7) = 0.5;
            aXiHat(1, 7) = 0.0;
            aXiHat(2, 7) = 0.0;

            aXiHat(0, 8) = 0.0;
            aXiHat(1, 8) = 0.5;
            aXiHat(2, 8) = 0.0;

            aXiHat(0, 9) = 0.0;
            aXiHat(1, 9) = 0.0;
            aXiHat(2, 9) = 0.5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );
            const real  tau = 1.0 - xi - eta - zeta;

            aN.set_size( 1, 10 );

            aN( 0,  0 ) = xi*(2.0*xi-1.0);
            aN( 0,  1 ) = eta*(2.0*eta-1.0);
            aN( 0,  2 ) = zeta*(2.0*zeta-1.0);
            aN( 0,  3 ) = tau*(1.0-2.0*(xi+eta+zeta));
            aN( 0,  4 ) = 4.0*xi*eta;
            aN( 0,  5 ) = 4.0*eta*zeta;
            aN( 0,  6 ) = 4.0*xi*zeta;
            aN( 0,  7 ) = 4.0*xi*tau;
            aN( 0,  8 ) = 4.0*eta*tau;
            aN( 0,  9 ) = 4.0*zeta*tau;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );

            adNdXi.set_size( 3, 10 );

            adNdXi( 0, 0 ) = 4.0*xi-1.0;
            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 2, 0 ) = 0.0;

            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 1, 1 ) = 4.0*eta-1.0;
            adNdXi( 2, 1 ) = 0.0;

            adNdXi( 0, 2 ) = 0.0;
            adNdXi( 1, 2 ) = 0.0;
            adNdXi( 2, 2 ) = 4.0*zeta-1.0;

            adNdXi( 0, 3 ) = 4.0*(xi+eta+zeta)-3.0;
            adNdXi( 1, 3 ) = 4.0*(xi+eta+zeta)-3.0;
            adNdXi( 2, 3 ) = 4.0*(xi+eta+zeta)-3.0;

            adNdXi( 0, 4 ) = 4.0*eta;
            adNdXi( 1, 4 ) = 4.0*xi;
            adNdXi( 2, 4 ) = 0.0;

            adNdXi( 0, 5 ) = 0.0;
            adNdXi( 1, 5 ) = 4.0*zeta;
            adNdXi( 2, 5 ) = 4.0*eta;

            adNdXi( 0, 6 ) = 4.0*zeta;
            adNdXi( 1, 6 ) = 0.0;
            adNdXi( 2, 6 ) = 4.0*xi;

            adNdXi( 0, 7 ) = 4.0*(1.0-2.0*xi-eta-zeta);
            adNdXi( 1, 7 ) = -4.0*xi;
            adNdXi( 2, 7 ) = -4.0*xi;

            adNdXi( 0, 8 ) = -4.0*eta;
            adNdXi( 1, 8 ) = 4.0*(1.0-xi-2.0*eta-zeta);
            adNdXi( 2, 8 ) = -4.0*eta;

            adNdXi( 0, 9 ) = -4.0*zeta;
            adNdXi( 1, 9 ) = -4.0*zeta;
            adNdXi( 2, 9 ) = 4.0*(1.0-xi-eta-2.0*zeta);
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 6, 10 );

            ad2NdXi2( 0, 0 ) =  4.0;
            ad2NdXi2( 1, 0 ) =  0.0;
            ad2NdXi2( 2, 0 ) =  0.0;
            ad2NdXi2( 3, 0 ) =  0.0;
            ad2NdXi2( 4, 0 ) =  0.0;
            ad2NdXi2( 5, 0 ) =  0.0;

            ad2NdXi2( 0, 1 ) =  0.0;
            ad2NdXi2( 1, 1 ) =  4.0;
            ad2NdXi2( 2, 1 ) =  0.0;
            ad2NdXi2( 3, 1 ) =  0.0;
            ad2NdXi2( 4, 1 ) =  0.0;
            ad2NdXi2( 5, 1 ) =  0.0;

            ad2NdXi2( 0, 2 ) =  0.0;
            ad2NdXi2( 1, 2 ) =  0.0;
            ad2NdXi2( 2, 2 ) =  4.0;
            ad2NdXi2( 3, 2 ) =  0.0;
            ad2NdXi2( 4, 2 ) =  0.0;
            ad2NdXi2( 5, 2 ) =  0.0;

            ad2NdXi2( 0, 3 ) =  4.0;
            ad2NdXi2( 1, 3 ) =  4.0;
            ad2NdXi2( 2, 3 ) =  4.0;
            ad2NdXi2( 3, 3 ) =  4.0;
            ad2NdXi2( 4, 3 ) =  4.0;
            ad2NdXi2( 5, 3 ) =  4.0;

            ad2NdXi2( 0, 4 ) =  0.0;
            ad2NdXi2( 1, 4 ) =  0.0;
            ad2NdXi2( 2, 4 ) =  0.0;
            ad2NdXi2( 3, 4 ) =  0.0;
            ad2NdXi2( 4, 4 ) =  0.0;
            ad2NdXi2( 5, 4 ) =  4.0;

            ad2NdXi2( 0, 5 ) =  0.0;
            ad2NdXi2( 1, 5 ) =  0.0;
            ad2NdXi2( 2, 5 ) =  0.0;
            ad2NdXi2( 3, 5 ) =  4.0;
            ad2NdXi2( 4, 5 ) =  0.0;
            ad2NdXi2( 5, 5 ) =  0.0;

            ad2NdXi2( 0, 6 ) =  0.0;
            ad2NdXi2( 1, 6 ) =  0.0;
            ad2NdXi2( 2, 6 ) =  0.0;
            ad2NdXi2( 3, 6 ) =  0.0;
            ad2NdXi2( 4, 6 ) =  4.0;
            ad2NdXi2( 5, 6 ) =  0.0;

            ad2NdXi2( 0, 7 ) = -8.0;
            ad2NdXi2( 1, 7 ) =  0.0;
            ad2NdXi2( 2, 7 ) =  0.0;
            ad2NdXi2( 3, 7 ) =  0.0;
            ad2NdXi2( 4, 7 ) = -4.0;
            ad2NdXi2( 5, 7 ) = -4.0;

            ad2NdXi2( 0, 8 ) =  0.0;
            ad2NdXi2( 1, 8 ) = -8.0;
            ad2NdXi2( 2, 8 ) =  0.0;
            ad2NdXi2( 3, 8 ) = -4.0;
            ad2NdXi2( 4, 8 ) =  0.0;
            ad2NdXi2( 5, 8 ) = -4.0;

            ad2NdXi2( 0, 9 ) =  0.0;
            ad2NdXi2( 1, 9 ) =  0.0;
            ad2NdXi2( 2, 9 ) = -8.0;
            ad2NdXi2( 3, 9 ) = -4.0;
            ad2NdXi2( 4, 9 ) = -4.0;
            ad2NdXi2( 5, 9 ) =  0.0;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TET10_HPP
