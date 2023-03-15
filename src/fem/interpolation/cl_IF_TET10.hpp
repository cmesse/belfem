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

            const real   xi2 = xi + xi ;
            const real  eta2 = eta + eta ;
            const real zeta2 = zeta + zeta ;
            const real  tau2 = tau + tau ;

            aN.set_size( 1, 10 );

            aN( 0,  0 ) = xi*(xi2-1.0);
            aN( 0,  1 ) = zeta*(zeta2-1.0);
            aN( 0,  2 ) = eta*(eta2-1.0);
            aN( 0,  3 ) = tau*(tau2-1.0);
            aN( 0,  4 ) = xi2*zeta2;
            aN( 0,  5 ) = eta2*zeta2;
            aN( 0,  6 ) = xi2*eta2;
            aN( 0,  7 ) = xi2*tau2;
            aN( 0,  8 ) = zeta2*tau2;
            aN( 0,  9 ) = eta2*tau2;
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

            const real xi4   = xi+xi+xi+xi ;
            const real eta4  = eta+eta+eta+eta ;
            const real zeta4 = zeta+zeta+zeta+zeta ;

            adNdXi.set_size( 3, 10 );

            adNdXi( 0, 0 ) = xi4-1.0;
            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 2, 0 ) = 0.0;

            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 1, 1 ) = 0.0 ;
            adNdXi( 2, 1 ) = zeta4-1.0;

            adNdXi( 0, 2 ) = 0.0;
            adNdXi( 1, 2 ) = eta4-1.0;
            adNdXi( 2, 2 ) = 0.0 ;

            adNdXi( 0, 3 ) = (xi4+eta4+zeta4)-3.0;
            adNdXi( 1, 3 ) = (xi4+eta4+zeta4)-3.0;
            adNdXi( 2, 3 ) = (xi4+eta4+zeta4)-3.0;

            adNdXi( 0, 4 ) = zeta4 ;
            adNdXi( 1, 4 ) = 0.0 ;
            adNdXi( 2, 4 ) = xi4 ;

            adNdXi( 0, 5 ) = 0.0;
            adNdXi( 1, 5 ) = zeta4 ;
            adNdXi( 2, 5 ) = eta4 ;

            adNdXi( 0, 6 ) = eta4 ;
            adNdXi( 1, 6 ) = xi4 ;
            adNdXi( 2, 6 ) = 0.0 ;

            adNdXi(0,  7 ) = 4. - xi4 - xi4 - eta4 - zeta4 ;
            adNdXi( 1, 7 ) = -xi4 ;
            adNdXi( 2, 7 ) = -xi4 ;

            adNdXi( 0, 8 ) = -zeta4 ;
            adNdXi( 1, 8 ) = -zeta4 ;
            adNdXi( 2, 8 ) = 4. - xi4 -  eta4 - zeta4- zeta4 ;

            adNdXi( 0, 9 ) = -eta4 ;
            adNdXi( 1, 9 ) = 4. - xi4 -  eta4 - eta4- zeta4 ;
            adNdXi( 2, 9 ) = -eta4 ;
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
            ad2NdXi2( 1, 1 ) =  0.0;
            ad2NdXi2( 2, 1 ) =  4.0;
            ad2NdXi2( 3, 1 ) =  0.0;
            ad2NdXi2( 4, 1 ) =  0.0;
            ad2NdXi2( 5, 1 ) =  0.0;

            ad2NdXi2( 0, 2 ) =  0.0;
            ad2NdXi2( 1, 2 ) =  4.0;
            ad2NdXi2( 2, 2 ) =  0.0;
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
            ad2NdXi2( 4, 4 ) =  4.0;
            ad2NdXi2( 5, 4 ) =  0.0;

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
            ad2NdXi2( 4, 6 ) =  0.0;
            ad2NdXi2( 5, 6 ) =  4.0;

            ad2NdXi2( 0, 7 ) = -8.0;
            ad2NdXi2( 1, 7 ) =  0.0;
            ad2NdXi2( 2, 7 ) =  0.0;
            ad2NdXi2( 3, 7 ) =  0.0;
            ad2NdXi2( 4, 7 ) = -4.0;
            ad2NdXi2( 5, 7 ) = -4.0;

            ad2NdXi2( 0, 8 ) =  0.0;
            ad2NdXi2( 1, 8 ) =  0.0;
            ad2NdXi2( 2, 8 ) = -8.0;
            ad2NdXi2( 3, 8 ) = -4.0;
            ad2NdXi2( 4, 8 ) = -4.0;
            ad2NdXi2( 5, 8 ) =  0.0;

            ad2NdXi2( 0, 9 ) =  0.0;
            ad2NdXi2( 1, 9 ) = -8.0;
            ad2NdXi2( 2, 9 ) =  0.0;
            ad2NdXi2( 3, 9 ) = -4.0;
            ad2NdXi2( 4, 9 ) =  0.0;
            ad2NdXi2( 5, 9 ) = -4.0;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TET10_HPP
