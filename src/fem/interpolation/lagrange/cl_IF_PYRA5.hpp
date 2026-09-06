/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_IF_PYRA5_HPP
#define BELFEM_CL_IF_PYRA5_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >
        ::interpolation_order() const
        {
            return InterpolationOrder::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >
        ::element_type() const
        {
            return ElementType::PYRA5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 3, 5 );

            aXiHat( 0, 0 ) = -1.0;
            aXiHat( 1, 0 ) = -1.0;
            aXiHat( 2, 0 ) =  0.0;

            aXiHat( 0, 1 ) =  1.0;
            aXiHat( 1, 1 ) = -1.0;
            aXiHat( 2, 1 ) =  0.0;

            aXiHat( 0, 2 ) =  1.0;
            aXiHat( 1, 2 ) =  1.0;
            aXiHat( 2, 2 ) =  0.0;

            aXiHat( 0, 3 ) = -1.0;
            aXiHat( 1, 3 ) =  1.0;
            aXiHat( 2, 3 ) =  0.0;

            aXiHat( 0, 4 ) =  0.0;
            aXiHat( 1, 4 ) =  0.0;
            aXiHat( 2, 4 ) =  1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >
                ::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
            const real zeta = aXi( 2 );
                  real xieta = xi*eta ;
                  real tau = 1.0 - zeta ;

            aN.set_size( 1, 5 );

            aN( 0, 0 ) = 0.25 * ( tau + xieta - xi - eta );
            aN( 0, 1 ) = 0.25 * ( tau - xieta + xi - eta );
            aN( 0, 2 ) = 0.25 * ( tau + xieta + xi + eta );
            aN( 0, 3 ) = 0.25 * ( tau - xieta - xi + eta );
            aN( 0, 4 ) = zeta ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >
                ::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );

            adNdXi.set_size( 3, 5 );

            adNdXi( 0, 0 ) = 0.25*eta - 0.25;
            adNdXi( 1, 0 ) = 0.25*xi - 0.25 ;
            adNdXi( 2, 0 ) = -0.25 ;

            adNdXi( 0, 1 ) =  0.25 - 0.25*eta ;
            adNdXi( 1, 1 ) = -0.25*xi - 0.25 ;
            adNdXi( 2, 1 ) = -0.25 ;

            adNdXi( 0, 2 ) = 0.25*eta + 0.25 ;
            adNdXi( 1, 2 ) = 0.25*xi + 0.250 ;
            adNdXi( 2, 2 ) = -0.25 ;

            adNdXi( 0, 3 ) = -0.25*eta - 0.25 ;
            adNdXi( 1, 3 ) = 0.25 - 0.25*xi ;
            adNdXi( 2, 3 ) = -0.25 ;

            adNdXi( 0, 4 ) = 0.0 ;
            adNdXi( 1, 4 ) = 0.0 ;
            adNdXi( 2, 4 ) = 1.0 ;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            // xi2 eta2 zeta2 eta*zeta xi*zeta xi*eta

            ad2NdXi2.set_size( 6, 5, 0.0 );

            ad2NdXi2( 5, 0 ) =  0.25 ;
            ad2NdXi2( 5, 1 ) = -0.25 ;
            ad2NdXi2( 5, 2 ) =  0.25 ;
            ad2NdXi2( 5, 3 ) = -0.25 ;
        }

//------------------------------------------------------------------------------
    }
}


#endif //BELFEM_CL_IF_PYRA5_HPP
