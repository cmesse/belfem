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

#ifndef BELFEM_CL_IF_TRI6_HPP
#define BELFEM_CL_IF_TRI6_HPP


#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

    	/**
    	 * Lagrange function for quadratic Triangle element,
    	 * following the EXODUSII node numbering
    	 * \f{align}{ N_1 & = \xi \,\left(2\,\xi -1\right)  \\ 
    	 *            N_2 & = \eta \,\left(2\,\eta -1\right) \\
    	 *            N_3 & = \zeta \,\left(2\,\zeta -1\right) \\
    	 *            N_4 & = 4 \, \xi \, \eta \\
    	 *            N_5 & = 4 \, \eta \, \zeta \\
    	 *            N_6 & = 4 \, \xi \, \zeta \f}.
    	 */
    	 
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >
        ::element_type() const
        {
            return ElementType::TRI6;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 6 );

            aXiHat( 0, 0 ) = 1.0;
            aXiHat( 1, 0 ) = 0.0;

            aXiHat( 0, 1 ) = 0.0;
            aXiHat( 1, 1 ) = 1.0;

            aXiHat( 0, 2 ) = 0.0;
            aXiHat( 1, 2 ) = 0.0;

            aXiHat( 0, 3 ) = 0.5;
            aXiHat( 1, 3 ) = 0.5;

            aXiHat( 0, 4 ) = 0.0;
            aXiHat( 1, 4 ) = 0.5;

            aXiHat( 0, 5 ) = 0.5;
            aXiHat( 1, 5 ) = 0.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real   xi = aXi( 0 );
            const real  eta = aXi( 1 );
			const real zeta = 1.0 - xi - eta ;
			
            aN.set_size( 1, 6 );

            aN( 0, 0 ) =   xi * ( 2.0 *   xi - 1.0 );
            aN( 0, 1 ) =  eta * ( 2.0 *  eta - 1.0 );
            aN( 0, 2 ) = zeta * ( 2.0 * zeta - 1.0 );
            aN( 0, 3 ) = 4.0 * xi * eta;
            aN( 0, 4 ) = 4.0 * eta * zeta ;
            aN( 0, 5 ) = 4.0 * xi * zeta;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );
			const real zeta = 1.0 - xi - eta ;

            adNdXi.set_size( 2, 6 );

			
            adNdXi( 0, 0 ) = xi * 4.0 - 1.0;  // dN1/dxi
            adNdXi( 1, 0 ) = 0.0;             // dN1/deta


            adNdXi( 0, 1 ) = 0.0;             // dN2/dxi
            adNdXi( 1, 1 ) = eta * 4.0 - 1.0; // dN2/deta

            adNdXi( 0, 2 ) = 4.0 * ( xi + eta ) - 3.0; // dN3/dxi
            adNdXi( 1, 2 ) = 4.0 * ( xi + eta ) - 3.0; // dN3/deta

            adNdXi( 0, 3 ) = 4.0 * eta; // dN4/dxi
            adNdXi( 1, 3 ) = 4.0 * xi;  // dN4/deta

            adNdXi( 0, 4 ) =  - 4.0 * eta;            // dN5/dxi
            adNdXi( 1, 4 ) =    4.0 * ( zeta - eta ); // dN5/deta

            adNdXi( 0, 5 ) =  4.0 * ( zeta - xi ); // dN6/dxi
            adNdXi( 1, 5 ) = - 4.0 * xi;           // dN6/deta

        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            ad2NdXi2.set_size( 3, 6 );

            ad2NdXi2( 0, 0 ) =  4.0; // d²N1/dxi²
            ad2NdXi2( 1, 0 ) =  0.0; // d²N1/deta²
            ad2NdXi2( 2, 0 ) =  0.0; // d²N1/(dxi*deta)

            ad2NdXi2( 0, 1 ) =  0.0; // d²N2/dxi²
            ad2NdXi2( 1, 1 ) =  4.0; // d²N2/deta²
            ad2NdXi2( 2, 1 ) =  0.0; // d²N2/(dxi*deta)

            ad2NdXi2( 0, 2 ) =  4.0; // d²N3/dxi²
            ad2NdXi2( 1, 2 ) =  4.0; // d²N3/deta²
            ad2NdXi2( 2, 2 ) =  4.0; // d²N3/(dxi*deta)

            ad2NdXi2( 0, 3 ) =  0.0; // d²N4/dxi²
            ad2NdXi2( 1, 3 ) =  0.0; // d²N4/deta²
            ad2NdXi2( 2, 3 ) =  4.0; // d²N4/(dxi*deta)

            ad2NdXi2( 0, 4 ) =  0.0; // d²N5/dxi²
            ad2NdXi2( 1, 4 ) = -8.0; // d²N5/deta²
            ad2NdXi2( 2, 4 ) = -4.0; // d²N5/(dxi*deta)

            ad2NdXi2( 0, 5 ) = -8.0; // d²N6/dxi²
            ad2NdXi2( 1, 5 ) =  0.0; // d²N6/deta²
            ad2NdXi2( 2, 5 ) = -4.0; // d²N6/(dxi*deta)
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_TRI6_HPP
