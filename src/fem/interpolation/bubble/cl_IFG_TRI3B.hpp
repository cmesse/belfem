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

#ifndef BELFEM_CL_IFG_TRI3B_HPP
#define BELFEM_CL_IFG_TRI3B_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------
    	/**
    	 * Bubble function for TRI3 element, side 1
    	 * \f$4\,\eta \,\zeta\f$
    	 */
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                 GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >
        ::element_type() const
        {
            return ElementType::TRI3;
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 2, 1 );

			aXiHat( 0, 0 ) =  0.0;
			aXiHat( 1, 0 ) =  0.5;
        }
        
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {
            const real  xi  = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = 1.0 - xi - eta ;
			
            aN.set_size( 1, 1 );
            aN( 0, 0 ) = 4.0 * eta * zeta ;
        }
        

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
            const real  xi  = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = 1.0 - xi - eta ;
            adNdXi.set_size( 2, 1 );
            adNdXi( 0, 0 ) = -4.0 * eta ;           // dN/dxi
            adNdXi( 1, 0 ) =  4.0 * (zeta-eta);     // dN/deta
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
            ad2NdXi2.set_size( 3, 1 );
 			ad2NdXi2( 0, 0 ) =  0.0 ; // d²N/dxi²
            ad2NdXi2( 1, 0 ) = -8.0 ; // d²N/deta²
            ad2NdXi2( 2, 0 ) = -4.0 ; // d²N/(dxi*deta)
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IFG_LINE2_HPP
