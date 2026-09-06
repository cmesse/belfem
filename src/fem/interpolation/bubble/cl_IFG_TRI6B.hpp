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

#ifndef BELFEM_CL_IFG_TRI6B_HPP
#define BELFEM_CL_IFG_TRI6B_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

    	 
    	 /**
    	 *  Bubble function for TRI6 element, side 1,
    	 *  following the EXODUSII numbering scheme
    	 * \f{align}{ N_1 & =\frac{32}{3} \, \eta \, \zeta \, \left(\eta^2-\zeta^2\right)  \\[2ex] 
    	 *            N_2 & =\frac{32}{3} \, \eta \, \zeta \, \left(\zeta^2-\eta^2\right) \f}.
    	 */
    	 
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                 GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >
        ::element_type() const
        {
            return ElementType::TRI6;
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 2, 2 );

			aXiHat( 0, 0 ) =  0.00 ;
			aXiHat( 0, 1 ) =  0.75 ;
			
			aXiHat( 1, 0 ) =  0.00 ;
			aXiHat( 1, 1 ) =  0.25 ;
        }
        
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {
            const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = 1. - xi - eta ;
			
            aN.set_size( 1, 2 );
            aN( 0, 0 ) = 32./3 * eta * zeta * ( eta*eta - zeta*zeta );
            aN( 0, 1 ) = -aN( 0, 0 );
        }
        

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
            const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = 1. - xi - eta ;
			
			real K    = 32./3. ;
	
			
			adNdXi.set_size( 2, 2 );
			
			//dN1/dxi
            adNdXi( 0, 0 ) = K * eta * ( 3. * zeta * zeta - eta * eta ) ;
            
             //dN1/deta
            adNdXi( 0, 1 ) = K * ((eta+zeta)*(2.*eta*zeta-(eta-zeta)*(eta-zeta))) ;
            
            adNdXi( 1, 0 ) = -adNdXi( 0, 0 );
            adNdXi( 1, 1 ) = -adNdXi( 0, 1 );
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
        	const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = 1. - xi - eta ;
			
            ad2NdXi2.set_size( 3, 2 );
            
            // d²N1/dxi²
 			ad2NdXi2( 0, 0 ) =  -64. * eta * zeta ;
 			
 			// d²N1/deta²
            ad2NdXi2( 1, 0 ) =   64. *(xi - 1.)*( eta - zeta );
            
            // d²N1/(dxi*deta)
            ad2NdXi2( 2, 0 ) =   32. * ( (zeta-eta)*(zeta-eta) - 2. * eta*eta);
            
            ad2NdXi2( 0, 1 ) = -ad2NdXi2( 0, 0 ) ;
            ad2NdXi2( 1, 1 ) = -ad2NdXi2( 1, 0 ) ;
            ad2NdXi2( 2, 1 ) = -ad2NdXi2( 2, 0 ) ;
            
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IFG_TRI6A_HPP
