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

#ifndef BELFEM_CL_IFG_TET4A_HPP
#define BELFEM_CL_IFG_TET4A_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

    	 /**
    	 *  Bubble function for TET4 element, side 0,
    	 *  following the EXODUSII numbering scheme
    	 * \f{align}{ N & = 27 \, \xi \, \eta \, ( 1 - \xi - \eta ) \f}.
    	 */
    	 
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                 GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >
        ::element_type() const
        {
            return ElementType::TET4;
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 1, 3 );

			aXiHat( 0, 0 ) =  1./3. ;
			aXiHat( 0, 1 ) =  1./3. ;
			aXiHat( 0, 2 ) =  0.0 ;

        }
        
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {
            const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );

			
            aN.set_size( 1, 1 );
            aN( 0, 0 ) = 27. * xi * eta * ( 1. - xi - eta );
        }
        

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
            const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );

            adNdXi.set_size( 3, 1 );
            adNdXi( 0, 0 ) = 27. * eta * (1 - 2. * xi - eta );  // dN/dxi
            adNdXi( 1, 0 ) = 27. * xi * (1 - xi - 2*eta );  // dN/deta
            adNdXi( 2, 0 ) = 0.0 ;  // dN/dzeta
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
        	const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			
            ad2NdXi2.set_size( 6,1 );
            
            // d²N/dxi²
 			ad2NdXi2( 0, 0 ) =  -54. * xi ;
 			
 			// d²N/deta²
            ad2NdXi2( 1, 0 ) =  -54. * eta ;
            
            // d²N/dzeta²
            ad2NdXi2( 2, 0 ) =  0.0 ;
            
            // d²N/(deta*dzeta)
            ad2NdXi2( 3, 0 )  = 0.0 ;
            
  			// d²N/(dxi*dzeta)
            ad2NdXi2( 4, 0 )  = 0.0 ;
            
        	// d²N/(dxi*deta)
            ad2NdXi2( 5, 0 )  = 27. - 54. * ( xi + eta );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IFG_TET4A_HPP
