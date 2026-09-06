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

#ifndef BELFEM_CL_IFG_TET4D_HPP
#define BELFEM_CL_IFG_TET4D_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

    	 /**
    	 *  Bubble function for TET4 element, side 3,
    	 *  following the EXODUSII numbering scheme
    	 * \f{align}{ N & = 27 \, \xi \, \eta \, \zeta \f}.
    	 */
    	 
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >
        ::interpolation_order() const
        {
            return InterpolationOrder::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                 GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >
        ::element_type() const
        {
            return ElementType::TET4;
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 1, 3 );

			aXiHat( 0, 0 ) =  1./3. ;
			aXiHat( 0, 1 ) =  1./3. ;
			aXiHat( 0, 2 ) =  1./3. ;

        }
        
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {
            const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = aXi( 2 );
			
            aN.set_size( 1, 1 );
            aN( 0, 0 ) = 27. * xi * eta *zeta ;
        }
        

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
            const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = aXi( 2 );

            adNdXi.set_size( 3, 1 );
            adNdXi( 0, 0 ) = 27. * eta * zeta ;  // dN/dxi
            adNdXi( 1, 0 ) = 27. * xi * zeta ;  // dN/deta
            adNdXi( 2, 0 ) = 27. * xi * eta ;  // dN/dzeta
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
        	const real   xi = aXi( 0 );
			const real  eta = aXi( 1 );
			const real zeta = aXi( 2 );
			
            ad2NdXi2.set_size( 6,1 );
            
            // d²N/dxi²
 			ad2NdXi2( 0, 0 ) =  0.0 ;
 			
 			// d²N/deta²
            ad2NdXi2( 1, 0 ) =  0.0 ;
            
            // d²N/dzeta²
            ad2NdXi2( 2, 0 ) =  0.0 ;
            
            // d²N/(deta*dzeta)
            ad2NdXi2( 3, 0 )  = 27. * xi ;
            
  			// d²N/(dxi*dzeta)
            ad2NdXi2( 4, 0 )  = 27. * eta ;
            
        	// d²N/(dxi*deta)
            ad2NdXi2( 5, 0 )  = 27. * zeta ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IFG_TET4D_HPP
