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

#ifndef BELFEM_CL_IFG_TET10D_HPP
#define BELFEM_CL_IFG_TET10D_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

    	 /**
    	 *  Bubble function for TET10 element, side 3,
    	 *  following the EXODUSII numbering scheme
    	 * \f{align}{ N_1 & = 64.8 \, \xi \, \eta \, \zeta \, \left( 2 \, \zeta^2 - \xi^2 - \eta^2 \right)\\[2ex]
    	 * N_2 & = 64.8 \, \xi \, \eta \, \zeta  \, \left( 2 \, \xi^2 - \zeta^2 - \eta^2 \right) \\[2ex]
    	 * N_3 & = 64.8 \, \xi \, \eta \, \zeta  \,  \left( 2 \, \eta^2 - \xi^2 - \zeta^2 \right) \f}.
    	 */
    	 
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                 GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >
        ::element_type() const
        {
            return ElementType::TET10;
        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >
        ::param_coords( Matrix <real> & aXiHat ) const
        {
            aXiHat.set_size( 3, 3 );

			// point 0
			aXiHat( 0, 0 ) =  2./3. ;
			aXiHat( 0, 1 ) =  1./6. ;
			aXiHat( 0, 2 ) =  1./6. ;

			// point 1
			aXiHat( 1, 0 ) =  1./6. ;
			aXiHat( 1, 1 ) =  1./6. ;
			aXiHat( 1, 2 ) =  2./3. ;
			
			// point 2
			aXiHat( 2, 0 ) =  1./6. ;
			aXiHat( 2, 1 ) =  2./3. ;
			aXiHat( 2, 2 ) =  1./6. ;
        }
        
//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >::N(
                const Vector <real> & aXi,
                Matrix <real> & aN ) const
        {
            const real    xi = aXi( 0 );
			const real   eta = aXi( 1 );
			const real  zeta = aXi( 2 );
			
			real K = 64.8 ;
			
			real  xi2  = xi*xi ;
			real  eta2 = eta*eta ;
			real zeta2 = zeta*zeta ;
			
			real a = 2*zeta2-xi2-eta2 ;
			real b = 2*xi2 - eta2-zeta2 ;
			real c = 2*eta2 - xi2-zeta2 ;
			real d = xi*eta*zeta ;
			
            aN.set_size( 1, 3 );
            aN( 0, 0 ) = K * d * a ;
            aN( 0, 1 ) = K * d * b ;
            aN( 0, 2 ) = K * d * c ;
        }
        

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >::dNdXi(
                const Vector <real> & aXi,
                Matrix <real> & adNdXi ) const
        {
			const real    xi = aXi( 0 );
			const real   eta = aXi( 1 );
			const real  zeta = aXi( 2 );
			
			real K = 64.8 ;
						
			real  xi2  = xi*xi ;
			real  eta2 = eta*eta ;
			real zeta2 = zeta*zeta ;
			
			real a = 2.*zeta2-xi2-eta2 ;
			real b = 2.*xi2 - eta2-zeta2 ;
			real c = 2.*eta2 - xi2-zeta2 ;
			real d = xi*eta*zeta ;
			
			real a_xi   = -2.*xi ;
			real a_eta  = -2.*eta ;
			real a_zeta =  4.*zeta ;
			
			real b_xi   = 4.* xi ;
			real b_eta  = -2.* eta ;
			real b_zeta = -2.*zeta ;
			
			real c_xi   = -2.*xi ;
			real c_eta  =  4.*eta ;
			real c_zeta =  -2.*zeta;
			
			real d_xi   = eta*zeta;
			real d_eta  = xi*zeta ;
			real d_zeta = xi*eta;
			
            adNdXi.set_size( 3, 3 );
            
            adNdXi( 0, 0 ) = K * ( d_xi   * a + d * a_xi   );  // dN1/dxi
            adNdXi( 0, 1 ) = K * ( d_eta  * a + d * a_eta  );  // dN1/deta
            adNdXi( 0, 2 ) = K * ( d_zeta * a + d * a_zeta );  // dN1/dzeta
            
            adNdXi( 1, 0 ) = K * ( d_xi   * b + d * b_xi   );  // dN2/dxi
            adNdXi( 1, 1 ) = K * ( d_eta  * b + d * b_eta  );  // dN2/deta
            adNdXi( 1, 2 ) = K * ( d_zeta * b + d * b_zeta );  // dN2/dzeta
          
	        adNdXi( 2, 0 ) = K * ( d_xi   * c + d * c_xi   );  // dN3/dxi
            adNdXi( 2, 1 ) = K * ( d_eta  * c + d * c_eta  );  // dN3/deta
            adNdXi( 2, 2 ) = K * ( d_zeta * c + d * c_zeta );  // dN3/dzeta

        }

//------------------------------------------------------------------------------


        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >::d2NdXi2(
                const Vector <real> & aXi,
                Matrix <real> & ad2NdXi2 ) const
        {
			const real    xi = aXi( 0 );
			const real   eta = aXi( 1 );
			const real  zeta = aXi( 2 );
			
			real K = 64.8 ;
			
			real  xi2  = xi*xi ;
			real  eta2 = eta*eta ;
			real zeta2 = zeta*zeta ;
			
			real a = 2.*zeta2-xi2-eta2 ;
			real b = 2.*xi2 - eta2-zeta2 ;
			real c = 2.*eta2 - xi2-zeta2 ;
			real d = xi*eta*zeta ;
			
			real a_xi   = -2.*xi ;
			real a_eta  = -2.*eta ;
			real a_zeta =  4.*zeta ;
			
			real b_xi   = 4.* xi ;
			real b_eta  = -2.* eta ;
			real b_zeta = -2.*zeta ;
			
			real c_xi   = -2.*xi ;
			real c_eta  =  4.*eta ;
			real c_zeta =  -2.*zeta;
			
			real d_xi   = eta*zeta;
			real d_eta  = xi*zeta ;
			real d_zeta = xi*eta;
			
            ad2NdXi2.set_size( 6,3 );
            
            // d²N1/dxi²
 			ad2NdXi2( 0, 0 ) =  2.*K*(a_xi*d_xi-d);
 			
 			// d²N1/deta²
            ad2NdXi2( 1, 0 ) =  2.*K*(a_eta*d_eta-d);
            
            // d²N1/dzeta²
            ad2NdXi2( 2, 0 ) =  2.*K*(a_zeta*d_zeta+2.*d);

            // d²N1/(deta*dzeta)
            ad2NdXi2( 3, 0 )  = K*(a_eta*d_zeta + a_zeta*d_eta + a*xi);
            
  			// d²N1/(dxi*dzeta)
            ad2NdXi2( 4, 0 )  = K*(a_xi*d_zeta + a_zeta*d_xi + a*eta);
            
        	// d²N1/(dxi*deta)
            ad2NdXi2( 5, 0 )  = K*(a_eta*d_xi + a_xi*d_eta + a*zeta);
            
            // d²N2/dxi²
 			ad2NdXi2( 0, 1 ) = 2.*K*(b_xi*d_xi+2.*d);
 			
 			// d²N2/deta²
            ad2NdXi2( 1, 1 ) = 2.*K*(b_eta*d_eta-d);
            
            // d²N2/dzeta²
            ad2NdXi2( 2, 1 ) = 2.*K*(b_zeta*d_zeta-d);

            // d²N2/(deta*dzeta)
            ad2NdXi2( 3, 1 )  = K*(b_eta*d_zeta + b_zeta*d_eta + b*xi);
            
  			// d²N2/(dxi*dzeta)
            ad2NdXi2( 4, 1 )  = K*(b_xi*d_zeta + b_zeta*d_xi + b*eta);
            
        	// d²N2/(dxi*deta)
            ad2NdXi2( 5, 1 )  = K*(b_eta*d_xi + b_xi*d_eta + b*zeta);

			// d²N3/dxi²
 			ad2NdXi2( 0, 2 ) =  2.*K*(c_xi*d_xi-d);
 			
 			// d²N3/deta²
            ad2NdXi2( 1, 2 ) =  2.*K*(c_eta*d_eta+2.*d);
            
            // d²N3/dzeta²
            ad2NdXi2( 2, 2 ) =  2.*K*(c_zeta*d_zeta-d);

            // d²N3/(deta*dzeta)
            ad2NdXi2( 3, 2 )  = K*(c_eta*d_zeta + c_zeta*d_eta + c*xi);
            
  			// d²N3/(dxi*dzeta)
            ad2NdXi2( 4, 2 )  = K*(c_xi*d_zeta + c_zeta*d_xi + c*eta);
            
        	// d²N3/(dxi*deta)
            ad2NdXi2( 5, 2 )  = K*(c_eta*d_xi + c_xi*d_eta + c*zeta);

        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IFG_TET4D_HPP
