/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_HESSIAN_HPP
#define BELFEM_FN_HESSIAN_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * Assemble the 9x9 chain-rule matrix H that maps physical first and second
     * shape-function derivatives onto reference-space first and second
     * derivatives for a 3D isoparametric element:
     *
     *     | dN/dxi   |       | dN/dx   |
     *     | d2N/dxi2 |  =  H | d2N/dx2 |
     *
     * Row ordering of H (reference space):
     *   - rows 0..2 : dN/dxi,   dN/deta,   dN/dzeta
     *   - rows 3..5 : d2N/dxi2, d2N/deta2, d2N/dzeta2
     *   - rows 6..8 : d2N/(deta dzeta), d2N/(dxi dzeta), d2N/(dxi deta)
     *
     * Column ordering of H (physical space):
     *   - cols 0..2 : dN/dx,   dN/dy,   dN/dz
     *   - cols 3..5 : d2N/dx2, d2N/dy2, d2N/dz2
     *   - cols 6..8 : d2N/(dy dz), d2N/(dx dz), d2N/(dx dy)
     *
     * Block structure:
     *   - upper-left  3x3 block : J
     *   - lower-left  6x3 block : K
     *   - upper-right 3x6 block : 0
     *   - lower-right 6x6 block : quadratic-in-J operator built from products
     *                             of Jacobian entries (chain rule)
     *
     * Physical derivatives are obtained from reference derivatives by solving
     * the corresponding 9x9 linear system with H (e.g. via H \ rhs).
     *
     * @param[in]  J   3x3 Jacobian with layout J(i,j) = dx_j / dxi_i,
     *                 where xi_0, xi_1, xi_2 = (xi, eta, zeta) and
     *                 x_0, x_1, x_2 = (x, y, z).
     * @param[in]  K   6x3 matrix of second derivatives of the geometry
     *                 mapping, K(i,j) = d2 x_j / d(xi-pair)_i, with row
     *                 pairing (xi2, eta2, zeta2, eta*zeta, xi*zeta, xi*eta).
     * @param[out] H   9x9 assembled transformation matrix (resized internally).
     */
    template < typename T >
    void
    hessian( const Matrix< T > & J, const Matrix< T > & K, Matrix< T > & H )
    {
        H.set_size( 9,9 );

        H( 0, 0 ) = J( 0, 0 );
        H( 1, 0 ) = J( 1, 0 );
        H( 2, 0 ) = J( 2, 0 );

        H( 3, 0 ) = K( 0, 0 );
        H( 4, 0 ) = K( 1, 0 );
        H( 5, 0 ) = K( 2, 0 );

        H( 6, 0 ) = K( 3, 0 );
        H( 7, 0 ) = K( 4, 0 );
        H( 8, 0 ) = K( 5, 0 );

        H( 0, 1 ) = J( 0, 1 );
        H( 1, 1 ) = J( 1, 1 );
        H( 2, 1 ) = J( 2, 1 );

        H( 3, 1 ) = K( 0, 1 );
        H( 4, 1 ) = K( 1, 1 );
        H( 5, 1 ) = K( 2, 1 );

        H( 6, 1 ) = K( 3, 1 );
        H( 7, 1 ) = K( 4, 1 );
        H( 8, 1 ) = K( 5, 1 );

        H( 0, 2 ) = J( 0, 2 );
        H( 1, 2 ) = J( 1, 2 );
        H( 2, 2 ) = J( 2, 2 );

        H( 3, 2 ) = K( 0, 2 );
        H( 4, 2 ) = K( 1, 2 );
        H( 5, 2 ) = K( 2, 2 );

        H( 6, 2 ) = K( 3, 2 );
        H( 7, 2 ) = K( 4, 2 );
        H( 8, 2 ) = K( 5, 2 );

        H( 0, 3 ) = 0. ;
        H( 1, 3 ) = 0. ;
        H( 2, 3 ) = 0. ;

        H( 3, 3 ) = J( 0, 0 ) * J( 0, 0 ) ;
        H( 4, 3 ) = J( 1, 0 ) * J( 1, 0 ) ;
        H( 5, 3 ) = J( 2, 0 ) * J( 2, 0 ) ;

        H( 6, 3 ) = J( 1, 0 ) * J( 2, 0 ) ;
        H( 7, 3 ) = J( 0, 0 ) * J( 2, 0 ) ;
        H( 8, 3 ) = J( 0, 0 ) * J( 1, 0 ) ;

        H( 0, 4 ) = 0. ;
        H( 1, 4 ) = 0. ;
        H( 2, 4 ) = 0. ;

        H( 3, 4 ) = J( 0, 1 ) * J( 0, 1 ) ;
        H( 4, 4 ) = J( 1, 1 ) * J( 1, 1 ) ;
        H( 5, 4 ) = J( 2, 1 ) * J( 2, 1 ) ;

        H( 6, 4 ) = J( 1, 1 ) * J( 2, 1 ) ;
        H( 7, 4 ) = J( 2, 1 ) * J( 0, 1 ) ;
        H( 8, 4 ) = J( 0, 1 ) * J( 1, 1 ) ;

        H( 0, 5 ) = 0. ;
        H( 1, 5 ) = 0. ;
        H( 2, 5 ) = 0. ;
        
        H( 3, 5 ) = J( 0, 2 ) * J( 0, 2 ) ;
        H( 4, 5 ) = J( 1, 2 ) * J( 1, 2 ) ;
        H( 5, 5 ) = J( 2, 2 ) * J( 2, 2 ) ;
        
        H( 6, 5 ) = J( 1, 2 ) * J( 2, 2 ) ;
        H( 7, 5 ) = J( 2, 2 ) * J( 0, 2 ) ;
        H( 8, 5 ) = J( 0, 2 ) * J( 1, 2 ) ;
        
        H( 0, 6 ) = 0. ;
        H( 1, 6 ) = 0. ;
        H( 2, 6 ) = 0. ;
        
        H( 3, 6 ) = 2. * J(0,1 ) * J ( 0, 2 ) ;
        H( 4, 6 ) = 2. * J(1,1 ) * J ( 1, 2 ) ;
        H( 5, 6 ) = 2. * J(2,1 ) * J ( 2, 2 ) ;

        H( 6, 6 ) = J(1,1 ) * J ( 2, 2 ) + J(2,1 ) * J ( 1, 2 );
        H( 7, 6 ) = J(0,1 ) * J ( 2, 2 ) + J(2,1 ) * J ( 0, 2 );
        H( 8, 6 ) = J(0,1 ) * J ( 1, 2 ) + J(1,1 ) * J ( 0, 2 );

        H( 0, 7 ) = 0. ;
        H( 1, 7 ) = 0. ;
        H( 2, 7 ) = 0. ;

        H( 3, 7 ) = 2. * J(0,0 ) * J ( 0, 2 ) ;
        H( 4, 7 ) = 2. * J(1,0 ) * J ( 1, 2 ) ;
        H( 5, 7 ) = 2. * J(2,0 ) * J ( 2, 2 ) ;

        H( 6, 7 ) = J(1,0 ) * J ( 2, 2 ) + J(2,0 ) * J ( 1, 2 );
        H( 7, 7 ) = J(0,0 ) * J ( 2, 2 ) + J(2,0 ) * J ( 0, 2 );
        H( 8, 7 ) = J(0,0 ) * J ( 1, 2 ) + J(1,0 ) * J ( 0, 2 );

        H( 0, 8 ) = 0. ;
        H( 1, 8 ) = 0. ;
        H( 2, 8 ) = 0. ;

        H( 3, 8 ) = 2. * J(0,0 ) * J ( 0, 1 ) ;
        H( 4, 8 ) = 2. * J(1,0 ) * J ( 1, 1 ) ;
        H( 5, 8 ) = 2. * J(2,0 ) * J ( 2, 1 ) ;

        H( 6, 8 ) = J(1,0 ) * J ( 2, 1 ) + J(2,0 ) * J ( 1, 1 );
        H( 7, 8 ) = J(0,0 ) * J ( 2, 1 ) + J(2,0 ) * J ( 0, 1 );
        H( 8, 8 ) = J(0,0 ) * J ( 1, 1 ) + J(1,0 ) * J ( 0, 1 );

    }
}
#endif //BELFEM_FN_HESSIAN_HPP
