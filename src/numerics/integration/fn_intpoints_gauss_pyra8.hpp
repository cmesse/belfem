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

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PYRA8_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PYRA8_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"


namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * Conical-product Gauss rule on the reference pyramid ( base [-1,1]^2 at
         * zeta = 0, apex ( 0, 0, 1 ), as in cl_IF_PYRA5.hpp ): 2-point Gauss-Legendre
         * in xi and eta, scaled by ( 1 - zeta ), times 2-point Gauss-Jacobi(2,0) in
         * zeta, which absorbs the ( 1 - zeta )^2 volume factor. Exact for polynomials
         * up to degree 3; the weights sum to 4/3. Regenerated 2026-09-03: the
         * stored coordinates did not integrate degree 1 ( the weights already were
         * these ). Construction as in 10.1108/02644400410554362
         */
        inline void
        gauss_pyra8(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 3, 8 );
            aWeights.set_size( 8 );

            aPoints( 0,   0 ) = -0.50661630334978724 ;
            aPoints( 1,   0 ) = -0.50661630334978724 ;
            aPoints( 2,   0 ) =  0.1225148226554415 ;
            aWeights(   0 ) =  0.23254745125350801 ;

            aPoints( 0,   1 ) =  0.50661630334978724 ;
            aPoints( 1,   1 ) = -0.50661630334978724 ;
            aPoints( 2,   1 ) =  0.1225148226554415 ;
            aWeights(   1 ) =  0.23254745125350801 ;

            aPoints( 0,   2 ) = -0.50661630334978724 ;
            aPoints( 1,   2 ) =  0.50661630334978724 ;
            aPoints( 2,   2 ) =  0.1225148226554415 ;
            aWeights(   2 ) =  0.23254745125350801 ;

            aPoints( 0,   3 ) =  0.50661630334978724 ;
            aPoints( 1,   3 ) =  0.50661630334978724 ;
            aPoints( 2,   3 ) =  0.1225148226554415 ;
            aWeights(   3 ) =  0.23254745125350801 ;

            aPoints( 0,   4 ) = -0.26318405556971358 ;
            aPoints( 1,   4 ) = -0.26318405556971358 ;
            aPoints( 2,   4 ) =  0.54415184401122529 ;
            aWeights(   4 ) =  0.10078588207982532 ;

            aPoints( 0,   5 ) =  0.26318405556971358 ;
            aPoints( 1,   5 ) = -0.26318405556971358 ;
            aPoints( 2,   5 ) =  0.54415184401122529 ;
            aWeights(   5 ) =  0.10078588207982532 ;

            aPoints( 0,   6 ) = -0.26318405556971358 ;
            aPoints( 1,   6 ) =  0.26318405556971358 ;
            aPoints( 2,   6 ) =  0.54415184401122529 ;
            aWeights(   6 ) =  0.10078588207982532 ;

            aPoints( 0,   7 ) =  0.26318405556971358 ;
            aPoints( 1,   7 ) =  0.26318405556971358 ;
            aPoints( 2,   7 ) =  0.54415184401122529 ;
            aWeights(   7 ) =  0.10078588207982532 ;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PYRA8_HPP
