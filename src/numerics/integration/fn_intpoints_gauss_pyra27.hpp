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

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PYRA27_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PYRA27_HPP

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
         * zeta = 0, apex ( 0, 0, 1 ), as in cl_IF_PYRA5.hpp ): 3-point Gauss-Legendre
         * in xi and eta, scaled by ( 1 - zeta ), times 3-point Gauss-Jacobi(2,0) in
         * zeta, which absorbs the ( 1 - zeta )^2 volume factor. Exact for polynomials
         * up to degree 5; the weights sum to 4/3. Regenerated 2026-09-03: the
         * stored coordinates did not integrate degree 1 ( the weights already were
         * these ). Construction as in 10.1108/02644400410554362
         */
        inline void
        gauss_pyra27(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 3, 27 );
            aWeights.set_size( 27 );

            aPoints( 0,   0 ) = -0.71805574131988903 ;
            aPoints( 1,   0 ) = -0.71805574131988903 ;
            aPoints( 2,   0 ) =  0.072994024073149699 ;
            aWeights(   0 ) =  0.04849887687187858 ;

            aPoints( 0,   1 ) =  0 ;
            aPoints( 1,   1 ) = -0.71805574131988903 ;
            aPoints( 2,   1 ) =  0.072994024073149699 ;
            aWeights(   1 ) =  0.077598202995005638 ;

            aPoints( 0,   2 ) =  0.71805574131988903 ;
            aPoints( 1,   2 ) = -0.71805574131988903 ;
            aPoints( 2,   2 ) =  0.072994024073149699 ;
            aWeights(   2 ) =  0.04849887687187858 ;

            aPoints( 0,   3 ) = -0.71805574131988903 ;
            aPoints( 1,   3 ) =  0 ;
            aPoints( 2,   3 ) =  0.072994024073149699 ;
            aWeights(   3 ) =  0.077598202995005638 ;

            aPoints( 0,   4 ) =  0 ;
            aPoints( 1,   4 ) =  0 ;
            aPoints( 2,   4 ) =  0.072994024073149699 ;
            aWeights(   4 ) =  0.1241571247920089 ;

            aPoints( 0,   5 ) =  0.71805574131988903 ;
            aPoints( 1,   5 ) =  0 ;
            aPoints( 2,   5 ) =  0.072994024073149699 ;
            aWeights(   5 ) =  0.077598202995005638 ;

            aPoints( 0,   6 ) = -0.71805574131988903 ;
            aPoints( 1,   6 ) =  0.71805574131988903 ;
            aPoints( 2,   6 ) =  0.072994024073149699 ;
            aWeights(   6 ) =  0.04849887687187858 ;

            aPoints( 0,   7 ) =  0 ;
            aPoints( 1,   7 ) =  0.71805574131988903 ;
            aPoints( 2,   7 ) =  0.072994024073149699 ;
            aWeights(   7 ) =  0.077598202995005638 ;

            aPoints( 0,   8 ) =  0.71805574131988903 ;
            aPoints( 1,   8 ) =  0.71805574131988903 ;
            aPoints( 2,   8 ) =  0.072994024073149699 ;
            aWeights(   8 ) =  0.04849887687187858 ;

            aPoints( 0,   9 ) = -0.50580870785392507 ;
            aPoints( 1,   9 ) = -0.50580870785392507 ;
            aPoints( 2,   9 ) =  0.34700376603835181 ;
            aWeights(   9 ) =  0.045137737425884644 ;

            aPoints( 0,  10 ) =  0 ;
            aPoints( 1,  10 ) = -0.50580870785392507 ;
            aPoints( 2,  10 ) =  0.34700376603835181 ;
            aWeights(  10 ) =  0.072220379881415345 ;

            aPoints( 0,  11 ) =  0.50580870785392507 ;
            aPoints( 1,  11 ) = -0.50580870785392507 ;
            aPoints( 2,  11 ) =  0.34700376603835181 ;
            aWeights(  11 ) =  0.045137737425884644 ;

            aPoints( 0,  12 ) = -0.50580870785392507 ;
            aPoints( 1,  12 ) =  0 ;
            aPoints( 2,  12 ) =  0.34700376603835181 ;
            aWeights(  12 ) =  0.072220379881415345 ;

            aPoints( 0,  13 ) =  0 ;
            aPoints( 1,  13 ) =  0 ;
            aPoints( 2,  13 ) =  0.34700376603835181 ;
            aWeights(  13 ) =  0.11555260781026443 ;

            aPoints( 0,  14 ) =  0.50580870785392507 ;
            aPoints( 1,  14 ) =  0 ;
            aPoints( 2,  14 ) =  0.34700376603835181 ;
            aWeights(  14 ) =  0.072220379881415345 ;

            aPoints( 0,  15 ) = -0.50580870785392507 ;
            aPoints( 1,  15 ) =  0.50580870785392507 ;
            aPoints( 2,  15 ) =  0.34700376603835181 ;
            aWeights(  15 ) =  0.045137737425884644 ;

            aPoints( 0,  16 ) =  0 ;
            aPoints( 1,  16 ) =  0.50580870785392507 ;
            aPoints( 2,  16 ) =  0.34700376603835181 ;
            aWeights(  16 ) =  0.072220379881415345 ;

            aPoints( 0,  17 ) =  0.50580870785392507 ;
            aPoints( 1,  17 ) =  0.50580870785392507 ;
            aPoints( 2,  17 ) =  0.34700376603835181 ;
            aWeights(  17 ) =  0.045137737425884644 ;

            aPoints( 0,  18 ) = -0.22850430565396737 ;
            aPoints( 1,  18 ) = -0.22850430565396737 ;
            aPoints( 2,  18 ) =  0.70500220988849838 ;
            aWeights(  18 ) =  0.0092440441384508461 ;

            aPoints( 0,  19 ) =  0 ;
            aPoints( 1,  19 ) = -0.22850430565396737 ;
            aPoints( 2,  19 ) =  0.70500220988849838 ;
            aWeights(  19 ) =  0.014790470621521336 ;

            aPoints( 0,  20 ) =  0.22850430565396737 ;
            aPoints( 1,  20 ) = -0.22850430565396737 ;
            aPoints( 2,  20 ) =  0.70500220988849838 ;
            aWeights(  20 ) =  0.0092440441384508461 ;

            aPoints( 0,  21 ) = -0.22850430565396737 ;
            aPoints( 1,  21 ) =  0 ;
            aPoints( 2,  21 ) =  0.70500220988849838 ;
            aWeights(  21 ) =  0.014790470621521336 ;

            aPoints( 0,  22 ) =  0 ;
            aPoints( 1,  22 ) =  0 ;
            aPoints( 2,  22 ) =  0.70500220988849838 ;
            aWeights(  22 ) =  0.023664752994434112 ;

            aPoints( 0,  23 ) =  0.22850430565396737 ;
            aPoints( 1,  23 ) =  0 ;
            aPoints( 2,  23 ) =  0.70500220988849838 ;
            aWeights(  23 ) =  0.014790470621521336 ;

            aPoints( 0,  24 ) = -0.22850430565396737 ;
            aPoints( 1,  24 ) =  0.22850430565396737 ;
            aPoints( 2,  24 ) =  0.70500220988849838 ;
            aWeights(  24 ) =  0.0092440441384508461 ;

            aPoints( 0,  25 ) =  0 ;
            aPoints( 1,  25 ) =  0.22850430565396737 ;
            aPoints( 2,  25 ) =  0.70500220988849838 ;
            aWeights(  25 ) =  0.014790470621521336 ;

            aPoints( 0,  26 ) =  0.22850430565396737 ;
            aPoints( 1,  26 ) =  0.22850430565396737 ;
            aPoints( 2,  26 ) =  0.70500220988849838 ;
            aWeights(  26 ) =  0.0092440441384508461 ;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PYRA27_HPP
