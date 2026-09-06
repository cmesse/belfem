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

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET10_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET10_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
// ----------------------------------------------------------------------------
namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        // 4th order in the Shunn and Ham numbering ( leading error term delta^4,
        // their Table 1 ): exact for polynomials up to degree 3. Weights are the
        // paper's ( Appendix F ) divided by 6, the reference tetrahedron volume.
        // UNUSED: fn_intpoints.cpp calls no tet10 table.
        // source Shunn and Ham, 10.1016/j.cam.2012.03.032
        inline void
        gauss_tet10(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {

            aPoints.set_size( 4, 10 );

            aPoints( 0 , 0 ) =  0.7784952948213300 ;
            aPoints( 1 , 0 ) =  0.0738349017262234 ;
            aPoints( 2 , 0 ) =  0.0738349017262234 ;

            aPoints( 0 , 1 ) =  0.0738349017262234 ;
            aPoints( 1 , 1 ) =  0.7784952948213300 ;
            aPoints( 2 , 1 ) =  0.0738349017262234 ;

            aPoints( 0 , 2 ) =  0.0738349017262234 ;
            aPoints( 1 , 2 ) =  0.0738349017262234 ;
            aPoints( 2 , 2 ) =  0.7784952948213300 ;

            aPoints( 0 , 3 ) =  0.0738349017262234 ;
            aPoints( 1 , 3 ) =  0.0738349017262234 ;
            aPoints( 2 , 3 ) =  0.0738349017262234 ;

            aPoints( 0 , 4 ) =  0.4062443438840510 ;
            aPoints( 1 , 4 ) =  0.4062443438840510 ;
            aPoints( 2 , 4 ) =  0.0937556561159491 ;

            aPoints( 0 , 5 ) =  0.4062443438840510 ;
            aPoints( 1 , 5 ) =  0.0937556561159491 ;
            aPoints( 2 , 5 ) =  0.4062443438840510 ;

            aPoints( 0 , 6 ) =  0.4062443438840510 ;
            aPoints( 1 , 6 ) =  0.0937556561159491 ;
            aPoints( 2 , 6 ) =  0.0937556561159491 ;

            aPoints( 0 , 7 ) =  0.0937556561159491 ;
            aPoints( 1 , 7 ) =  0.4062443438840510 ;
            aPoints( 2 , 7 ) =  0.4062443438840510 ;

            aPoints( 0 , 8 ) =  0.0937556561159491 ;
            aPoints( 1 , 8 ) =  0.4062443438840510 ;
            aPoints( 2 , 8 ) =  0.0937556561159491 ;

            aPoints( 0 , 9 ) =  0.0937556561159491 ;
            aPoints( 1 , 9 ) =  0.0937556561159491 ;
            aPoints( 2 , 9 ) =  0.4062443438840510 ;

            for( uint k=0; k<10; ++k )
            {
                aPoints( 3, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k )
                                   - aPoints( 2, k );
            }

            aWeights.set_size( 10 );

            aWeights( 0 ) = 0.0079388558072014826 ;
            aWeights( 1 ) = 0.0079388558072014826 ;
            aWeights( 2 ) = 0.0079388558072014826 ;
            aWeights( 3 ) = 0.0079388558072014826 ;
            aWeights( 4 ) = 0.022485207239643503 ;
            aWeights( 5 ) = 0.022485207239643503 ;
            aWeights( 6 ) = 0.022485207239643503 ;
            aWeights( 7 ) = 0.022485207239643503 ;
            aWeights( 8 ) = 0.022485207239643503 ;
            aWeights( 9 ) = 0.022485207239643503 ;
        }

// ----------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET10_HPP
