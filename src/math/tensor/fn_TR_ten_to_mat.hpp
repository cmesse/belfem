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

#ifndef FN_BELFEM_TR_TO_MAT_HPP
#define FN_BELFEM_TR_TO_MAT_HPP

#include "cl_Matrix.hpp"

namespace belfem
{
    namespace tensor
    {
//----------------------------------------------------------------------------

        /**
         * elasticity matrix conversion
         * Elasticity tensor A => elasticity matrix C
         */
        template < typename T > void
        ten_to_mat( const T * A, T * C )
        {
            C[  0 ] = A[  0 ];
            C[  1 ] = A[  4 ];
            C[  2 ] = A[  8 ];
            C[  3 ] = A[  7 ];
            C[  4 ] = A[  2 ];
            C[  5 ] = A[  3 ];
            C[  6 ] = A[ 36 ];
            C[  7 ] = A[ 40 ];
            C[  8 ] = A[ 44 ];
            C[  9 ] = A[ 43 ];
            C[ 10 ] = A[ 38 ];
            C[ 11 ] = A[ 39 ];
            C[ 12 ] = A[ 72 ];
            C[ 13 ] = A[ 76 ];
            C[ 14 ] = A[ 80 ];
            C[ 15 ] = A[ 79 ];
            C[ 16 ] = A[ 74 ];
            C[ 17 ] = A[ 75 ];
            C[ 18 ] = A[ 63 ];
            C[ 19 ] = A[ 67 ];
            C[ 20 ] = A[ 71 ];
            C[ 21 ] = A[ 70 ];
            C[ 22 ] = A[ 65 ];
            C[ 23 ] = A[ 66 ];
            C[ 24 ] = A[ 18 ];
            C[ 25 ] = A[ 22 ];
            C[ 26 ] = A[ 26 ];
            C[ 27 ] = A[ 25 ];
            C[ 28 ] = A[ 20 ];
            C[ 29 ] = A[ 21 ];
            C[ 30 ] = A[ 27 ];
            C[ 31 ] = A[ 31 ];
            C[ 32 ] = A[ 35 ];
            C[ 33 ] = A[ 34 ];
            C[ 34 ] = A[ 29 ];
            C[ 35 ] = A[ 30 ];
        }

        template < typename T > void
        ten_to_mat( const T * A, Matrix< T > & C )
        {
#ifdef BELFEM_ARMADILLO
            ten_to_mat( A, C.data() );
#else
            C( 0, 0 ) = A[  0 ];
            C( 1, 0 ) = A[  4 ];
            C( 2, 0 ) = A[  8 ];
            C( 3, 0 ) = A[  7 ];
            C( 4, 0 ) = A[  2 ];
            C( 5, 0 ) = A[  3 ];
            C( 0, 1 ) = A[ 36 ];
            C( 1, 1 ) = A[ 40 ];
            C( 2, 1 ) = A[ 44 ];
            C( 3, 1 ) = A[ 43 ];
            C( 4, 1 ) = A[ 38 ];
            C( 5, 1 ) = A[ 39 ];
            C( 0, 2 ) = A[ 72 ];
            C( 1, 2 ) = A[ 76 ];
            C( 2, 2 ) = A[ 80 ];
            C( 3, 2 ) = A[ 79 ];
            C( 4, 2 ) = A[ 74 ];
            C( 5, 2 ) = A[ 75 ];
            C( 0, 3 ) = A[ 63 ];
            C( 1, 3 ) = A[ 67 ];
            C( 2, 3 ) = A[ 71 ];
            C( 3, 3 ) = A[ 70 ];
            C( 4, 3 ) = A[ 65 ];
            C( 5, 3 ) = A[ 66 ];
            C( 0, 4 ) = A[ 18 ];
            C( 1, 4 ) = A[ 22 ];
            C( 2, 4 ) = A[ 26 ];
            C( 3, 4 ) = A[ 25 ];
            C( 4, 4 ) = A[ 20 ];
            C( 5, 4 ) = A[ 21 ];
            C( 0, 5 ) = A[ 27 ];
            C( 1, 5 ) = A[ 31 ];
            C( 2, 5 ) = A[ 35 ];
            C( 3, 5 ) = A[ 34 ];
            C( 4, 5 ) = A[ 29 ];
            C( 5, 5 ) = A[ 30 ];
#endif
        }
//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */

#endif //FN_BELFEM_TR_TO_MAT_HPP
