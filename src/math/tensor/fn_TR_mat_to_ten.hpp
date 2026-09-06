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

#ifndef BELFEM_FN_TR_MAT_TO_TEN_HPP
#define BELFEM_FN_TR_MAT_TO_TEN_HPP

#include "cl_Matrix.hpp"
namespace belfem
{
    namespace tensor
    {
//----------------------------------------------------------------------------

        /**
         * elasticity matrix conversion
         * elasticity matrix C => tensor A
         */
        template < typename T > void
        mat_to_ten( const T * C, T * A )
        {
            A[  0 ] = C[  0 ];
            A[  1 ] = C[  5 ];
            A[  2 ] = C[  4 ];
            A[  3 ] = C[  5 ];
            A[  4 ] = C[  1 ];
            A[  5 ] = C[  3 ];
            A[  6 ] = C[  4 ];
            A[  7 ] = C[  3 ];
            A[  8 ] = C[  2 ];
            A[  9 ] = C[ 30 ];
            A[ 10 ] = C[ 35 ];
            A[ 11 ] = C[ 34 ];
            A[ 12 ] = C[ 35 ];
            A[ 13 ] = C[ 31 ];
            A[ 14 ] = C[ 33 ];
            A[ 15 ] = C[ 34 ];
            A[ 16 ] = C[ 33 ];
            A[ 17 ] = C[ 32 ];
            A[ 18 ] = C[ 24 ];
            A[ 19 ] = C[ 29 ];
            A[ 20 ] = C[ 28 ];
            A[ 21 ] = C[ 29 ];
            A[ 22 ] = C[ 25 ];
            A[ 23 ] = C[ 27 ];
            A[ 24 ] = C[ 28 ];
            A[ 25 ] = C[ 27 ];
            A[ 26 ] = C[ 26 ];
            A[ 27 ] = C[ 30 ];
            A[ 28 ] = C[ 35 ];
            A[ 29 ] = C[ 34 ];
            A[ 30 ] = C[ 35 ];
            A[ 31 ] = C[ 31 ];
            A[ 32 ] = C[ 33 ];
            A[ 33 ] = C[ 34 ];
            A[ 34 ] = C[ 33 ];
            A[ 35 ] = C[ 32 ];
            A[ 36 ] = C[  6 ];
            A[ 37 ] = C[ 11 ];
            A[ 38 ] = C[ 10 ];
            A[ 39 ] = C[ 11 ];
            A[ 40 ] = C[  7 ];
            A[ 41 ] = C[  9 ];
            A[ 42 ] = C[ 10 ];
            A[ 43 ] = C[  9 ];
            A[ 44 ] = C[  8 ];
            A[ 45 ] = C[ 18 ];
            A[ 46 ] = C[ 23 ];
            A[ 47 ] = C[ 22 ];
            A[ 48 ] = C[ 23 ];
            A[ 49 ] = C[ 19 ];
            A[ 50 ] = C[ 21 ];
            A[ 51 ] = C[ 22 ];
            A[ 52 ] = C[ 21 ];
            A[ 53 ] = C[ 20 ];
            A[ 54 ] = C[ 24 ];
            A[ 55 ] = C[ 29 ];
            A[ 56 ] = C[ 28 ];
            A[ 57 ] = C[ 29 ];
            A[ 58 ] = C[ 25 ];
            A[ 59 ] = C[ 27 ];
            A[ 60 ] = C[ 28 ];
            A[ 61 ] = C[ 27 ];
            A[ 62 ] = C[ 26 ];
            A[ 63 ] = C[ 18 ];
            A[ 64 ] = C[ 23 ];
            A[ 65 ] = C[ 22 ];
            A[ 66 ] = C[ 23 ];
            A[ 67 ] = C[ 19 ];
            A[ 68 ] = C[ 21 ];
            A[ 69 ] = C[ 22 ];
            A[ 70 ] = C[ 21 ];
            A[ 71 ] = C[ 20 ];
            A[ 72 ] = C[ 12 ];
            A[ 73 ] = C[ 17 ];
            A[ 74 ] = C[ 16 ];
            A[ 75 ] = C[ 17 ];
            A[ 76 ] = C[ 13 ];
            A[ 77 ] = C[ 15 ];
            A[ 78 ] = C[ 16 ];
            A[ 79 ] = C[ 15 ];
            A[ 80 ] = C[ 14 ];
        }

//----------------------------------------------------------------------------

        template < typename T > void
        mat_to_ten( const Matrix< T > & C , T * A )
        {
#ifdef BELFEM_ARMADILLO
            mat_to_ten( C.data(), A );
#else
            A[  0 ] = C( 0, 0 );
            A[  1 ] = C( 5, 0 );
            A[  2 ] = C( 4, 0 );
            A[  3 ] = C( 5, 0 );
            A[  4 ] = C( 1, 0 );
            A[  5 ] = C( 3, 0 );
            A[  6 ] = C( 4, 0 );
            A[  7 ] = C( 3, 0 );
            A[  8 ] = C( 2, 0 );
            A[  9 ] = C( 0, 5 );
            A[ 10 ] = C( 5, 5 );
            A[ 11 ] = C( 4, 5 );
            A[ 12 ] = C( 5, 5 );
            A[ 13 ] = C( 1, 5 );
            A[ 14 ] = C( 3, 5 );
            A[ 15 ] = C( 4, 5 );
            A[ 16 ] = C( 3, 5 );
            A[ 17 ] = C( 2, 5 );
            A[ 18 ] = C( 0, 4 );
            A[ 19 ] = C( 5, 4 );
            A[ 20 ] = C( 4, 4 );
            A[ 21 ] = C( 5, 4 );
            A[ 22 ] = C( 1, 4 );
            A[ 23 ] = C( 3, 4 );
            A[ 24 ] = C( 4, 4 );
            A[ 25 ] = C( 3, 4 );
            A[ 26 ] = C( 2, 4 );
            A[ 27 ] = C( 0, 5 );
            A[ 28 ] = C( 5, 5 );
            A[ 29 ] = C( 4, 5 );
            A[ 30 ] = C( 5, 5 );
            A[ 31 ] = C( 1, 5 );
            A[ 32 ] = C( 3, 5 );
            A[ 33 ] = C( 4, 5 );
            A[ 34 ] = C( 3, 5 );
            A[ 35 ] = C( 2, 5 );
            A[ 36 ] = C( 0, 1 );
            A[ 37 ] = C( 5, 1 );
            A[ 38 ] = C( 4, 1 );
            A[ 39 ] = C( 5, 1 );
            A[ 40 ] = C( 1, 1 );
            A[ 41 ] = C( 3, 1 );
            A[ 42 ] = C( 4, 1 );
            A[ 43 ] = C( 3, 1 );
            A[ 44 ] = C( 2, 1 );
            A[ 45 ] = C( 0, 3 );
            A[ 46 ] = C( 5, 3 );
            A[ 47 ] = C( 4, 3 );
            A[ 48 ] = C( 5, 3 );
            A[ 49 ] = C( 1, 3 );
            A[ 50 ] = C( 3, 3 );
            A[ 51 ] = C( 4, 3 );
            A[ 52 ] = C( 3, 3 );
            A[ 53 ] = C( 2, 3 );
            A[ 54 ] = C( 0, 4 );
            A[ 55 ] = C( 5, 4 );
            A[ 56 ] = C( 4, 4 );
            A[ 57 ] = C( 5, 4 );
            A[ 58 ] = C( 1, 4 );
            A[ 59 ] = C( 3, 4 );
            A[ 60 ] = C( 4, 4 );
            A[ 61 ] = C( 3, 4 );
            A[ 62 ] = C( 2, 4 );
            A[ 63 ] = C( 0, 3 );
            A[ 64 ] = C( 5, 3 );
            A[ 65 ] = C( 4, 3 );
            A[ 66 ] = C( 5, 3 );
            A[ 67 ] = C( 1, 3 );
            A[ 68 ] = C( 3, 3 );
            A[ 69 ] = C( 4, 3 );
            A[ 70 ] = C( 3, 3 );
            A[ 71 ] = C( 2, 3 );
            A[ 72 ] = C( 0, 2 );
            A[ 73 ] = C( 5, 2 );
            A[ 74 ] = C( 4, 2 );
            A[ 75 ] = C( 5, 2 );
            A[ 76 ] = C( 1, 2 );
            A[ 77 ] = C( 3, 2 );
            A[ 78 ] = C( 4, 2 );
            A[ 79 ] = C( 3, 2 );
            A[ 80 ] = C( 2, 2 );
#endif
        }
    } /* namespace tensor */
} /* namespace belfem */

#endif //BELFEM_FN_TR_MAT_TO_TEN_HPP
