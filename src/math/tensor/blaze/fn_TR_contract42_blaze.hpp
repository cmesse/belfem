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

#ifndef BELFEM_TR_CONTRACT42_BLAZE_HPP
#define BELFEM_TR_CONTRACT42_BLAZE_HPP

#include "cl_Matrix.hpp"
namespace belfem
{
    namespace tensor
    {
//----------------------------------------------------------------------------

        /**
         * tensor contraction A_ijkl * B_kl = C_ij
         */
        template < typename T >
        inline void
        contract42( const T * A, const Matrix< T > & B, Matrix< T > & C )
        {
            C( 0, 0 ) =   A[  0 ] * B( 0, 0 )
                         + A[  9 ] * B( 1, 0 )
                         + A[ 18 ] * B( 2, 0 )
                         + A[ 27 ] * B( 0, 1 )
                         + A[ 36 ] * B( 1, 1 )
                         + A[ 45 ] * B( 2, 1 )
                         + A[ 54 ] * B( 0, 2 )
                         + A[ 63 ] * B( 1, 2 )
                         + A[ 72 ] * B( 2, 2 );

            C( 1, 0 ) =    A[  1 ] * B( 0, 0 )
                         + A[ 10 ] * B( 1, 0 )
                         + A[ 19 ] * B( 2, 0 )
                         + A[ 28 ] * B( 0, 1 )
                         + A[ 37 ] * B( 1, 1 )
                         + A[ 46 ] * B( 2, 1 )
                         + A[ 55 ] * B( 0, 2 )
                         + A[ 64 ] * B( 1, 2 )
                         + A[ 73 ] * B( 2, 2 );

            C( 2, 0 ) =    A[  2 ] * B( 0, 0 )
                         + A[ 11 ] * B( 1, 0 )
                         + A[ 20 ] * B( 2, 0 )
                         + A[ 29 ] * B( 0, 1 )
                         + A[ 38 ] * B( 1, 1 )
                         + A[ 47 ] * B( 2, 1 )
                         + A[ 56 ] * B( 0, 2 )
                         + A[ 65 ] * B( 1, 2 )
                         + A[ 74 ] * B( 2, 2 );

            C( 0, 1 ) =    A[  3 ] * B( 0, 0 )
                         + A[ 12 ] * B( 1, 0 )
                         + A[ 21 ] * B( 2, 0 )
                         + A[ 30 ] * B( 0, 1 )
                         + A[ 39 ] * B( 1, 1 )
                         + A[ 48 ] * B( 2, 1 )
                         + A[ 57 ] * B( 0, 2 )
                         + A[ 66 ] * B( 1, 2 )
                         + A[ 75 ] * B( 2, 2 );

            C( 1, 1 ) =    A[  4 ] * B( 0, 0 )
                         + A[ 13 ] * B( 1, 0 )
                         + A[ 22 ] * B( 2, 0 )
                         + A[ 31 ] * B( 0, 1 )
                         + A[ 40 ] * B( 1, 1 )
                         + A[ 49 ] * B( 2, 1 )
                         + A[ 58 ] * B( 0, 2 )
                         + A[ 67 ] * B( 1, 2 )
                         + A[ 76 ] * B( 2, 2 );

            C( 2, 1 )=     A[  5 ] * B( 0, 0 )
                         + A[ 14 ] * B( 1, 0 )
                         + A[ 23 ] * B( 2, 0 )
                         + A[ 32 ] * B( 0, 1 )
                         + A[ 41 ] * B( 1, 1 )
                         + A[ 50 ] * B( 2, 1 )
                         + A[ 59 ] * B( 0, 2 )
                         + A[ 68 ] * B( 1, 2 )
                         + A[ 77 ] * B( 2, 2 );

            C( 0, 2 ) =    A[  6 ] * B( 0, 0 )
                         + A[ 15 ] * B( 1, 0 )
                         + A[ 24 ] * B( 2, 0 )
                         + A[ 33 ] * B( 0, 1 )
                         + A[ 42 ] * B( 1, 1 )
                         + A[ 51 ] * B( 2, 1 )
                         + A[ 60 ] * B( 0, 2 )
                         + A[ 69 ] * B( 1, 2 )
                         + A[ 78 ] * B( 2, 2 );

            C( 1, 2 ) =    A[  7 ] * B( 0, 0 )
                         + A[ 16 ] * B( 1, 0 )
                         + A[ 25 ] * B( 2, 0 )
                         + A[ 34 ] * B( 0, 1 )
                         + A[ 43 ] * B( 1, 1 )
                         + A[ 52 ] * B( 2, 1 )
                         + A[ 61 ] * B( 0, 2 )
                         + A[ 70 ] * B( 1, 2 )
                         + A[ 79 ] * B( 2, 2 );

            C( 2, 2 ) =    A[  8 ] * B( 0, 0 )
                         + A[ 17 ] * B( 1, 0 )
                         + A[ 26 ] * B( 2, 0 )
                         + A[ 35 ] * B( 0, 1 )
                         + A[ 44 ] * B( 1, 1 )
                         + A[ 53 ] * B( 2, 1 )
                         + A[ 62 ] * B( 0, 2 )
                         + A[ 71 ] * B( 1, 2 )
                         + A[ 80 ] * B( 2, 2 );
        }
//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */
#endif //BELFEM_FN_TR_CONTRACT42_BLAZE_HPP
