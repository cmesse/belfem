//
// Created by Christian Messe on 17.12.20.
//

#ifndef FN_BELFEM_TR_TO_MAT_HPP
#define FN_BELFEM_TR_TO_MAT_HPP
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
        ten_to_mat( T * C, const T * A )
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

//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */

#endif //FN_BELFEM_TR_TO_MAT_HPP
