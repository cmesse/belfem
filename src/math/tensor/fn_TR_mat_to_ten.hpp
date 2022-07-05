//
// Created by Christian Messe on 17.12.20.
//

#ifndef BELFEM_FN_TR_MAT_TO_TEN_HPP
#define BELFEM_FN_TR_MAT_TO_TEN_HPP
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
        mat_to_ten( T * A, const T * C  )
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
    } /* namespace tensor */
} /* namespace belfem */

#endif //BELFEM_FN_TR_MAT_TO_TEN_HPP
