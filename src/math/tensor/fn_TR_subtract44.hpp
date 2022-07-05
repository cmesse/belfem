//
// Created by Christian Messe on 17.12.20.
//

#ifndef BELFEM_FN_TR_SUBTRACT44_HPP
#define BELFEM_FN_TR_SUBTRACT44_HPP

namespace belfem
{
    namespace tensor
    {
//----------------------------------------------------------------------------
        /**
          * simple subtract
          * A = A - B
          */
        template < typename T > void
        subtract( T * A, const T * B )
        {
            A[  0 ] -= B[  0 ];
            A[  1 ] -= B[  1 ];
            A[  2 ] -= B[  2 ];
            A[  3 ] -= B[  3 ];
            A[  4 ] -= B[  4 ];
            A[  5 ] -= B[  5 ];
            A[  6 ] -= B[  6 ];
            A[  7 ] -= B[  7 ];
            A[  8 ] -= B[  8 ];
            A[  9 ] -= B[  9 ];
            A[ 10 ] -= B[ 10 ];
            A[ 11 ] -= B[ 11 ];
            A[ 12 ] -= B[ 12 ];
            A[ 13 ] -= B[ 13 ];
            A[ 14 ] -= B[ 14 ];
            A[ 15 ] -= B[ 15 ];
            A[ 16 ] -= B[ 16 ];
            A[ 17 ] -= B[ 17 ];
            A[ 18 ] -= B[ 18 ];
            A[ 19 ] -= B[ 19 ];
            A[ 20 ] -= B[ 20 ];
            A[ 21 ] -= B[ 21 ];
            A[ 22 ] -= B[ 22 ];
            A[ 23 ] -= B[ 23 ];
            A[ 24 ] -= B[ 24 ];
            A[ 25 ] -= B[ 25 ];
            A[ 26 ] -= B[ 26 ];
            A[ 27 ] -= B[ 27 ];
            A[ 28 ] -= B[ 28 ];
            A[ 29 ] -= B[ 29 ];
            A[ 30 ] -= B[ 30 ];
            A[ 31 ] -= B[ 31 ];
            A[ 32 ] -= B[ 32 ];
            A[ 33 ] -= B[ 33 ];
            A[ 34 ] -= B[ 34 ];
            A[ 35 ] -= B[ 35 ];
            A[ 36 ] -= B[ 36 ];
            A[ 37 ] -= B[ 37 ];
            A[ 38 ] -= B[ 38 ];
            A[ 39 ] -= B[ 39 ];
            A[ 40 ] -= B[ 40 ];
            A[ 41 ] -= B[ 41 ];
            A[ 42 ] -= B[ 42 ];
            A[ 43 ] -= B[ 43 ];
            A[ 44 ] -= B[ 44 ];
            A[ 45 ] -= B[ 45 ];
            A[ 46 ] -= B[ 46 ];
            A[ 47 ] -= B[ 47 ];
            A[ 48 ] -= B[ 48 ];
            A[ 49 ] -= B[ 49 ];
            A[ 50 ] -= B[ 50 ];
            A[ 51 ] -= B[ 51 ];
            A[ 52 ] -= B[ 52 ];
            A[ 53 ] -= B[ 53 ];
            A[ 54 ] -= B[ 54 ];
            A[ 55 ] -= B[ 55 ];
            A[ 56 ] -= B[ 56 ];
            A[ 57 ] -= B[ 57 ];
            A[ 58 ] -= B[ 58 ];
            A[ 59 ] -= B[ 59 ];
            A[ 60 ] -= B[ 60 ];
            A[ 61 ] -= B[ 61 ];
            A[ 62 ] -= B[ 62 ];
            A[ 63 ] -= B[ 63 ];
            A[ 64 ] -= B[ 64 ];
            A[ 65 ] -= B[ 65 ];
            A[ 66 ] -= B[ 66 ];
            A[ 67 ] -= B[ 67 ];
            A[ 68 ] -= B[ 68 ];
            A[ 69 ] -= B[ 69 ];
            A[ 70 ] -= B[ 70 ];
            A[ 71 ] -= B[ 71 ];
            A[ 72 ] -= B[ 72 ];
            A[ 73 ] -= B[ 73 ];
            A[ 74 ] -= B[ 74 ];
            A[ 75 ] -= B[ 75 ];
            A[ 76 ] -= B[ 76 ];
            A[ 77 ] -= B[ 77 ];
            A[ 78 ] -= B[ 78 ];
            A[ 79 ] -= B[ 79 ];
            A[ 80 ] -= B[ 80 ];
        }
//----------------------------------------------------------------------------

        /**
         * simple subtract
         *
         * A = B - C
         */
        template < typename T > void
        subtract( T * A, const T * B, const T * C )
        {
            A[  0 ] = B[  0 ] - C[  0 ];
            A[  1 ] = B[  1 ] - C[  1 ];
            A[  2 ] = B[  2 ] - C[  2 ];
            A[  3 ] = B[  3 ] - C[  3 ];
            A[  4 ] = B[  4 ] - C[  4 ];
            A[  5 ] = B[  5 ] - C[  5 ];
            A[  6 ] = B[  6 ] - C[  6 ];
            A[  7 ] = B[  7 ] - C[  7 ];
            A[  8 ] = B[  8 ] - C[  8 ];
            A[  9 ] = B[  9 ] - C[  9 ];
            A[ 10 ] = B[ 10 ] - C[ 10 ];
            A[ 11 ] = B[ 11 ] - C[ 11 ];
            A[ 12 ] = B[ 12 ] - C[ 12 ];
            A[ 13 ] = B[ 13 ] - C[ 13 ];
            A[ 14 ] = B[ 14 ] - C[ 14 ];
            A[ 15 ] = B[ 15 ] - C[ 15 ];
            A[ 16 ] = B[ 16 ] - C[ 16 ];
            A[ 17 ] = B[ 17 ] - C[ 17 ];
            A[ 18 ] = B[ 18 ] - C[ 18 ];
            A[ 19 ] = B[ 19 ] - C[ 19 ];
            A[ 20 ] = B[ 20 ] - C[ 20 ];
            A[ 21 ] = B[ 21 ] - C[ 21 ];
            A[ 22 ] = B[ 22 ] - C[ 22 ];
            A[ 23 ] = B[ 23 ] - C[ 23 ];
            A[ 24 ] = B[ 24 ] - C[ 24 ];
            A[ 25 ] = B[ 25 ] - C[ 25 ];
            A[ 26 ] = B[ 26 ] - C[ 26 ];
            A[ 27 ] = B[ 27 ] - C[ 27 ];
            A[ 28 ] = B[ 28 ] - C[ 28 ];
            A[ 29 ] = B[ 29 ] - C[ 29 ];
            A[ 30 ] = B[ 30 ] - C[ 30 ];
            A[ 31 ] = B[ 31 ] - C[ 31 ];
            A[ 32 ] = B[ 32 ] - C[ 32 ];
            A[ 33 ] = B[ 33 ] - C[ 33 ];
            A[ 34 ] = B[ 34 ] - C[ 34 ];
            A[ 35 ] = B[ 35 ] - C[ 35 ];
            A[ 36 ] = B[ 36 ] - C[ 36 ];
            A[ 37 ] = B[ 37 ] - C[ 37 ];
            A[ 38 ] = B[ 38 ] - C[ 38 ];
            A[ 39 ] = B[ 39 ] - C[ 39 ];
            A[ 40 ] = B[ 40 ] - C[ 40 ];
            A[ 41 ] = B[ 41 ] - C[ 41 ];
            A[ 42 ] = B[ 42 ] - C[ 42 ];
            A[ 43 ] = B[ 43 ] - C[ 43 ];
            A[ 44 ] = B[ 44 ] - C[ 44 ];
            A[ 45 ] = B[ 45 ] - C[ 45 ];
            A[ 46 ] = B[ 46 ] - C[ 46 ];
            A[ 47 ] = B[ 47 ] - C[ 47 ];
            A[ 48 ] = B[ 48 ] - C[ 48 ];
            A[ 49 ] = B[ 49 ] - C[ 49 ];
            A[ 50 ] = B[ 50 ] - C[ 50 ];
            A[ 51 ] = B[ 51 ] - C[ 51 ];
            A[ 52 ] = B[ 52 ] - C[ 52 ];
            A[ 53 ] = B[ 53 ] - C[ 53 ];
            A[ 54 ] = B[ 54 ] - C[ 54 ];
            A[ 55 ] = B[ 55 ] - C[ 55 ];
            A[ 56 ] = B[ 56 ] - C[ 56 ];
            A[ 57 ] = B[ 57 ] - C[ 57 ];
            A[ 58 ] = B[ 58 ] - C[ 58 ];
            A[ 59 ] = B[ 59 ] - C[ 59 ];
            A[ 60 ] = B[ 60 ] - C[ 60 ];
            A[ 61 ] = B[ 61 ] - C[ 61 ];
            A[ 62 ] = B[ 62 ] - C[ 62 ];
            A[ 63 ] = B[ 63 ] - C[ 63 ];
            A[ 64 ] = B[ 64 ] - C[ 64 ];
            A[ 65 ] = B[ 65 ] - C[ 65 ];
            A[ 66 ] = B[ 66 ] - C[ 66 ];
            A[ 67 ] = B[ 67 ] - C[ 67 ];
            A[ 68 ] = B[ 68 ] - C[ 68 ];
            A[ 69 ] = B[ 69 ] - C[ 69 ];
            A[ 70 ] = B[ 70 ] - C[ 70 ];
            A[ 71 ] = B[ 71 ] - C[ 71 ];
            A[ 72 ] = B[ 72 ] - C[ 72 ];
            A[ 73 ] = B[ 73 ] - C[ 73 ];
            A[ 74 ] = B[ 74 ] - C[ 74 ];
            A[ 75 ] = B[ 75 ] - C[ 75 ];
            A[ 76 ] = B[ 76 ] - C[ 76 ];
            A[ 77 ] = B[ 77 ] - C[ 77 ];
            A[ 78 ] = B[ 78 ] - C[ 78 ];
            A[ 79 ] = B[ 79 ] - C[ 79 ];
            A[ 80 ] = B[ 80 ] - C[ 80 ];
        }
//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */

#endif //BELFEM_FN_TR_SUBTRACT44_HPP
