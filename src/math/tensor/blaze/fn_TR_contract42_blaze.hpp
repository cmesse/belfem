//
// Created by Christian Messe on 17.12.20.
//
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
        contract42( const T * A, const T * B, T * C )
        {
            C[ 0 ] =  A[  0 ] * B[ 0 ]
                    + A[  9 ] * B[ 1 ]
                    + A[ 18 ] * B[ 2 ]
                    + A[ 27 ] * B[ 4 ]
                    + A[ 36 ] * B[ 5 ]
                    + A[ 45 ] * B[ 6 ]
                    + A[ 54 ] * B[ 8 ]
                    + A[ 63 ] * B[ 9 ]
                    + A[ 72 ] * B[ 10 ];

            C[ 1 ] =  A[  1 ] * B[ 0 ]
                    + A[ 10 ] * B[ 1 ]
                    + A[ 19 ] * B[ 2 ]
                    + A[ 28 ] * B[ 4 ]
                    + A[ 37 ] * B[ 5 ]
                    + A[ 46 ] * B[ 6 ]
                    + A[ 55 ] * B[ 8 ]
                    + A[ 64 ] * B[ 9 ]
                    + A[ 73 ] * B[ 10 ];

           C[ 2 ] =  A[  2 ] * B[ 0 ]
                    + A[ 11 ] * B[ 1 ]
                    + A[ 20 ] * B[ 2 ]
                    + A[ 29 ] * B[ 4 ]
                    + A[ 38 ] * B[ 5 ]
                    + A[ 47 ] * B[ 6 ]
                    + A[ 56 ] * B[ 8 ]
                    + A[ 65 ] * B[ 9 ]
                    + A[ 74 ] * B[ 10 ];

            C[ 4 ] =  A[  3 ] * B[ 0 ]
                    + A[ 12 ] * B[ 1 ]
                    + A[ 21 ] * B[ 2 ]
                    + A[ 30 ] * B[ 4 ]
                    + A[ 39 ] * B[ 5 ]
                    + A[ 48 ] * B[ 6 ]
                    + A[ 57 ] * B[ 8 ]
                    + A[ 66 ] * B[ 9 ]
                    + A[ 75 ] * B[ 10 ];

            C[ 5 ] =  A[  4 ] * B[ 0 ]
                    + A[ 13 ] * B[ 1 ]
                    + A[ 22 ] * B[ 2 ]
                    + A[ 31 ] * B[ 4 ]
                    + A[ 40 ] * B[ 5 ]
                    + A[ 49 ] * B[ 6 ]
                    + A[ 58 ] * B[ 8 ]
                    + A[ 67 ] * B[ 9 ]
                    + A[ 76 ] * B[ 10 ];

            C[ 6 ] =  A[  5 ] * B[ 0 ]
                    + A[ 14 ] * B[ 1 ]
                    + A[ 23 ] * B[ 2 ]
                    + A[ 32 ] * B[ 4 ]
                    + A[ 41 ] * B[ 5 ]
                    + A[ 50 ] * B[ 6 ]
                    + A[ 59 ] * B[ 8 ]
                    + A[ 68 ] * B[ 9 ]
                    + A[ 77 ] * B[ 10 ];

            C[ 8 ] =  A[  6 ] * B[ 0 ]
                    + A[ 15 ] * B[ 1 ]
                    + A[ 24 ] * B[ 2 ]
                    + A[ 33 ] * B[ 4 ]
                    + A[ 42 ] * B[ 5 ]
                    + A[ 51 ] * B[ 6 ]
                    + A[ 60 ] * B[ 8 ]
                    + A[ 69 ] * B[ 9 ]
                    + A[ 78 ] * B[ 10 ];

            C[ 9 ] =  A[  7 ] * B[ 0 ]
                    + A[ 16 ] * B[ 1 ]
                    + A[ 25 ] * B[ 2 ]
                    + A[ 34 ] * B[ 4 ]
                    + A[ 43 ] * B[ 5 ]
                    + A[ 52 ] * B[ 6 ]
                    + A[ 61 ] * B[ 8 ]
                    + A[ 70 ] * B[ 9 ]
                    + A[ 79 ] * B[ 10 ];

            C[ 10 ] =  A[  8 ] * B[ 0 ]
                     + A[ 17 ] * B[ 1 ]
                     + A[ 26 ] * B[ 2 ]
                     + A[ 35 ] * B[ 4 ]
                     + A[ 44 ] * B[ 5 ]
                     + A[ 53 ] * B[ 6 ]
                     + A[ 62 ] * B[ 8 ]
                     + A[ 71 ] * B[ 9 ]
                     + A[ 80 ] * B[ 10 ];


        }
//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */
#endif //BELFEM_FN_TR_CONTRACT42_BLAZE_HPP
