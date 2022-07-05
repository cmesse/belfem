//
// Created by Christian Messe on 17.12.20.
//
#ifndef BELFEM_TR_CONTRACT42_ARMA_HPP
#define BELFEM_TR_CONTRACT42_ARMA_HPP

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

            C[ 0 ] =   A[  0 ] * B[ 0 ]
                     + A[  9 ] * B[ 1 ]
                     + A[ 18 ] * B[ 2 ]
                     + A[ 27 ] * B[ 3 ]
                     + A[ 36 ] * B[ 4 ]
                     + A[ 45 ] * B[ 5 ]
                     + A[ 54 ] * B[ 6 ]
                     + A[ 63 ] * B[ 7 ]
                     + A[ 72 ] * B[ 8 ];

            C[ 1 ] =   A[  1 ] * B[ 0 ]
                     + A[ 10 ] * B[ 1 ]
                     + A[ 19 ] * B[ 2 ]
                     + A[ 28 ] * B[ 3 ]
                     + A[ 37 ] * B[ 4 ]
                     + A[ 46 ] * B[ 5 ]
                     + A[ 55 ] * B[ 6 ]
                     + A[ 64 ] * B[ 7 ]
                     + A[ 73 ] * B[ 8 ];

            C[ 2 ] =   A[  2 ] * B[ 0 ]
                     + A[ 11 ] * B[ 1 ]
                     + A[ 20 ] * B[ 2 ]
                     + A[ 29 ] * B[ 3 ]
                     + A[ 38 ] * B[ 4 ]
                     + A[ 47 ] * B[ 5 ]
                     + A[ 56 ] * B[ 6 ]
                     + A[ 65 ] * B[ 7 ]
                     + A[ 74 ] * B[ 8 ];

            C[ 3 ] =   A[  3 ] * B[ 0 ]
                     + A[ 12 ] * B[ 1 ]
                     + A[ 21 ] * B[ 2 ]
                     + A[ 30 ] * B[ 3 ]
                     + A[ 39 ] * B[ 4 ]
                     + A[ 48 ] * B[ 5 ]
                     + A[ 57 ] * B[ 6 ]
                     + A[ 66 ] * B[ 7 ]
                     + A[ 75 ] * B[ 8 ];

            C[ 4 ] =   A[  4 ] * B[ 0 ]
                     + A[ 13 ] * B[ 1 ]
                     + A[ 22 ] * B[ 2 ]
                     + A[ 31 ] * B[ 3 ]
                     + A[ 40 ] * B[ 4 ]
                     + A[ 49 ] * B[ 5 ]
                     + A[ 58 ] * B[ 6 ]
                     + A[ 67 ] * B[ 7 ]
                     + A[ 76 ] * B[ 8 ];

            C[ 5 ] =   A[  5 ] * B[ 0 ]
                     + A[ 14 ] * B[ 1 ]
                     + A[ 23 ] * B[ 2 ]
                     + A[ 32 ] * B[ 3 ]
                     + A[ 41 ] * B[ 4 ]
                     + A[ 50 ] * B[ 5 ]
                     + A[ 59 ] * B[ 6 ]
                     + A[ 68 ] * B[ 7 ]
                     + A[ 77 ] * B[ 8 ];

            C[ 6 ] =   A[  6 ] * B[ 0 ]
                     + A[ 15 ] * B[ 1 ]
                     + A[ 24 ] * B[ 2 ]
                     + A[ 33 ] * B[ 3 ]
                     + A[ 42 ] * B[ 4 ]
                     + A[ 51 ] * B[ 5 ]
                     + A[ 60 ] * B[ 6 ]
                     + A[ 69 ] * B[ 7 ]
                     + A[ 78 ] * B[ 8 ];

            C[ 7 ] =   A[  7 ] * B[ 0 ]
                     + A[ 16 ] * B[ 1 ]
                     + A[ 25 ] * B[ 2 ]
                     + A[ 34 ] * B[ 3 ]
                     + A[ 43 ] * B[ 4 ]
                     + A[ 52 ] * B[ 5 ]
                     + A[ 61 ] * B[ 6 ]
                     + A[ 70 ] * B[ 7 ]
                     + A[ 79 ] * B[ 8 ];

            C[ 8 ] =   A[  8 ] * B[ 0 ]
                     + A[ 17 ] * B[ 1 ]
                     + A[ 26 ] * B[ 2 ]
                     + A[ 35 ] * B[ 3 ]
                     + A[ 44 ] * B[ 4 ]
                     + A[ 53 ] * B[ 5 ]
                     + A[ 62 ] * B[ 6 ]
                     + A[ 71 ] * B[ 7 ]
                     + A[ 80 ] * B[ 8 ];
        }
//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */
#endif //BELFEM_TR_CONTRACT42_ARMA_HPP
