/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_TR_FILL_HPP
#define BELFEM_FN_TR_FILL_HPP


namespace belfem
{
    namespace tensor
    {
//---------------------------------------------------------------------------

        /*
         * fills a tensor in the fashion of
         *
         * A_ijkl = a * δ_ij * δ_kl + b * (δ_ik * δ_jl + δ_il * δ_jk - 2/3 * δ_ij * δ_kl )
          */
        template < typename T >
        void
        fill( T * A, const T a, const T b )
        {
            T  z = ( T ) 0 ;
            T  c = a + b*4.0/3.0;
            T  d = a - b*2.0/3.0 ;

            A[  0 ] = c;
            A[  1 ] = z;
            A[  2 ] = z;
            A[  3 ] = z;
            A[  4 ] = d;
            A[  5 ] = z;
            A[  6 ] = z;
            A[  7 ] = z;
            A[  8 ] = d;
            A[  9 ] = z;
            A[ 10 ] = b;
            A[ 11 ] = z;
            A[ 12 ] = b;
            A[ 13 ] = z;
            A[ 14 ] = z;
            A[ 15 ] = z;
            A[ 16 ] = z;
            A[ 17 ] = z;
            A[ 18 ] = z;
            A[ 19 ] = z;
            A[ 20 ] = b;
            A[ 21 ] = z;
            A[ 22 ] = z;
            A[ 23 ] = z;
            A[ 24 ] = b;
            A[ 25 ] = z;
            A[ 26 ] = z;
            A[ 27 ] = z;
            A[ 28 ] = b;
            A[ 29 ] = z;
            A[ 30 ] = b;
            A[ 31 ] = z;
            A[ 32 ] = z;
            A[ 33 ] = z;
            A[ 34 ] = z;
            A[ 35 ] = z;
            A[ 36 ] = d;
            A[ 37 ] = z;
            A[ 38 ] = z;
            A[ 39 ] = z;
            A[ 40 ] = c;
            A[ 41 ] = z;
            A[ 42 ] = z;
            A[ 43 ] = z;
            A[ 44 ] = d;
            A[ 45 ] = z;
            A[ 46 ] = z;
            A[ 47 ] = z;
            A[ 48 ] = z;
            A[ 49 ] = z;
            A[ 50 ] = b;
            A[ 51 ] = z;
            A[ 52 ] = b;
            A[ 53 ] = z;
            A[ 54 ] = z;
            A[ 55 ] = z;
            A[ 56 ] = b;
            A[ 57 ] = z;
            A[ 58 ] = z;
            A[ 59 ] = z;
            A[ 60 ] = b;
            A[ 61 ] = z;
            A[ 62 ] = z;
            A[ 63 ] = z;
            A[ 64 ] = z;
            A[ 65 ] = z;
            A[ 66 ] = z;
            A[ 67 ] = z;
            A[ 68 ] = b;
            A[ 69 ] = z;
            A[ 70 ] = b;
            A[ 71 ] = z;
            A[ 72 ] = d;
            A[ 73 ] = z;
            A[ 74 ] = z;
            A[ 75 ] = z;
            A[ 76 ] = d;
            A[ 77 ] = z;
            A[ 78 ] = z;
            A[ 79 ] = z;
            A[ 80 ] = c;

        }
    }
}
#endif //BELFEM_FN_TR_FILL_HPP