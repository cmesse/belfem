//
// Created by Christian Messe on 18.12.20.
//

#ifndef BELFEM_INVERT_SYMMETRIC_HPP
#define BELFEM_INVERT_SYMMETRIC_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "fn_LAPACK_getrf.hpp"
#include "fn_LAPACK_getri.hpp"
#include "fn_LAPACK_gemm.hpp"


#include "cl_Tensor.hpp"
#include "cl_Vector.hpp"



namespace belfem
{
//----------------------------------------------------------------------------

    namespace tensor
    {
        template < typename T >
        inline void
        invert_symmetric(
                Tensor< T >    & aTensor,
                Vector< T >    & aWork,
                Vector< int  > & aPivot )
        {
            BELFEM_ASSERT( aTensor.is_3333(), "The tensor must be a 3x3x3x3 tensor" );

            BELFEM_ASSERT( aWork.length() >= 72,
                          "The work vector must have a length of at least 72." );
            
            BELFEM_ASSERT( aPivot.length() >= 36,
                          "The pivot vector must have a length of at least 36" );
            
            // pointer to tensor data
            T * tData = aTensor.data();

            // pointer to matrix data
            T * tWorkA = aWork.data();

            // pointer to Indentity data
            T * tWorkB = aWork.data() + 36 ;

            // pointer to Indentity data, can use tensor here
            T * tWorkI = aTensor.data();

            // copy data into WorkA vector
            tWorkA[  0 ] = tData[  0 ];
            tWorkA[  1 ] = tData[  1 ];
            tWorkA[  2 ] = tData[  2 ];
            tWorkA[  3 ] = tData[  4 ];
            tWorkA[  4 ] = tData[  5 ];
            tWorkA[  5 ] = tData[  8 ];
            tWorkA[  6 ] = tData[  9 ] + tData[  9 ];
            tWorkA[  7 ] = tData[ 10 ] + tData[ 10 ];
            tWorkA[  8 ] = tData[ 11 ] + tData[ 11 ];
            tWorkA[  9 ] = tData[ 13 ] + tData[ 13 ];
            tWorkA[ 10 ] = tData[ 14 ] + tData[ 14 ];
            tWorkA[ 11 ] = tData[ 17 ] + tData[ 17 ];
            tWorkA[ 12 ] = tData[ 18 ] + tData[ 18 ];
            tWorkA[ 13 ] = tData[ 19 ] + tData[ 19 ];
            tWorkA[ 14 ] = tData[ 20 ] + tData[ 20 ];
            tWorkA[ 15 ] = tData[ 22 ] + tData[ 22 ];
            tWorkA[ 16 ] = tData[ 23 ] + tData[ 23 ];
            tWorkA[ 17 ] = tData[ 26 ] + tData[ 26 ];
            tWorkA[ 18 ] = tData[ 36 ];
            tWorkA[ 19 ] = tData[ 37 ];
            tWorkA[ 20 ] = tData[ 38 ];
            tWorkA[ 21 ] = tData[ 40 ];
            tWorkA[ 22 ] = tData[ 41 ];
            tWorkA[ 23 ] = tData[ 44 ];
            tWorkA[ 24 ] = tData[ 45 ] + tData[ 45 ];
            tWorkA[ 25 ] = tData[ 46 ] + tData[ 46 ];
            tWorkA[ 26 ] = tData[ 47 ] + tData[ 47 ];
            tWorkA[ 27 ] = tData[ 49 ] + tData[ 49 ];
            tWorkA[ 28 ] = tData[ 50 ] + tData[ 50 ];
            tWorkA[ 29 ] = tData[ 53 ] + tData[ 53 ];
            tWorkA[ 30 ] = tData[ 72 ];
            tWorkA[ 31 ] = tData[ 73 ];
            tWorkA[ 32 ] = tData[ 74 ];
            tWorkA[ 33 ] = tData[ 76 ];
            tWorkA[ 34 ] = tData[ 77 ];
            tWorkA[ 35 ] = tData[ 80 ];

            // help values
            int    tInfo   = 0 ;

            // create an LU factorization for A
            lapack::getrf( 6, 6, tWorkA, 6, aPivot.data(), tInfo );

            BELFEM_ASSERT( tInfo == 0, "LAPACK getrf has thrown error %i", tInfo );

            // create the inverse of A using the LU factorization
            lapack::getri( 6, tWorkA, 6, aPivot.data(), tWorkB, 36, tInfo );

            BELFEM_ASSERT( tInfo == 0, "LAPACK getri has thrown error %i", tInfo );

            // initialize identity matrix
            tWorkI[  0 ] = 1.0 ;
            tWorkI[  1 ] = 0.0 ;
            tWorkI[  2 ] = 0.0 ;
            tWorkI[  3 ] = 0.0 ;
            tWorkI[  4 ] = 0.0 ;
            tWorkI[  5 ] = 0.0 ;

            tWorkI[  6 ] = 0.0 ;
            tWorkI[  7 ] = 0.5 ;
            tWorkI[  8 ] = 0.0 ;
            tWorkI[  9 ] = 0.0 ;
            tWorkI[ 10 ] = 0.0 ;
            tWorkI[ 11 ] = 0.0 ;

            tWorkI[ 12 ] = 0.0 ;
            tWorkI[ 13 ] = 0.0 ;
            tWorkI[ 14 ] = 0.5 ;
            tWorkI[ 15 ] = 0.0 ;
            tWorkI[ 16 ] = 0.0 ;
            tWorkI[ 17 ] = 0.0 ;

            tWorkI[ 18 ] = 0.0 ;
            tWorkI[ 19 ] = 0.0 ;
            tWorkI[ 20 ] = 0.0 ;
            tWorkI[ 21 ] = 1.0 ;
            tWorkI[ 22 ] = 0.0 ;
            tWorkI[ 23 ] = 0.0 ;

            tWorkI[ 24 ] = 0.0 ;
            tWorkI[ 25 ] = 0.0 ;
            tWorkI[ 26 ] = 0.0 ;
            tWorkI[ 27 ] = 0.0 ;
            tWorkI[ 28 ] = 0.5 ;
            tWorkI[ 29 ] = 0.0 ;

            tWorkI[ 30 ] = 0.0 ;
            tWorkI[ 31 ] = 0.0 ;
            tWorkI[ 32 ] = 0.0 ;
            tWorkI[ 33 ] = 0.0 ;
            tWorkI[ 34 ] = 0.0 ;
            tWorkI[ 35 ] = 1.0 ;

            // call lapack for matrix-matrix multiplication
            lapack::gemm( 'N', 'N',
                  6, 6, 6,
                  1.0, tWorkA, 6,
                  tWorkI, 6, 0.0,
                  tWorkB, 6 );

            // write data back
            tData[  0 ] = tWorkB[  0 ];
            tData[  1 ] = tWorkB[  1 ];
            tData[  2 ] = tWorkB[  2 ];
            tData[  3 ] = tWorkB[  1 ];
            tData[  4 ] = tWorkB[  3 ];
            tData[  5 ] = tWorkB[  4 ];
            tData[  6 ] = tWorkB[  2 ];
            tData[  7 ] = tWorkB[  4 ];
            tData[  8 ] = tWorkB[  5 ];
            tData[  9 ] = tWorkB[  6 ];
            tData[ 10 ] = tWorkB[  7 ];
            tData[ 11 ] = tWorkB[  8 ];
            tData[ 12 ] = tWorkB[  7 ];
            tData[ 13 ] = tWorkB[  9 ];
            tData[ 14 ] = tWorkB[ 10 ];
            tData[ 15 ] = tWorkB[  8 ];
            tData[ 16 ] = tWorkB[ 10 ];
            tData[ 17 ] = tWorkB[ 11 ];
            tData[ 18 ] = tWorkB[ 12 ];
            tData[ 19 ] = tWorkB[ 13 ];
            tData[ 20 ] = tWorkB[ 14 ];
            tData[ 21 ] = tWorkB[ 13 ];
            tData[ 22 ] = tWorkB[ 15 ];
            tData[ 23 ] = tWorkB[ 16 ];
            tData[ 24 ] = tWorkB[ 14 ];
            tData[ 25 ] = tWorkB[ 16 ];
            tData[ 26 ] = tWorkB[ 17 ];
            tData[ 27 ] = tWorkB[  6 ];
            tData[ 28 ] = tWorkB[  7 ];
            tData[ 29 ] = tWorkB[  8 ];
            tData[ 30 ] = tWorkB[  7 ];
            tData[ 31 ] = tWorkB[  9 ];
            tData[ 32 ] = tWorkB[ 10 ];
            tData[ 33 ] = tWorkB[  8 ];
            tData[ 34 ] = tWorkB[ 10 ];
            tData[ 35 ] = tWorkB[ 11 ];
            tData[ 36 ] = tWorkB[ 18 ];
            tData[ 37 ] = tWorkB[ 19 ];
            tData[ 38 ] = tWorkB[ 20 ];
            tData[ 39 ] = tWorkB[ 19 ];
            tData[ 40 ] = tWorkB[ 21 ];
            tData[ 41 ] = tWorkB[ 22 ];
            tData[ 42 ] = tWorkB[ 20 ];
            tData[ 43 ] = tWorkB[ 22 ];
            tData[ 44 ] = tWorkB[ 23 ];
            tData[ 45 ] = tWorkB[ 24 ];
            tData[ 46 ] = tWorkB[ 25 ];
            tData[ 47 ] = tWorkB[ 26 ];
            tData[ 48 ] = tWorkB[ 25 ];
            tData[ 49 ] = tWorkB[ 27 ];
            tData[ 50 ] = tWorkB[ 28 ];
            tData[ 51 ] = tWorkB[ 26 ];
            tData[ 52 ] = tWorkB[ 28 ];
            tData[ 53 ] = tWorkB[ 29 ];
            tData[ 54 ] = tWorkB[ 12 ];
            tData[ 55 ] = tWorkB[ 13 ];
            tData[ 56 ] = tWorkB[ 14 ];
            tData[ 57 ] = tWorkB[ 13 ];
            tData[ 58 ] = tWorkB[ 15 ];
            tData[ 59 ] = tWorkB[ 16 ];
            tData[ 60 ] = tWorkB[ 14 ];
            tData[ 61 ] = tWorkB[ 16 ];
            tData[ 62 ] = tWorkB[ 17 ];
            tData[ 63 ] = tWorkB[ 24 ];
            tData[ 64 ] = tWorkB[ 25 ];
            tData[ 65 ] = tWorkB[ 26 ];
            tData[ 66 ] = tWorkB[ 25 ];
            tData[ 67 ] = tWorkB[ 27 ];
            tData[ 68 ] = tWorkB[ 28 ];
            tData[ 69 ] = tWorkB[ 26 ];
            tData[ 70 ] = tWorkB[ 28 ];
            tData[ 71 ] = tWorkB[ 29 ];
            tData[ 72 ] = tWorkB[ 30 ];
            tData[ 73 ] = tWorkB[ 31 ];
            tData[ 74 ] = tWorkB[ 32 ];
            tData[ 75 ] = tWorkB[ 31 ];
            tData[ 76 ] = tWorkB[ 33 ];
            tData[ 77 ] = tWorkB[ 34 ];
            tData[ 78 ] = tWorkB[ 32 ];
            tData[ 79 ] = tWorkB[ 34 ];
            tData[ 80 ] = tWorkB[ 35 ];

        }
    }
//----------------------------------------------------------------------------

}

#endif //BELFEM_INVERT_SYMMETRIC_HPP
