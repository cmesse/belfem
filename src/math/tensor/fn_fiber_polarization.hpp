//
// Created by Christian Messe on 23.12.20.
//

#ifndef BELFEM_FN_ESHELBY_FIBER_HPP
#define BELFEM_FN_ESHELBY_FIBER_HPP
#include "assert.h"
#include "cl_Tensor.hpp"

namespace belfem
{
    /**
     * create the eshelby tensor for fibers
     * see doi: 10.1016/S0020-7683(02)00369-4.
     */
    template< typename T >
    inline void
    fiber_polarization( const Matrix <T> & aC,
                        Tensor <T> & aP )
    {
        BELFEM_ASSERT( aC.n_cols() == 6 && aC.n_rows() == 6,
                      "elasticity matrix must be 6x6" );

        BELFEM_ASSERT( aP.is_3333(),
                      "polarization tensor must be 3x3x3x3" );

        // get data containers

        T * P = aP.data();

        const T & C11 = aC( 2,2 );
        const T & C12 = aC( 1,2 );
        const T & C22 = aC( 1,1 );
        //const T & C33 = aC( 0,0 );

        const T & C44 = aC( 5,5 );
        const T & C55 = aC( 4,4 );
        const T & C66 = aC( 3,3 );

        T g = std::sqrt( ( C11 * C22 - C12 * C12
                           + 2. * C66 * ( std::sqrt( C11 * C22 ) - C12 ) ) / ( C22 * C66 ) );

        T b = std::sqrt( C55 / C44 );
        T h = std::sqrt( C11 / C22 ) ;
        T s = g * g - h - h ;

        T n0 = 2. * ( h + b * ( h + g ) ) + ( g + b ) + b * b * ( g + 1.) + s * ( b + 1.) + h * g ;
        T n2 = g + b + 1. ;
        T n4 = b * h + h + b * g ;
        T n6 = 2. * h * b * (1. + g + b ) + ( b * s + b * b * g + h * g ) + b * b * ( s + h * g ) + h * h * ( b + 1. );

        T d1 = g * ( h + b * g + b * b ) * ( 1. + g + h ) * ( b + 1. );
        T d2 = h * b ;
        T r  = 1.0 /  ( d1 * C22 * C44 * C66 );

        T P33 = r * ( C55 * C66 * n0 / d2 + ( C22 * C55 + C44 * C66 ) * n2 + C22 * C44 * n4 );
        T P12 = 0.0 ;
        T P13 = 0.0 ;

        T P22 = r * ( C11 * C55 * n2 + ( C55 * C66 + C11 * C44 ) * n4 + C44 * C66 * n6 / d2 ) ;
        T P23 = - r * ( C55 * ( C12 + C66 ) * n2 + C44 * ( C12 + C66 ) * n4 ) ;
        T P11 = 0.0 ;

        T P44 = r * ( C11 * C55 * n0 / d2 + ( C11 * C44 - 2.0 * C12 * C55 ) * n2
                      + ( C22 * C55 - 2. * C12 * C44 ) * n4 + C22 * C44 * n6 );

        T P55 = r * ( C11 * C66 * n0 / d2
                      + ( C11 * C22 + C66 * C66 - ( C12 + C66 ) * ( C12 + C66 ) ) * n2 + C22 * C66 * n4 );

        T P66 = r * ( C11 * C66 * n2 + ( C11 * C22 + C66 * C66 - ( C12 + C66 ) * ( C12 + C66 ) ) * n4
                      + C22 * C66 * n6 );

        // populate tensor
        P[ 0 ] = P11; // 1111
        P[ 1 ] = 0.0; // 2111
        P[ 2 ] = 0.0; // 3111

        P[ 3 ] = 0.0; // 1211
        P[ 4 ] = P12; // 2211
        P[ 5 ] = 0.0; // 3211

        P[ 6 ] = 0.0; // 1311
        P[ 7 ] = 0.0; // 2311
        P[ 8 ] = P13; // 3311

        P[ 9 ] = 0.0; // 1121
        P[ 10 ] = P66; // 2121
        P[ 11 ] = 0.0; // 3121

        P[ 12 ] = P66; // 1221
        P[ 13 ] = 0.0; // 2221
        P[ 14 ] = 0.0; // 3221

        P[ 15 ] = 0.0; // 1321
        P[ 16 ] = 0.0; // 2321
        P[ 17 ] = 0.0; // 3321

        P[ 18 ] = 0.0; // 1131
        P[ 19 ] = 0.0; // 2131
        P[ 20 ] = P55; // 3131

        P[ 21 ] = 0.0; // 1231
        P[ 22 ] = 0.0; // 2231
        P[ 23 ] = 0.0; // 3231

        P[ 24 ] = P55; // 1331
        P[ 25 ] = 0.0; // 2331
        P[ 26 ] = 0.0; // 3331

        P[ 27 ] = 0.0; // 1112
        P[ 28 ] = P66; // 2112
        P[ 29 ] = 0.0; // 3112

        P[ 30 ] = P66; // 1212
        P[ 31 ] = 0.0; // 2212
        P[ 32 ] = 0.0; // 3212

        P[ 33 ] = 0.0; // 1312
        P[ 34 ] = 0.0; // 2312
        P[ 35 ] = 0.0; // 3312

        P[ 36 ] = P12; // 1122
        P[ 37 ] = 0.0; // 2122
        P[ 38 ] = 0.0; // 3122

        P[ 39 ] = 0.0; // 1222
        P[ 40 ] = P22; // 2222
        P[ 41 ] = 0.0; // 3222

        P[ 42 ] = 0.0; // 1322
        P[ 43 ] = 0.0; // 2322
        P[ 44 ] = P23; // 3322

        P[ 45 ] = 0.0; // 1132
        P[ 46 ] = 0.0; // 2132
        P[ 47 ] = 0.0; // 3132

        P[ 48 ] = 0.0; // 1232
        P[ 49 ] = 0.0; // 2232
        P[ 50 ] = P44; // 3232

        P[ 51 ] = 0.0; // 1332
        P[ 52 ] = P44; // 2332
        P[ 53 ] = 0.0; // 3332

        P[ 54 ] = 0.0; // 1113
        P[ 55 ] = 0.0; // 2113
        P[ 56 ] = P55; // 3113

        P[ 57 ] = 0.0; // 1213
        P[ 58 ] = 0.0; // 2213
        P[ 59 ] = 0.0; // 3213

        P[ 60 ] = P55; // 1313
        P[ 61 ] = 0.0; // 2313
        P[ 62 ] = 0.0; // 3313

        P[ 63 ] = 0.0; // 1123
        P[ 64 ] = 0.0; // 2123
        P[ 65 ] = 0.0; // 3123

        P[ 66 ] = 0.0; // 1223
        P[ 67 ] = 0.0; // 2223
        P[ 68 ] = P44; // 3223

        P[ 69 ] = 0.0; // 1323
        P[ 70 ] = P44; // 2323
        P[ 71 ] = 0.0; // 3323

        P[ 72 ] = P13; // 1133
        P[ 73 ] = 0.0; // 2133
        P[ 74 ] = 0.0; // 3133

        P[ 75 ] = 0.0; // 1233
        P[ 76 ] = P23; // 2233
        P[ 77 ] = 0.0; // 3233

        P[ 78 ] = 0.0; // 1333
        P[ 79 ] = 0.0; // 2333
        P[ 80 ] = P33; // 3333
    }
}
#endif //BELFEM_FN_ESHELBY_FIBER_HPP
