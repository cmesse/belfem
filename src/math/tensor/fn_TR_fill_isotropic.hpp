//
// Created by Christian Messe on 19.12.20.
//

#ifndef BELFEM_FN_TR_FILL_ISOTROPIC_HPP
#define BELFEM_FN_TR_FILL_ISOTROPIC_HPP

namespace belfem
{
    namespace tensor
    {
//---------------------------------------------------------------------------

        template < typename T >
        inline void
        fill_isotropic( T * A, const T & kappa, const T & gammadivthree )
        {

            T alpha = gammadivthree + gammadivthree ;
            T beta = alpha + gammadivthree ;
            T delta = alpha + alpha ;
            T zero = ( T ) 0 ;

            A[  0 ] = kappa +  delta ;
            A[  1 ] =           zero ;
            A[  2 ] =           zero ;
            A[  3 ] =           zero ;
            A[  4 ] = kappa -  alpha ;
            A[  5 ] =           zero ;
            A[  6 ] =           zero ;
            A[  7 ] =           zero ;
            A[  8 ] = kappa -  alpha ;
            A[  9 ] =           zero ;
            A[ 10 ] =           beta ;
            A[ 11 ] =           zero ;
            A[ 12 ] =           beta ;
            A[ 13 ] =           zero ;
            A[ 14 ] =           zero ;
            A[ 15 ] =           zero ;
            A[ 16 ] =           zero ;
            A[ 17 ] =           zero ;
            A[ 18 ] =           zero ;
            A[ 19 ] =           zero ;
            A[ 20 ] =           beta ;
            A[ 21 ] =           zero ;
            A[ 22 ] =           zero ;
            A[ 23 ] =           zero ;
            A[ 24 ] =           beta ;
            A[ 25 ] =           zero ;
            A[ 26 ] =           zero ;
            A[ 27 ] =           zero ;
            A[ 28 ] =           beta ;
            A[ 29 ] =           zero ;
            A[ 30 ] =           beta ;
            A[ 31 ] =           zero ;
            A[ 32 ] =           zero ;
            A[ 33 ] =           zero ;
            A[ 34 ] =           zero ;
            A[ 35 ] =           zero ;
            A[ 36 ] = kappa -  alpha ;
            A[ 37 ] =           zero ;
            A[ 38 ] =           zero ;
            A[ 39 ] =           zero ;
            A[ 40 ] = kappa +  delta ;
            A[ 41 ] =           zero ;
            A[ 42 ] =           zero ;
            A[ 43 ] =           zero ;
            A[ 44 ] = kappa -  alpha ;
            A[ 45 ] =           zero ;
            A[ 46 ] =           zero ;
            A[ 47 ] =           zero ;
            A[ 48 ] =           zero ;
            A[ 49 ] =           zero ;
            A[ 50 ] =           beta ;
            A[ 51 ] =           zero ;
            A[ 52 ] =           beta ;
            A[ 53 ] =           zero ;
            A[ 54 ] =           zero ;
            A[ 55 ] =           zero ;
            A[ 56 ] =           beta ;
            A[ 57 ] =           zero ;
            A[ 58 ] =           zero ;
            A[ 59 ] =           zero ;
            A[ 60 ] =           beta ;
            A[ 61 ] =           zero ;
            A[ 62 ] =           zero ;
            A[ 63 ] =           zero ;
            A[ 64 ] =           zero ;
            A[ 65 ] =           zero ;
            A[ 66 ] =           zero ;
            A[ 67 ] =           zero ;
            A[ 68 ] =           beta ;
            A[ 69 ] =           zero ;
            A[ 70 ] =           beta ;
            A[ 71 ] =           zero ;
            A[ 72 ] = kappa -  alpha ;
            A[ 73 ] =           zero ;
            A[ 74 ] =           zero ;
            A[ 75 ] =           zero ;
            A[ 76 ] = kappa -  alpha ;
            A[ 77 ] =           zero ;
            A[ 78 ] =           zero ;
            A[ 79 ] =           zero ;
            A[ 80 ] = kappa +  delta ;
        }
        
//---------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_TR_FILL_ISOTROPIC_HPP
