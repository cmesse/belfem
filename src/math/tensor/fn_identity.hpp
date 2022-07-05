//
// Created by Christian Messe on 18.12.20.
//

#ifndef BELFEM_FN_IDENTITY_HPP
#define BELFEM_FN_IDENTITY_HPP

#include "cl_Tensor.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace tensor
    {
//---------------------------------------------------------------------------

        /**
         * populate the symmetric identity tensor
         *
         * I_ijkl = 0.5 * ( d_ik * d_jl + d_il * d_jk )
         * @param aTensor
         */
        template< typename T >
        void
        identity( Tensor< T > & aTensor )
        {
            BELFEM_ASSERT( aTensor.is_3333(),
                          "Tensor must be of size 3x3x3x3" );

            // get pointer to tensor
            T * tData = aTensor.data();

            // populate
            tData[  0 ] = 1.0 ;
            tData[  1 ] = 0.0 ;
            tData[  2 ] = 0.0 ;
            tData[  3 ] = 0.0 ;
            tData[  4 ] = 0.0 ;
            tData[  5 ] = 0.0 ;
            tData[  6 ] = 0.0 ;
            tData[  7 ] = 0.0 ;
            tData[  8 ] = 0.0 ;
            tData[  9 ] = 0.0 ;
            tData[ 10 ] = 0.5 ;
            tData[ 11 ] = 0.0 ;
            tData[ 12 ] = 0.5 ;
            tData[ 13 ] = 0.0 ;
            tData[ 14 ] = 0.0 ;
            tData[ 15 ] = 0.0 ;
            tData[ 16 ] = 0.0 ;
            tData[ 17 ] = 0.0 ;
            tData[ 18 ] = 0.0 ;
            tData[ 19 ] = 0.0 ;
            tData[ 20 ] = 0.5 ;
            tData[ 21 ] = 0.0 ;
            tData[ 22 ] = 0.0 ;
            tData[ 23 ] = 0.0 ;
            tData[ 24 ] = 0.5 ;
            tData[ 25 ] = 0.0 ;
            tData[ 26 ] = 0.0 ;
            tData[ 27 ] = 0.0 ;
            tData[ 28 ] = 0.5 ;
            tData[ 29 ] = 0.0 ;
            tData[ 30 ] = 0.5 ;
            tData[ 31 ] = 0.0 ;
            tData[ 32 ] = 0.0 ;
            tData[ 33 ] = 0.0 ;
            tData[ 34 ] = 0.0 ;
            tData[ 35 ] = 0.0 ;
            tData[ 36 ] = 0.0 ;
            tData[ 37 ] = 0.0 ;
            tData[ 38 ] = 0.0 ;
            tData[ 39 ] = 0.0 ;
            tData[ 40 ] = 1.0 ;
            tData[ 41 ] = 0.0 ;
            tData[ 42 ] = 0.0 ;
            tData[ 43 ] = 0.0 ;
            tData[ 44 ] = 0.0 ;
            tData[ 45 ] = 0.0 ;
            tData[ 46 ] = 0.0 ;
            tData[ 47 ] = 0.0 ;
            tData[ 48 ] = 0.0 ;
            tData[ 49 ] = 0.0 ;
            tData[ 50 ] = 0.5 ;
            tData[ 51 ] = 0.0 ;
            tData[ 52 ] = 0.5 ;
            tData[ 53 ] = 0.0 ;
            tData[ 54 ] = 0.0 ;
            tData[ 55 ] = 0.0 ;
            tData[ 56 ] = 0.5 ;
            tData[ 57 ] = 0.0 ;
            tData[ 58 ] = 0.0 ;
            tData[ 59 ] = 0.0 ;
            tData[ 60 ] = 0.5 ;
            tData[ 61 ] = 0.0 ;
            tData[ 62 ] = 0.0 ;
            tData[ 63 ] = 0.0 ;
            tData[ 64 ] = 0.0 ;
            tData[ 65 ] = 0.0 ;
            tData[ 66 ] = 0.0 ;
            tData[ 67 ] = 0.0 ;
            tData[ 68 ] = 0.5 ;
            tData[ 69 ] = 0.0 ;
            tData[ 70 ] = 0.5 ;
            tData[ 71 ] = 0.0 ;
            tData[ 72 ] = 0.0 ;
            tData[ 73 ] = 0.0 ;
            tData[ 74 ] = 0.0 ;
            tData[ 75 ] = 0.0 ;
            tData[ 76 ] = 0.0 ;
            tData[ 77 ] = 0.0 ;
            tData[ 78 ] = 0.0 ;
            tData[ 79 ] = 0.0 ;
            tData[ 80 ] = 1.0 ;
        }
//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */
#endif //BELFEM_FN_IDENTITY_HPP
