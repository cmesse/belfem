//
// Created by christian on 5/9/23.
//

#ifndef BELFEM_FN_INV3_HPP
#define BELFEM_FN_INV3_HPP

#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    inline real
    inv3( const Matrix< real > & aA, Matrix< real > & aB )
    {
        // compute determinant
        real aDetJ =   aA( 0, 0 ) * ( aA( 1, 1 ) * aA( 2, 2 ) - aA( 1, 2 ) * aA( 2, 1 ) )
                     + aA( 0, 1 ) * ( aA( 1, 2 ) * aA( 2, 0 ) - aA( 1, 0 ) * aA( 2, 2 ) )
                     + aA( 0, 2 ) * ( aA( 1, 0 ) * aA( 2, 1 ) - aA( 1, 1 ) * aA( 2, 0 ) );

        aB( 0, 0 ) = aA( 1, 1 ) * aA( 2, 2 ) - aA( 1, 2 ) * aA( 2, 1 ) ;
        aB( 1, 0 ) = aA( 1, 2 ) * aA( 2, 0 ) - aA( 1, 0 ) * aA( 2, 2 ) ;
        aB( 2, 0 ) = aA( 1, 0 ) * aA( 2, 1 ) - aA( 1, 1 ) * aA( 2, 0 ) ;

        aB( 0, 1 ) = aA( 0, 2 ) * aA( 2, 1 ) - aA( 0, 1 ) * aA( 2, 2 ) ;
        aB( 1, 1 ) = aA( 0, 0 ) * aA( 2, 2 ) - aA( 0, 2 ) * aA( 2, 0 ) ;
        aB( 2, 1 ) = aA( 0, 1 ) * aA( 2, 0 ) - aA( 0, 0 ) * aA( 2, 1 ) ;

        aB( 0, 2 ) = aA( 0, 1 ) * aA( 1, 2 ) - aA( 0, 2 ) * aA( 1, 1 ) ;
        aB( 1, 2 ) = aA( 0, 2 ) * aA( 1, 0 ) - aA( 0, 0 ) * aA( 1, 2 ) ;
        aB( 2, 2 ) = aA( 0, 0 ) * aA( 1, 1 ) - aA( 0, 1 ) * aA( 1, 0 ) ;

        // divide by determinant
        aB /= aDetJ ;

        // return determinant
        return aDetJ ;
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_INV3_HPP
