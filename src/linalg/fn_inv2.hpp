//
// Created by christian on 5/9/23.
//

#ifndef BELFEM_FN_INV2_HPP
#define BELFEM_FN_INV2_HPP

#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    inline real
    inv2( const Matrix< real > & aA, Matrix< real > & aB )
    {
        // compute determinant
        real aDetA = aA( 0, 0 )*aA( 1,1 )
                   - aA( 1, 0 )*aA( 0,1 );

        aB( 0, 0 ) =   aA( 1, 1 );
        aB( 1, 0 ) = - aA( 1, 0 );
        aB( 0, 1 ) = - aA( 0, 1 );
        aB( 1, 1 ) =   aA( 0, 0 );

        aB /= aDetA ;

        return aDetA ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_INV2_HPP
