//
// Created by Christian Messe on 01.12.19.
//

#ifndef BELFEM_FN_QUADRATIC_GRADIENT_HPP
#define BELFEM_FN_QUADRATIC_GRADIENT_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "assert.hpp"

namespace belfem
{
    /**
     * compute the derivative of F to X at index
     * @tparam T
     * @param aF
     * @param aX
     * @param aIndex
     * @return
     */
    template < typename  T >
    T
    quadratic_gradient(
            const Vector< T > & aF,
            const Vector< T > & aX,
            const index_t aIndex=1 )
    {
        // BELFEM_ASSERT( aX.length() == aF.length(), "Lengths do not match" );
        BELFEM_ASSERT( aIndex <  aF.length(), "Invalid Index" );

        index_t n = aX.length();

        if( aIndex == 0 )
        {
            return   ( aF( 0 ) - aF( 1 ) ) / ( aX( 0 ) - aX( 1 ) )
                   + ( aF( 0 ) - aF( 2 ) ) / ( aX( 0 ) - aX( 2 ) )
                   + ( aF( 2 ) - aF( 1 ) ) / ( aX( 1 ) - aX( 2 ) );
        }
        else if ( aIndex == aX.length() - 1 )
        {
            return  ( aF( n - 2 ) - aF( n - 3 ) ) / ( aX( n - 3 ) - aX( n - 2 ) )
                  + ( aF( n - 3 ) - aF( n - 1 ) ) / ( aX( n - 3 ) - aX( n - 1 ) )
                  + ( aF( n - 2 ) - aF( n - 1 ) ) / ( aX( n - 2 ) - aX( n - 1 ) );
        }
        else
        {
            return  ( aF( aIndex - 1 ) - aF( aIndex     ) )/( aX( aIndex-1 ) - aX( aIndex   ) )
                  + ( aF( aIndex + 1 ) - aF( aIndex - 1 ) )/( aX( aIndex-1 ) - aX( aIndex + 1 ) )
                  + ( aF( aIndex     ) - aF( aIndex + 1 ) )/( aX( aIndex   ) - aX( aIndex + 1 ) );
        }

    }

}
#endif //BELFEM_FN_QUADRATIC_GRADIENT_HPP
