//
// Created by Christian Messe on 2019-01-11.
//

#ifndef BELFEM_FN_CREATE_CUBE_HPP
#define BELFEM_FN_CREATE_CUBE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
namespace belfem
{

    /**
     * create a polynomial so that
     * f = a*x^3 + b*x^2 + c*x + d
     * @tparam T
     * @param aX = { x0, x1 }
     * @param aF = { f0, dfdx0, f1, dfdx1 }
     * @param aC = { a, b, c, d }
     */
    template < typename  T >
    void
    create_truss_poly(
            const Vector< T > & aX,
            const Vector< T > & aF,
                  Vector< T > & aC )
    {
        // help variable
        real tG = aX( 0 ) - aX( 1 );
        real tH = aX( 0 ) + aX( 1 );

        aC.set_size( 4 );

        // a
        aC( 0 ) = ( 2.0*( aF( 2 ) - aF( 0 ) ) + (aF( 1 ) + aF( 3 ) )* tG )*tG;

        // b
        aC( 1 ) = tG * (3*tH*( aF( 0 ) - aF( 2 ) ) - tG*( aF( 1 )*aX( 0 )
                  + 2.0*( aF( 3 )*aX( 0 ) + aF( 1 )*aX( 1 ) ) + aF( 3 )*aX( 1 ) ) );

        // c
        aC( 2 ) = tG * ( aF( 3 )*std::pow( aX( 0 ), 3 ) + aX( 0 )*aX( 1 )*(
                  6.0 * ( aF( 2 ) - aF( 0 ) ) +
                  aX( 0 ) * ( 2.0 * aF( 1 ) + aF( 3 ) )
        - aX( 1 ) * (aF( 1 ) + 2.0 * aF( 3 ) ) ) - aF( 1 )*std::pow( aX( 1 ), 3 ) );

        // d
        aC( 3 ) = tG * ( std::pow( aX( 0 ), 2 )
                *( aF( 2 )*(aX( 0 ) - 3.0*aX( 1 )) - tG*aX( 1 )*aF( 3 ) )
                  - std::pow( aX( 1 ), 2 ) *( aF( 0 )*(aX( 1 ) - 3.0*aX( 0 ) )
                  + tG*aX( 0 )*aF( 1 ) ) );

        // scale parameters to correct size
        aC *= std::pow( tG, -4 );
    }
}
#endif //BELFEM_FN_CREATE_CUBE_HPP
