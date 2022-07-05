//
// Created by Christian Messe on 30.11.19.
//

#ifndef BELFEM_FN_CIRCLE_FROM_POINTS_HPP
#define BELFEM_FN_CIRCLE_FROM_POINTS_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"
#include "fn_norm.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    // aCircle : { Xm, Ym, R }
    inline void
    circle_from_points(
            const Vector< real > & aX,
            const Vector< real > & aY,
                  Vector< real > & aCircle )
    {
        BELFEM_ASSERT( aX.length() == 3, "aX must have length 3 " );
        BELFEM_ASSERT( aY.length() == 3, "aX must have length 3 " );
        // initial solution by intersecting two lines
        // calculated by creating lines between 0-1 and 1-2 and intersecting
        // the normals
        aCircle.set_size( 3 );


        // initial guess for Xm
        aCircle( 0 ) = aX(2)*aX(2)*(aY(0) - aY(1))
                     + (aX(0)*aX(0) + (aY(0) - aY(1))*(aY(0) - aY(2)))*(aY(1) - aY(2))
                     + aX(1)*aX(1)*(aY(2) + aY(0) );

        // initial guess for Ym
        aCircle( 1 ) =  -(aX(1)*aX(1)*aX(2)) + aX(0)*aX(0)*(-aX(1) + aX(2))
                        + aX(2)*(aY(0) - aY(1))*(aY(0) + aY(1)) +  aX(0)*(aX(1)*aX(1)
                        - aX(2)*aX(2) + aY(1)*aY(1) - aY(2)*aY(2)) + aX(1)
                        *(aX(2)*aX(2) - aY(0)*aY(0) + aY(2)*aY(2));

        real tDiv = 2.0*(aX(2)*(aY(0) - aY(1)) + aX(0)*(aY(1)
                - aY(2)) + aX(1)*(aY(2) + aY(0)));

        aCircle( 0 ) /= tDiv;
        aCircle( 1 ) /= tDiv;

        // initial guess for radius
        aCircle( 2 ) = ( std::sqrt(
                  std::pow( aX( 0 ) - aCircle( 0 ), 2 )
                + std::pow( aY( 0 ) - aCircle( 1 ), 2 ) )
                + std::sqrt(
                  std::pow( aX( 1 ) - aCircle( 0 ), 2 )
                + std::pow( aY( 1 ) - aCircle( 1 ), 2 ) )
                + std::sqrt(
                  std::pow( aX( 2 ) - aCircle( 0 ), 2 )
                + std::pow( aY( 2 ) - aCircle( 1 ), 2 ) )   ) / 3.0;


        real tResidual = BELFEM_REAL_MAX;

        // initial counter
        uint tCount = 0;

        // Jacobian
        Matrix< real > tJ( 3, 3 );

        // right hand side
        Vector< real > tF ( 3 );

        // for gesv
        Vector< int > tP( 3 );

        const real tEpsilon = 1e-9;

        while( tResidual > tEpsilon )
        {
            tJ( 0, 0 ) = aCircle( 0 ) - aX( 0 );
            tJ( 1, 0 ) = aCircle( 0 ) - aX( 1 );
            tJ( 2, 0 ) = aCircle( 0 ) - aX( 2 );
            tJ( 0, 1 ) = aCircle( 1 ) - aY( 0 );
            tJ( 1, 1 ) = aCircle( 1 ) - aY( 1 );
            tJ( 2, 1 ) = aCircle( 1 ) - aY( 2 );
            tJ( 0, 2 ) = -aCircle( 2 );
            tJ( 1, 2 ) = -aCircle( 2 );
            tJ( 2, 2 ) = -aCircle( 2 );

            tJ *= 2.0;

            tF( 0 ) =   std::pow( aX( 0 ) - aCircle( 0 ), 2 )
                             + std::pow( aY( 0 ) - aCircle( 1 ), 2 )
                             - aCircle( 2 ) *aCircle( 2 );

            tF( 1 ) =   std::pow( aX( 1 ) - aCircle( 0 ), 2 )
                             + std::pow( aY( 1 ) - aCircle( 1 ), 2 )
                             - aCircle( 2 ) *aCircle( 2 );

            tF( 2 ) =   std::pow( aX( 2 ) - aCircle( 0 ), 2 )
                             + std::pow( aY( 2 ) - aCircle( 1 ), 2 )
                             - aCircle( 2 ) *aCircle( 2 );

            tResidual = norm( tF );

            if( tResidual > tEpsilon )
            {

                gesv( tJ, tF, tP );

                aCircle -= 0.9 * tF;
            }
            else
            {
                break;
            }

            BELFEM_ERROR( tCount++ < 1000,
                    "Infinite loop while trying to find circle" );

        }
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_CIRCLE_FROM_POINTS_HPP
