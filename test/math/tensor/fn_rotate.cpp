//
// Created by Christian Messe on 19.12.20.
//


#include <gtest/gtest.h>
#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Tensor.hpp"
#include "cl_Matrix.hpp"
#include "fn_rotate.hpp"
#include "fn_norm.hpp"
#include "fn_rotation_matrix.hpp"
#include "fn_inv.hpp"
using namespace belfem ;

TEST( TENSOR, rotate )
{
    // define an axis
    Vector< real > tAxis = { 3.0, 5.0, 7.0 };
    tAxis /= norm( tAxis );

    // define an angle
    real tAngle = 42.0 * constant::deg ;

    // create the rotation matrix
    Matrix< real > tR( 3, 3 );
    rotation_matrix( tAxis, tAngle, tR );

    // example data for C/C-SiC - like

    real E1 = 60. ;
    real E2 = 65. ;
    real E3 = 20. ;
    real G23 = 8. ;
    real G13 = 8.5 ;
    real G12 = 9. ;
    real nu23 = 0.2 ;
    real nu13 = 0.25 ;
    real nu12 = 0.1 ;


    // create compliance matrix
    Matrix< real > tS( 6, 6, 0.0 );

    tS( 0, 0 ) = 1. / E1 ;
    tS( 1, 0 ) = -nu12/E1 ;
    tS( 2, 0 ) = -nu13 / E1 ;

    tS( 0, 1 ) =  -nu12/E1 ;
    tS( 1, 1 ) =  1./ E2 ;
    tS( 2, 1 ) = -nu23 / E2 ;

    tS( 0, 2 ) = -nu13 / E1 ;
    tS( 1, 2 ) = -nu23 / E2 ;
    tS( 2, 2 ) = 1. / E3 ;

    tS( 3, 3 ) = 1.0 / G23 ;
    tS( 4, 4 ) = 1.0 / G13 ;
    tS( 5, 5 ) = 1.0 / G12 ;

    // create elasiticity matrix
    Matrix< real > tE = inv( tS );

    // create the elasticity tensor
    Tensor< real > tC( tE );

    // manual rotation
    Tensor< real > tA( 3, 3, 3, 3 );

    tA.fill( 0.0 );

    for( uint i=0; i<3; ++i )
    {
        for( uint j=0; j<3; ++j )
        {
            for( uint k=0; k<3; ++k )
            {
                for( uint l=0; l<3; ++l )
                {
                    for( uint m=0; m<3; ++m )
                    {
                        for( uint n=0; n<3; ++n )
                        {
                            for( uint o=0; o<3; ++o )
                            {
                                for( uint p=0; p<3; ++p )
                                {
                                    tA( m, n, o, p )
                                        += tC( i, j, k, l )
                                                * tR( i, m )
                                                * tR( j, n )
                                                * tR( k, o )
                                                * tR( l, p );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // functional rotation
    Tensor< real > tB( 3, 3, 3, 3 );

    rotate( tC, tR, tB );

    // create a vector to check the error
    Vector < real > tError( 81 );

    real * tX = tA.data();
    real * tY = tB.data();

    for( uint k=0; k<81; ++k )
    {
        tError( k ) = tX[ k ] - tY[ k ] ;
    }
    // hope this works
    EXPECT_NEAR( norm( tError ),  0.0,  1e-12 );

}