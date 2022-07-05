//
// Created by Christian Messe on 18.12.20.
//

#include <gtest/gtest.h>
#include "typedefs.hpp"
#include "cl_Tensor.hpp"
#include "fn_ddot.hpp"
using namespace belfem ;

TEST( TENSOR, ddot_contract44 )
{
    Tensor< real > tA( 3, 3, 3, 3 );
    Tensor< real > tB( 3, 3, 3, 3 );
    Tensor< real > tC( 3, 3, 3, 3 );

    Tensor< real > tD( 3, 3, 3, 3 );

    real tValue = 1.0 ;

    for( uint i=0; i<3; ++i )
    {
        for ( uint j = 0; j < 3; ++j )
        {
            for ( uint k = 0; k < 3; ++k )
            {
                for ( uint l = 0; l < 3; ++l )
                {
                    tB( i, j, k, l ) = tValue ;
                    tValue += 1. ;
                }
            }
        }
    }
    for( uint i=0; i<3; ++i )
    {
        for ( uint j = 0; j < 3; ++j )
        {
            for ( uint k = 0; k < 3; ++k )
            {
                for ( uint l = 0; l < 3; ++l )
                {
                    tC( i, j, k, l ) = tValue ;
                    tValue += 1. ;
                }
            }
        }
    }

    // manual contraction
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
                            tA( i, j, k, l ) += tB( i, j, m, n ) * tC( m, n, k, l );
                        }
                    }
                }
            }
        }
    }

    // 4x4 tensor contraction
    tD = tB % tC ;

    EXPECT_TRUE( tA == tD );

    tD.fill( 0.0 );
    ddot( tB, tC, tD );
    EXPECT_TRUE( tA == tD );

}