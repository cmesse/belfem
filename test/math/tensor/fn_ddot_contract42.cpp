//
// Created by Christian Messe on 18.12.20.
//

#include <gtest/gtest.h>
#include "typedefs.hpp"
#include "cl_Tensor.hpp"
#include "cl_Matrix.hpp"
#include "fn_ddot.hpp"

using namespace belfem ;

TEST( TENSOR, ddot_contract42 )
{
    Matrix< real > tA( 3, 3 );
    Tensor< real > tB( 3, 3, 3, 3 );
    Matrix< real > tC ( 3, 3 );

    Matrix< real > tD( 3, 3 );

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
            tC( i, j  ) = tValue ;
            tValue += 1. ;
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
                    tA( i, j ) += tB( i, j, k, l ) * tC( k, l );
                }
            }
        }
    }

    // 4x2 tensor contraction
    tD = tB % tC ;

    EXPECT_TRUE( tA == tD );

    tD.fill( 0.0 );
    ddot( tB, tC, tD );

    EXPECT_TRUE( tA == tD );
}