//
// Created by Christian Messe on 2018-12-22.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp"
#include "fn_dot.hpp"

using namespace belfem;

Vector< real > tVector =  { 4, 8, 15, 16, 23, 42 };

TEST( LINALG, vector_initialization )
{
    // test length
    EXPECT_EQ( tVector.length(), (uint) 6 );

    // test content
    EXPECT_EQ( tVector( 0 ), 4 );
    EXPECT_EQ( tVector( 1 ), 8 );
    EXPECT_EQ( tVector( 2 ), 15 );
    EXPECT_EQ( tVector( 3 ), 16 );
    EXPECT_EQ( tVector( 4 ), 23 );
    EXPECT_EQ( tVector( 5 ), 42 );

    // test writing function
    tVector( 0 ) = -1;
    EXPECT_EQ( tVector( 0 ), -1 );
    tVector( 0 ) = 4;
}

TEST( LINALG, vector_minmax )
{
    EXPECT_EQ( min( tVector ), 4 );
    EXPECT_EQ( max( tVector ), 42 );
}

TEST( LINALG, vector_sum )
{
    EXPECT_EQ( sum( tVector ), 108 );
}

TEST( LINALG, vector_unique )
{
    Vector<real> tUniqueVector = {4, 4, 16, 15, 8, 42, 23, 42, 8, 8};

    // call unique funciton
    unique( tUniqueVector );

    // test content
    EXPECT_EQ( tUniqueVector( 0 ), 4 );
    EXPECT_EQ( tUniqueVector( 1 ), 8 );
    EXPECT_EQ( tUniqueVector( 2 ), 15 );
    EXPECT_EQ( tUniqueVector( 3 ), 16 );
    EXPECT_EQ( tUniqueVector( 4 ), 23 );
    EXPECT_EQ( tUniqueVector( 5 ), 42 );

    // test content
    EXPECT_EQ( tUniqueVector.length(), (uint) 6 );
}

TEST( LINALG, vector_equal_equal )
{
    Vector<real> tA = { 4, 3, 2, 1 };
    Vector<real> tB = { 4, 3, 2, 1 };
    Vector<real> tC = { 4, 4, 4, 4 };

    EXPECT_TRUE(  tA == tB );
    EXPECT_FALSE( tA == tC );
    EXPECT_TRUE(  tC == 4.0 );
    EXPECT_FALSE( tA == 4.0 );
}