//
// Created by Christian Messe on 15.09.19.
//
#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "fn_norm.hpp"
#include "fn_r2.hpp"
#include "fn_cardano.hpp"

using namespace belfem;

TEST( MATHTOOLS, cardano )
{
    // case with one real solution
    Vector< real > tA = { 1.0, 3.0, 9.0, 9.0 };
    Vector< real > tX;

    cardano( tA, tX );


    Vector< real > tExpect = { -1.327480002073326 };

    EXPECT_NEAR( r2( tX, tExpect ), 1.0, 1e-12 );

    // case with three real solutions
    tA = { 2.0, -3.0, -3.0, 2.0 };
    tExpect = { -1.0, 0.5, 2.0 };

    cardano( tA, tX );
    EXPECT_NEAR( r2( tX, tExpect ), 1.0, 1e-12 );
}