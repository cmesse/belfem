//
// Created by Christian Messe on 04.10.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_r2.hpp"
#include "fn_cubic_bezier.hpp"

using namespace belfem;

TEST( MATHTOOLS, cubic_bezier )
{
    // Control points of this curve
    Matrix< real > tPoints = { {0.0, 1./3., 1./3., 0.0},
                               {-0.5, -1./6., 1./6., 0.5}  };

    // Work vector
    Vector< real > tWork( 4 );

    // Parameter
    real tXi = -0.5;

    // Expected solution
    Vector< real > tExpect = { 0.1875, -0.25 };

    Vector< real > tResult( 2 );


    // test point function
    cubic_bezier(tPoints, tWork, tXi, tResult );

    EXPECT_NEAR( r2( tResult, tExpect ), 1.0, 1e-12 );

    // test derivative
    tExpect = { 0.25, 0.5 };
    cubic_bezier_derivative(tPoints, tWork, tXi, tResult );

    EXPECT_NEAR( r2( tResult, tExpect ), 1.0, 1e-12 );
}