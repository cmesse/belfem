/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the cubic Bezier class in src/numerics/bezier, pinning the
 * derivative chain dy/dx, d2y/dx2 and d3y/dx3 and its mirror in x. The group
 * is named BezierCurve because test_MathTools.cpp already owns the name
 * Bezier for the free function in fn_cubic_bezier.hpp.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Bezier.hpp"

namespace
{
    using belfem::real;

    // the third derivative is constant wherever x is affine in xi, so the
    // bisection in xi_by_x cannot contribute error here
    const real tExactTol = 1e-9;
}

// =============================================================================
// third derivative on an affine-x curve
// =============================================================================

TEST( BezierCurve, ThirdDerivativeIsConstantForAffineX )
{
    // equally spaced x control points give x( xi ) = 1.5 ( xi + 1 ), hence
    // d3y/dx3 = 6 ( P3 - 3 P2 + 3 P1 - P0 ) / L^3 = 6 * 8 / 27 = 16 / 9
    // everywhere. This is what pins the [-1,1] normalization of the basis:
    // a Bernstein third derivative written for [0,1] is eight times too large.
    const belfem::Vector< real > tX = {  0.0, 1.0,  2.0, 3.0 };
    const belfem::Vector< real > tY = {  0.0, 1.0, -1.0, 2.0 };

    const belfem::Bezier tCurve( tX, tY, belfem::BezierType::Horizontal );

    const real tExpect = 16.0 / 9.0;

    EXPECT_NEAR( tCurve.d3ydx3( 0.4 ), tExpect, tExactTol );
    EXPECT_NEAR( tCurve.d3ydx3( 1.3 ), tExpect, tExactTol );
    EXPECT_NEAR( tCurve.d3ydx3( 2.7 ), tExpect, tExactTol );
}

// =============================================================================
// the consumer's constructor against the Hermite closed form
// =============================================================================

TEST( BezierCurve, ThirdDerivativeMatchesHermiteForTheHorizontalConstructor )
{
    // compute_basis_xwise places the x control points at equal spacing, so
    // y( x ) is the cubic Hermite interpolant on L = x1 - x0 and its third
    // derivative is 12 ( y0 - y1 ) / L^3 + 6 ( s0 + s1 ) / L^2
    const real tX0 = 1.0;
    const real tY0 = 0.2;
    const real tS0 = 0.3;
    const real tX1 = 4.0;
    const real tY1 = 1.1;
    const real tS1 = -0.4;

    const belfem::Bezier tCurve( tX0, tY0, tS0, tX1, tY1, tS1,
                                 belfem::BezierType::Horizontal );

    const real tL      = tX1 - tX0;
    const real tExpect = 12.0 * ( tY0 - tY1 ) / ( tL * tL * tL )
                       +  6.0 * ( tS0 + tS1 ) / ( tL * tL );

    EXPECT_NEAR( tCurve.d3ydx3( 1.7 ), tExpect, tExactTol );
    EXPECT_NEAR( tCurve.d3ydx3( 3.2 ), tExpect, tExactTol );

    // the end slopes are the interpolation data; xi_by_x saturates there, so
    // both are reproduced without a bisection step
    EXPECT_NEAR( tCurve.dydx( tX0 ), tS0, tExactTol );
    EXPECT_NEAR( tCurve.dydx( tX1 ), tS1, tExactTol );
}

// =============================================================================
// the full chain rule on a non-affine curve
// =============================================================================

TEST( BezierCurve, DerivativeChainOnANonAffineCurve )
{
    // the production shape: the metals assign basis_x() directly with unequal
    // spacings, so x( xi ) is not affine and the x'', x''' terms of the chain
    // rule carry weight. dx/dxi runs from 0.75 at the ends to 1.875 in the
    // middle and never changes sign, so the curve is a graph over x.
    const belfem::Vector< real > tX = { 0.0, 0.5,  2.5, 3.0 };
    const belfem::Vector< real > tY = { 0.0, 1.0, -1.0, 2.0 };

    const belfem::Bezier tCurve( tX, tY, belfem::BezierType::Horizontal );

    // one step for all three orders. The measured error of the central
    // difference is 1.6e-8 for dy/dx, 6.7e-8 for d2y/dx2 and 4.4e-7 for
    // d3y/dx3 at these points; below h = 1e-4 the 1e-12 bisection tolerance
    // of xi_by_x starts to dominate, above it the h^2 truncation does.
    const real tStep  = 1e-4;
    const real tTolD1 = 1e-6;
    const real tTolD2 = 1e-5;
    const real tTolD3 = 1e-5;

    const real tPoints[ 3 ] = { 0.35, 1.2, 2.55 };

    for ( belfem::uint p = 0; p < 3; ++p )
    {
        const real x = tPoints[ p ];

        const real tFdD1 = ( tCurve.y( x + tStep ) - tCurve.y( x - tStep ) )
                           / ( 2.0 * tStep );
        EXPECT_NEAR( tCurve.dydx( x ), tFdD1, tTolD1 ) << "at x = " << x;

        const real tFdD2 = ( tCurve.dydx( x + tStep ) - tCurve.dydx( x - tStep ) )
                           / ( 2.0 * tStep );
        EXPECT_NEAR( tCurve.d2ydx2( x ), tFdD2, tTolD2 ) << "at x = " << x;

        const real tFdD3 = ( tCurve.d2ydx2( x + tStep ) - tCurve.d2ydx2( x - tStep ) )
                           / ( 2.0 * tStep );
        EXPECT_NEAR( tCurve.d3ydx3( x ), tFdD3, tTolD3 ) << "at x = " << x;
    }
}

// =============================================================================
// the mirror in x
// =============================================================================

TEST( BezierCurve, MirrorThirdDerivativeForVertical )
{
    // ThirdDerivativeIsConstantForAffineX with the axes swapped: y is now the
    // affine coordinate and d3x/dy3 is the constant. The BezierType argument
    // of the Vector constructor is not read -- swapping the two control-point
    // vectors is what makes this the mirror.
    const belfem::Vector< real > tX = { 0.0, 1.0, -1.0, 2.0 };
    const belfem::Vector< real > tY = { 0.0, 1.0,  2.0, 3.0 };

    const belfem::Bezier tCurve( tX, tY, belfem::BezierType::Vertical );

    const real tExpect = 16.0 / 9.0;

    EXPECT_NEAR( tCurve.d3xdy3( 0.4 ), tExpect, tExactTol );
    EXPECT_NEAR( tCurve.d3xdy3( 1.3 ), tExpect, tExactTol );
    EXPECT_NEAR( tCurve.d3xdy3( 2.7 ), tExpect, tExactTol );
}
