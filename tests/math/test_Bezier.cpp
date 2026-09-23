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
#include <algorithm>
#include <cmath>
#include <random>

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

// =============================================================================
// the inverse map xi_by_x / xi_by_y
// =============================================================================

namespace
{
    // The bisection in xi_by_x stops at an ABSOLUTE residual of 1e-12 in x, so
    // the round-trip error in xi is bounded by 1e-12 / min | dx/dxi |. On the
    // Silver net below the smallest slope is at the right end,
    // 1.5 * 0.0385 = 0.058, and on the random nets it is 1.5 * 0.05 = 0.075;
    // 2e-11 covers both.
    const real tInverseTol = 2e-11;

    // the Q half of the control net of Silver::create_cp
    // ( src/physics/materials/cl_Material_Silver.cpp:151-152 ). The math test
    // suite does not link the materials module, so the net is reproduced here.
    // Qx is monotone and Qy is NOT, which is why only the x half of the round
    // trip is pinned on this net.
    const belfem::Vector< real > tSilverQx =
            { 4.41755060119, 5.18492206955, 7.07069952285, 7.10918477700 };
    const belfem::Vector< real > tSilverQy =
            { 5.15013111673, 5.58938782216, 5.41551370364, 5.77804911703 };

    //! strictly increasing control coordinates with a minimum span are what
    //! make the cubic monotone, hence its inverse single valued
    bool
    spans_are_wide_enough( const real * aP, const real aMinSpan )
    {
        for( belfem::uint k = 0; k < 3; ++k )
        {
            if ( aP[ k + 1 ] - aP[ k ] < aMinSpan )
            {
                return false;
            }
        }
        return true;
    }

    //! four sorted draws, redrawn until every span clears aMinSpan
    void
    draw_monotone_net( std::mt19937 & aGen,
                       std::uniform_real_distribution< real > & aDraw,
                       belfem::Vector< real > & aP,
                       const real aMinSpan )
    {
        real tValues[ 4 ];

        do
        {
            for( belfem::uint k = 0; k < 4; ++k )
            {
                tValues[ k ] = aDraw( aGen );
            }
            std::sort( tValues, tValues + 4 );
        }
        while( ! spans_are_wide_enough( tValues, aMinSpan ) );

        aP.set_size( 4 );
        for( belfem::uint k = 0; k < 4; ++k )
        {
            aP( k ) = tValues[ k ];
        }
    }
}

TEST( BezierCurve, XiByXRoundTripsOnTheSilverQNet )
{
    // Mutation: the midpoint in the bisection written as
    // 0.5 * ( tXi1 - tXi0 ) instead of 0.5 * ( tXi0 + tXi1 ) turns this red.
    const belfem::Bezier tCurve( tSilverQx, tSilverQy );

    for( belfem::uint k = 0; k <= 200; ++k )
    {
        const real tXi = -1.0 + 2.0 * ( real ) k / 200.0;

        EXPECT_NEAR( tCurve.xi_by_x( tCurve.x_by_xi( tXi ) ), tXi, tInverseTol )
                << "at xi = " << tXi;
    }
}

TEST( BezierCurve, RoundTripsOnRandomMonotoneNets )
{
    // a fixed seed keeps the fixture reproducible on one standard library.
    // The engine is specified by the standard but the distribution is not, so
    // another library may draw other nets; the round trip holds on any net
    // this loop accepts
    std::mt19937 tGen( 20260922 );
    std::uniform_real_distribution< real > tDraw( 0.0, 10.0 );

    const real tMinSpan = 0.05;

    belfem::Vector< real > tX;
    belfem::Vector< real > tY;

    for( belfem::uint n = 0; n < 20; ++n )
    {
        draw_monotone_net( tGen, tDraw, tX, tMinSpan );
        draw_monotone_net( tGen, tDraw, tY, tMinSpan );

        const belfem::Bezier tCurve( tX, tY );

        for( belfem::uint k = 0; k <= 40; ++k )
        {
            const real tXi = -1.0 + 2.0 * ( real ) k / 40.0;

            EXPECT_NEAR( tCurve.xi_by_x( tCurve.x_by_xi( tXi ) ), tXi, tInverseTol )
                    << "net " << n << ", x, at xi = " << tXi;

            EXPECT_NEAR( tCurve.xi_by_y( tCurve.y_by_xi( tXi ) ), tXi, tInverseTol )
                    << "net " << n << ", y, at xi = " << tXi;
        }
    }
}

TEST( BezierCurve, InverseSaturatesOutsideTheRange )
{
    // outside [ P0, P3 ] the cubic can fold back on itself, so the inverse
    // saturates at the ends of the curve rather than searching there. The
    // comparisons are <= and >=, so the ends themselves take the same branch
    // and return exactly.
    // Mutation: swap the two return values and this goes red.
    const belfem::Bezier tCurve( tSilverQx, tSilverQy );

    EXPECT_EQ( tCurve.xi_by_x( tSilverQx( 0 ) - 1.0 ), -1.0 );
    EXPECT_EQ( tCurve.xi_by_x( tSilverQx( 3 ) + 1.0 ),  1.0 );
    EXPECT_EQ( tCurve.xi_by_x( tSilverQx( 0 ) ),       -1.0 );
    EXPECT_EQ( tCurve.xi_by_x( tSilverQx( 3 ) ),        1.0 );

    // the y half needs a net whose y control coordinates are monotone, which
    // the Silver net is not
    const belfem::Vector< real > tMonoX = { 0.0, 1.0, 2.0, 3.0 };
    const belfem::Vector< real > tMonoY = { 1.0, 2.5, 3.0, 4.25 };

    const belfem::Bezier tMono( tMonoX, tMonoY );

    EXPECT_EQ( tMono.xi_by_y( tMonoY( 0 ) - 1.0 ), -1.0 );
    EXPECT_EQ( tMono.xi_by_y( tMonoY( 3 ) + 1.0 ),  1.0 );
    EXPECT_EQ( tMono.xi_by_y( tMonoY( 0 ) ),       -1.0 );
    EXPECT_EQ( tMono.xi_by_y( tMonoY( 3 ) ),        1.0 );
}
