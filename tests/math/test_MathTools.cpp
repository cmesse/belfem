/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for standalone mathematical utility functions.
 * See: tests_05_mathtools.md §1–§9
 *
 * Includes regression tests for BUG-M1 (symratiospace midpoint) and
 * BUG-M2 (find_interval stale output).
 */

#include <gtest/gtest.h>
#include <cmath>
#include <limits>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_sign.hpp"
#include "fn_cardano.hpp"
#include "fn_circle_from_points.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_create_fifth_order_beam_poly.hpp"
#include "fn_cubic_bezier.hpp"
#include "fn_find_interval.hpp"
#include "fn_quadratic_gradient.hpp"
#include "fn_rotation_matrix.hpp"
#include "fn_symrationspace.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "fn_det.hpp"

namespace
{
    const belfem::real tEps = 1e-12;  // exact-in-theory results
    const belfem::real tTol = 1e-9;   // iterative solver results

    // Verify a cubic root by substitution: a*x³ + b*x² + c*x + d ≈ 0
    void verify_root( const belfem::Vector< belfem::real > & aA,
                      belfem::real aX )
    {
        belfem::real tVal = aA( 0 ) * aX * aX * aX
                          + aA( 1 ) * aX * aX
                          + aA( 2 ) * aX
                          + aA( 3 );
        EXPECT_NEAR( tVal, 0.0, tEps );
    }
}

// =============================================================================
// §1  Sign Function  [semantic]
// =============================================================================

TEST( Sign, SignPositive )
{
    EXPECT_EQ( belfem::sign( 5.0 ), 1 );
}

TEST( Sign, SignNegative )
{
    EXPECT_EQ( belfem::sign( -3.2 ), -1 );
}

TEST( Sign, SignZero )
{
    EXPECT_EQ( belfem::sign( 0.0 ), 0 );
}

TEST( Sign, SignNegativeZero )
{
    EXPECT_EQ( belfem::sign( -0.0 ), 0 );
}

TEST( Sign, SignInfinity )
{
    EXPECT_EQ( belfem::sign( std::numeric_limits< double >::infinity() ), 1 );
    EXPECT_EQ( belfem::sign( -std::numeric_limits< double >::infinity() ), -1 );
}

TEST( Sign, SignInteger )
{
    EXPECT_EQ( belfem::sign( 7 ), 1 );
    EXPECT_EQ( belfem::sign( -7 ), -1 );
    EXPECT_EQ( belfem::sign( 0 ), 0 );
}

// =============================================================================
// §2  Cardano  [semantic]
// =============================================================================

// --- §2.2 Degenerate cases (a=0) ---

TEST( Cardano, CardanoLinearFallback )
{
    // 0*x³ + 0*x² + 5*x + 7 = 0 → x = -7/5
    belfem::Vector< belfem::real > tA = { 0.0, 0.0, 5.0, 7.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 1u );
    verify_root( tA, tX( 0 ) );
}

TEST( Cardano, CardanoConstantNoRoots )
{
    belfem::Vector< belfem::real > tA = { 0.0, 0.0, 0.0, 1.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    EXPECT_EQ( tX.length(), 0u );
}

TEST( Cardano, CardanoQuadraticTwoRoots )
{
    // x² + 3x + 2 = (x+1)(x+2) → roots -2, -1
    belfem::Vector< belfem::real > tA = { 0.0, 1.0, 3.0, 2.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 2u );
    verify_root( tA, tX( 0 ) );
    verify_root( tA, tX( 1 ) );
    // sorted
    EXPECT_LE( tX( 0 ), tX( 1 ) );
}

TEST( Cardano, CardanoQuadraticNoRealRoots )
{
    // x² + 1 = 0 → no real roots
    belfem::Vector< belfem::real > tA = { 0.0, 1.0, 0.0, 1.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    EXPECT_EQ( tX.length(), 0u );
}

TEST( Cardano, CardanoQuadraticDoubleRoot )
{
    // (x-3)² = x² - 6x + 9 → one distinct root at 3
    // Contract: distinct roots only → length == 1
    belfem::Vector< belfem::real > tA = { 0.0, 1.0, -6.0, 9.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 1u );
    EXPECT_NEAR( tX( 0 ), 3.0, tEps );
}

TEST( Cardano, CardanoZeroAllCoeffs )
{
    belfem::Vector< belfem::real > tA = { 0.0, 0.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    EXPECT_EQ( tX.length(), 0u );
}

// --- §2.3 One real root (D > 0) ---

TEST( Cardano, CardanoOneRealRoot )
{
    // x³ + x + 1 = 0 → 1 real root
    belfem::Vector< belfem::real > tA = { 1.0, 0.0, 1.0, 1.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 1u );
    verify_root( tA, tX( 0 ) );
}

TEST( Cardano, CardanoOneRealRootKnown )
{
    // (x-2)³ = x³ - 6x² + 12x - 8
    belfem::Vector< belfem::real > tA = { 1.0, -6.0, 12.0, -8.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_GE( tX.length(), 1u );
    verify_root( tA, tX( 0 ) );
}

// --- §2.4 Three real roots (D < 0) ---

TEST( Cardano, CardanoThreeRealRoots )
{
    // (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6
    belfem::Vector< belfem::real > tA = { 1.0, -6.0, 11.0, -6.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 3u );
    for( size_t i = 0; i < 3; ++i )
    {
        verify_root( tA, tX( i ) );
    }
    // sorted
    EXPECT_LE( tX( 0 ), tX( 1 ) );
    EXPECT_LE( tX( 1 ), tX( 2 ) );
}

TEST( Cardano, CardanoThreeRealRootsIrreducible )
{
    // x³ - 3x - 1 = 0 → 3 real roots (casus irreducibilis)
    belfem::Vector< belfem::real > tA = { 1.0, 0.0, -3.0, -1.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 3u );
    for( size_t i = 0; i < 3; ++i )
    {
        verify_root( tA, tX( i ) );
    }
}

TEST( Cardano, CardanoThreeRealRootsNegativeLeading )
{
    // -1*(x-1)(x-2)(x-3) = -x³ + 6x² - 11x + 6
    belfem::Vector< belfem::real > tA = { -1.0, 6.0, -11.0, 6.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 3u );
    for( size_t i = 0; i < 3; ++i )
    {
        verify_root( tA, tX( i ) );
    }
}

// --- §2.5 Repeated root (D = 0) ---

TEST( Cardano, CardanoTripleRoot )
{
    // (x-1)³ = x³ - 3x² + 3x - 1
    // Contract: distinct roots only → 1 root
    belfem::Vector< belfem::real > tA = { 1.0, -3.0, 3.0, -1.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 1u );
    EXPECT_NEAR( tX( 0 ), 1.0, tEps );
}

TEST( Cardano, CardanoDoubleRoot )
{
    // (x-1)²(x+2) = x³ + 0x² - 3x + 2
    // Contract: return distinct real roots only, sorted → {-2, 1}
    belfem::Vector< belfem::real > tA = { 1.0, 0.0, -3.0, 2.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    ASSERT_EQ( tX.length(), 2u );
    EXPECT_NEAR( tX( 0 ), -2.0, tEps );
    EXPECT_NEAR( tX( 1 ),  1.0, tEps );

    // verify by substitution
    for( size_t i = 0; i < tX.length(); ++i )
    {
        verify_root( tA, tX( i ) );
    }
}

// --- §2.6 Robustness ---

TEST( Cardano, CardanoReturnedRootsSorted )
{
    belfem::Vector< belfem::real > tA = { 1.0, -6.0, 11.0, -6.0 };
    belfem::Vector< belfem::real > tX;
    belfem::cardano( tA, tX );

    for( size_t i = 1; i < tX.length(); ++i )
    {
        EXPECT_LE( tX( i - 1 ), tX( i ) );
    }
}

// =============================================================================
// §3  Circle from Points  [semantic]
// =============================================================================

TEST( Circle, CircleFromPointsUnitCircle )
{
    belfem::Vector< belfem::real > tX = { 1.0, 0.0, -1.0 };
    belfem::Vector< belfem::real > tY = { 0.0, 1.0, 0.0 };
    belfem::Vector< belfem::real > tCircle;
    belfem::circle_from_points( tX, tY, tCircle );

    EXPECT_NEAR( tCircle( 0 ), 0.0, tTol );   // Xm
    EXPECT_NEAR( tCircle( 1 ), 0.0, tTol );   // Ym
    EXPECT_NEAR( tCircle( 2 ), 1.0, tTol );   // R
}

TEST( Circle, CircleFromPointsTranslatedCircle )
{
    // circle with center (5, 3), radius 2
    belfem::real tCx = 5.0, tCy = 3.0, tR = 2.0;
    belfem::Vector< belfem::real > tX = {
        tCx + tR,
        tCx + tR * std::cos( 2.0 * M_PI / 3.0 ),
        tCx + tR * std::cos( 4.0 * M_PI / 3.0 ) };
    belfem::Vector< belfem::real > tY = {
        tCy,
        tCy + tR * std::sin( 2.0 * M_PI / 3.0 ),
        tCy + tR * std::sin( 4.0 * M_PI / 3.0 ) };

    belfem::Vector< belfem::real > tCircle;
    belfem::circle_from_points( tX, tY, tCircle );

    EXPECT_NEAR( tCircle( 0 ), tCx, tTol );
    EXPECT_NEAR( tCircle( 1 ), tCy, tTol );
    EXPECT_NEAR( tCircle( 2 ), tR,  tTol );
}

TEST( Circle, CircleFromPointsAllPointsEquidistant )
{
    belfem::Vector< belfem::real > tX = { 1.0, 0.0, -1.0 };
    belfem::Vector< belfem::real > tY = { 0.0, 1.0, 0.0 };
    belfem::Vector< belfem::real > tCircle;
    belfem::circle_from_points( tX, tY, tCircle );

    for( size_t i = 0; i < 3; ++i )
    {
        belfem::real tDist = std::sqrt(
            std::pow( tX( i ) - tCircle( 0 ), 2 ) +
            std::pow( tY( i ) - tCircle( 1 ), 2 ) );
        EXPECT_NEAR( tDist, tCircle( 2 ), tTol );
    }
}

// --- §3.2 debug ---

#ifndef NDEBUG

TEST( CircleDebug, CircleFromPointsWrongXLengthThrows )
{
    belfem::Vector< belfem::real > tX( 2, 1.0 );
    belfem::Vector< belfem::real > tY( 3, 1.0 );
    belfem::Vector< belfem::real > tCircle;
    EXPECT_THROW( belfem::circle_from_points( tX, tY, tCircle ),
                  std::runtime_error );
}

TEST( CircleDebug, CircleFromPointsWrongYLengthThrows )
{
    belfem::Vector< belfem::real > tX( 3, 1.0 );
    belfem::Vector< belfem::real > tY( 4, 1.0 );
    belfem::Vector< belfem::real > tCircle;
    EXPECT_THROW( belfem::circle_from_points( tX, tY, tCircle ),
                  std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §4.1  Cubic Beam Poly  [semantic]
// =============================================================================

TEST( BeamPoly, BeamPolyInterpolatesValues )
{
    // f(x) = x³ at x1=0 and x2=1
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0;     // f(0)=0, f'(0)=0
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 3.0;     // f(1)=1, f'(1)=3

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_beam_poly( tX1, tF1, tDF1, tX2, tF2, tDF2, tCoeffs );

    // evaluate at endpoints via polyval
    EXPECT_NEAR( belfem::polyval( tCoeffs, tX1 ), tF1, tTol );
    EXPECT_NEAR( belfem::polyval( tCoeffs, tX2 ), tF2, tTol );
}

TEST( BeamPoly, BeamPolyInterpolatesDerivatives )
{
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0;
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 3.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_beam_poly( tX1, tF1, tDF1, tX2, tF2, tDF2, tCoeffs );

    EXPECT_NEAR( belfem::dpolyval( tCoeffs, tX1 ), tDF1, tTol );
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, tX2 ), tDF2, tTol );
}

TEST( BeamPoly, BeamPolyKnownCubic )
{
    // f(x) = x³ → coeffs should be {1, 0, 0, 0}
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0;
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 3.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_beam_poly( tX1, tF1, tDF1, tX2, tF2, tDF2, tCoeffs );

    EXPECT_NEAR( tCoeffs( 0 ), 1.0, tTol );
    EXPECT_NEAR( tCoeffs( 1 ), 0.0, tTol );
    EXPECT_NEAR( tCoeffs( 2 ), 0.0, tTol );
    EXPECT_NEAR( tCoeffs( 3 ), 0.0, tTol );
}

TEST( BeamPoly, BeamPolyLinearFunction )
{
    // f(x) = 2x + 3, f'(x) = 2
    belfem::real tX1 = 1.0, tF1 = 5.0, tDF1 = 2.0;
    belfem::real tX2 = 3.0, tF2 = 9.0, tDF2 = 2.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_beam_poly( tX1, tF1, tDF1, tX2, tF2, tDF2, tCoeffs );

    // cubic and quadratic coeffs should be ~0
    EXPECT_NEAR( tCoeffs( 0 ), 0.0, tTol );
    EXPECT_NEAR( tCoeffs( 1 ), 0.0, tTol );
    // linear and constant
    EXPECT_NEAR( tCoeffs( 2 ), 2.0, tTol );
    EXPECT_NEAR( tCoeffs( 3 ), 3.0, tTol );
}

// =============================================================================
// §4.2  Fifth-Order Beam Poly  [semantic]
// =============================================================================

TEST( BeamPoly, FifthOrderInterpolatesValues )
{
    // f(x) = x⁵ at x1=0, x2=1
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0, tDDF1 = 0.0;
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 5.0, tDDF2 = 20.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_fifth_order_beam_poly(
        tX1, tF1, tDF1, tDDF1,
        tX2, tF2, tDF2, tDDF2,
        tCoeffs );

    EXPECT_NEAR( belfem::polyval( tCoeffs, tX1 ), tF1, tTol );
    EXPECT_NEAR( belfem::polyval( tCoeffs, tX2 ), tF2, tTol );
}

TEST( BeamPoly, FifthOrderInterpolatesFirstDerivatives )
{
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0, tDDF1 = 0.0;
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 5.0, tDDF2 = 20.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_fifth_order_beam_poly(
        tX1, tF1, tDF1, tDDF1,
        tX2, tF2, tDF2, tDDF2,
        tCoeffs );

    EXPECT_NEAR( belfem::dpolyval( tCoeffs, tX1 ), tDF1, tTol );
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, tX2 ), tDF2, tTol );
}

TEST( BeamPoly, FifthOrderInterpolatesSecondDerivatives )
{
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0, tDDF1 = 0.0;
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 5.0, tDDF2 = 20.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_fifth_order_beam_poly(
        tX1, tF1, tDF1, tDDF1,
        tX2, tF2, tDF2, tDDF2,
        tCoeffs );

    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, tX1 ), tDDF1, tTol );
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, tX2 ), tDDF2, tTol );
}

TEST( BeamPoly, FifthOrderNonZeroEndpoints )
{
    // f(x) = x⁵ at x1=1, x2=2 — exercises more Vandermonde entries
    // (idea from Gemini — non-zero x1 is a stronger test)
    belfem::real tX1 = 1.0, tF1 = 1.0, tDF1 = 5.0, tDDF1 = 20.0;
    belfem::real tX2 = 2.0, tF2 = 32.0, tDF2 = 80.0, tDDF2 = 160.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_fifth_order_beam_poly(
        tX1, tF1, tDF1, tDDF1,
        tX2, tF2, tDF2, tDDF2,
        tCoeffs );

    EXPECT_NEAR( belfem::polyval( tCoeffs, tX1 ), tF1, tTol );
    EXPECT_NEAR( belfem::polyval( tCoeffs, tX2 ), tF2, tTol );
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, tX1 ), tDF1, tTol );
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, tX2 ), tDF2, tTol );
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, tX1 ), tDDF1, tTol );
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, tX2 ), tDDF2, tTol );
}

TEST( BeamPoly, FifthOrderKnownQuintic )
{
    // f(x) = x⁵ → coeffs {1, 0, 0, 0, 0, 0}
    belfem::real tX1 = 0.0, tF1 = 0.0, tDF1 = 0.0, tDDF1 = 0.0;
    belfem::real tX2 = 1.0, tF2 = 1.0, tDF2 = 5.0, tDDF2 = 20.0;

    belfem::Vector< belfem::real > tCoeffs;
    belfem::create_fifth_order_beam_poly(
        tX1, tF1, tDF1, tDDF1,
        tX2, tF2, tDF2, tDDF2,
        tCoeffs );

    EXPECT_NEAR( tCoeffs( 0 ), 1.0, tTol );
    for( size_t i = 1; i < 6; ++i )
    {
        EXPECT_NEAR( tCoeffs( i ), 0.0, tTol );
    }
}

// =============================================================================
// §5  Cubic Bézier  [semantic]
// =============================================================================

TEST( Bezier, BezierEndpointMinus1 )
{
    // at ξ=-1, basis = {1,0,0,0} → selects column 0
    belfem::Matrix< belfem::real > tPts = { { 1.0, 2.0, 3.0, 4.0 },
                                              { 10.0, 20.0, 30.0, 40.0 } };
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 2 );

    belfem::cubic_bezier( tPts, tWork, -1.0, tPoint );

    EXPECT_NEAR( tPoint( 0 ), 1.0,  tEps );
    EXPECT_NEAR( tPoint( 1 ), 10.0, tEps );
}

TEST( Bezier, BezierEndpointPlus1 )
{
    // at ξ=+1, basis = {0,0,0,1} → selects column 3
    belfem::Matrix< belfem::real > tPts = { { 1.0, 2.0, 3.0, 4.0 },
                                              { 10.0, 20.0, 30.0, 40.0 } };
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 2 );

    belfem::cubic_bezier( tPts, tWork, 1.0, tPoint );

    EXPECT_NEAR( tPoint( 0 ), 4.0,  tEps );
    EXPECT_NEAR( tPoint( 1 ), 40.0, tEps );
}

TEST( Bezier, BezierStraightLine )
{
    // collinear: P0=(0,0), P1=(1,1), P2=(2,2), P3=(3,3)
    belfem::Matrix< belfem::real > tPts = { { 0.0, 1.0, 2.0, 3.0 },
                                              { 0.0, 1.0, 2.0, 3.0 } };
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 2 );

    // at ξ=0 midpoint
    belfem::cubic_bezier( tPts, tWork, 0.0, tPoint );
    EXPECT_NEAR( tPoint( 0 ), 1.5, tEps );
    EXPECT_NEAR( tPoint( 1 ), 1.5, tEps );
}

TEST( Bezier, BezierMidpointWeightedCombination )
{
    // at ξ=0, basis = {1, 3, 3, 1}/8 → known weighted combination
    // (idea from ChatGPT)
    belfem::Matrix< belfem::real > tPts = { { 1.0, 2.0, 3.0, 5.0 } };
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 1 );

    belfem::cubic_bezier( tPts, tWork, 0.0, tPoint );

    belfem::real tExpected = ( 1.0 + 3.0 * 2.0 + 3.0 * 3.0 + 5.0 ) / 8.0;
    EXPECT_NEAR( tPoint( 0 ), tExpected, tEps );
}

TEST( Bezier, BezierDustCleanup )
{
    // near-zero values should be cleaned to exactly 0.0 (idea from ChatGPT)
    belfem::Matrix< belfem::real > tPts = { { 1e-20, 1e-20, 1e-20, 1e-20 } };
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 1 );

    belfem::cubic_bezier( tPts, tWork, 0.0, tPoint );

    EXPECT_EQ( tPoint( 0 ), 0.0 );
}

TEST( Bezier, BezierDerivativeConsistency )
{
    belfem::Matrix< belfem::real > tPts = { { 0.0, 1.0, 3.0, 4.0 },
                                              { 0.0, 2.0, 1.0, 3.0 } };
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 2 );
    belfem::Vector< belfem::real > tDeriv( 2 );

    belfem::real tXi = 0.3;
    belfem::real tH  = 1e-6;

    // analytical derivative
    belfem::cubic_bezier_derivative( tPts, tWork, tXi, tDeriv );

    // finite difference
    belfem::Vector< belfem::real > tPplus( 2 ), tPminus( 2 );
    belfem::cubic_bezier( tPts, tWork, tXi + tH, tPplus );
    belfem::cubic_bezier( tPts, tWork, tXi - tH, tPminus );

    for( size_t i = 0; i < 2; ++i )
    {
        belfem::real tFD = ( tPplus( i ) - tPminus( i ) ) / ( 2.0 * tH );
        EXPECT_NEAR( tDeriv( i ), tFD, 1e-5 );
    }
}

// --- §5.2 debug ---

#ifndef NDEBUG

TEST( BezierDebug, BezierWrongPointLengthThrows )
{
    belfem::Matrix< belfem::real > tPts( 2, 4, 0.0 );
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 3 );   // wrong: should be 2
    EXPECT_THROW( belfem::cubic_bezier( tPts, tWork, 0.0, tPoint ),
                  std::runtime_error );
}

TEST( BezierDebug, BezierWrongColCountThrows )
{
    belfem::Matrix< belfem::real > tPts( 2, 3, 0.0 );   // wrong: need 4 cols
    belfem::Vector< belfem::real > tWork( 4 );
    belfem::Vector< belfem::real > tPoint( 2 );
    EXPECT_THROW( belfem::cubic_bezier( tPts, tWork, 0.0, tPoint ),
                  std::runtime_error );
}

TEST( BezierDebug, BezierWrongWorkLengthThrows )
{
    belfem::Matrix< belfem::real > tPts( 2, 4, 0.0 );
    belfem::Vector< belfem::real > tWork( 3 );   // wrong: need 4
    belfem::Vector< belfem::real > tPoint( 2 );
    EXPECT_THROW( belfem::cubic_bezier( tPts, tWork, 0.0, tPoint ),
                  std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §6  Find Interval  [semantic]
// =============================================================================

TEST( FindInterval, FindIntervalBelowRange )
{
    belfem::Vector< belfem::real > tData = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::index_t tIndex = 999;
    belfem::real tXi = -999.0;

    belfem::find_interval( tData, 0.5, tIndex, tXi );
    EXPECT_EQ( tIndex, 0u );
    EXPECT_NEAR( tXi, 0.5, tEps );
}

TEST( FindInterval, FindIntervalAboveRange )
{
    belfem::Vector< belfem::real > tData = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::index_t tIndex = 999;
    belfem::real tXi = -999.0;

    belfem::find_interval( tData, 3.5, tIndex, tXi );
    EXPECT_EQ( tIndex, 3u );
    EXPECT_NEAR( tXi, 0.5, tEps );
}

TEST( FindInterval, FindIntervalExactKnot )
{
    belfem::Vector< belfem::real > tData = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::index_t tIndex = 999;
    belfem::real tXi = -999.0;

    // value at first data point — handled by below-range branch
    belfem::find_interval( tData, 0.0, tIndex, tXi );
    EXPECT_EQ( tIndex, 0u );
    EXPECT_NEAR( tXi, 0.0, tEps );
}

TEST( FindInterval, FindIntervalRegressionBugM2 )
{
    // BUG-M2: after the while loop breaks at k-i==1, aIndex and aXi
    // may retain stale values from a previous iteration.
    // Correct: for value=2.5 in {0,1,2,3,4}, the interval is [2,3],
    // so aIndex should be 2 and aXi should be 0.5.
    belfem::Vector< belfem::real > tData = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::index_t tIndex = 999;
    belfem::real tXi = -999.0;

    belfem::find_interval( tData, 2.5, tIndex, tXi );

    EXPECT_EQ( tIndex, 2u );
    EXPECT_NEAR( tXi, 0.5, tEps );
}

// =============================================================================
// §7  Quadratic Gradient  [semantic]
// =============================================================================

TEST( QuadGrad, QuadraticGradientExactQuadratic )
{
    // f(x) = 3x² + 2x + 1, f'(x) = 6x + 2
    // sampled at x = {0, 1, 2, 3, 4}
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tF( 5 );
    for( size_t i = 0; i < 5; ++i )
    {
        belfem::real tXi = tX( i );
        tF( i ) = 3.0 * tXi * tXi + 2.0 * tXi + 1.0;
    }

    // interior points (index 1, 2, 3) should recover exact derivative
    for( belfem::index_t k = 1; k <= 3; ++k )
    {
        belfem::real tExpected = 6.0 * tX( k ) + 2.0;
        belfem::real tGrad = belfem::quadratic_gradient( tF, tX, k );
        EXPECT_NEAR( tGrad, tExpected, tTol );
    }
}

TEST( QuadGrad, QuadraticGradientLeftBoundary )
{
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tF( 5 );
    for( size_t i = 0; i < 5; ++i )
    {
        belfem::real tXi = tX( i );
        tF( i ) = 3.0 * tXi * tXi + 2.0 * tXi + 1.0;
    }

    belfem::real tExpected = 2.0;   // f'(0) = 6*0 + 2 = 2
    EXPECT_NEAR( belfem::quadratic_gradient( tF, tX, 0 ), tExpected, tTol );
}

TEST( QuadGrad, QuadraticGradientRightBoundary )
{
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tF( 5 );
    for( size_t i = 0; i < 5; ++i )
    {
        belfem::real tXi = tX( i );
        tF( i ) = 3.0 * tXi * tXi + 2.0 * tXi + 1.0;
    }

    belfem::real tExpected = 26.0;   // f'(4) = 6*4 + 2 = 26
    EXPECT_NEAR( belfem::quadratic_gradient( tF, tX, 4 ), tExpected, tTol );
}

TEST( QuadGrad, QuadraticGradientLinearFunction )
{
    // f(x) = 5x + 3 → f'(x) = 5 everywhere
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tF( 5 );
    for( size_t i = 0; i < 5; ++i )
    {
        tF( i ) = 5.0 * tX( i ) + 3.0;
    }

    for( belfem::index_t k = 0; k < 5; ++k )
    {
        EXPECT_NEAR( belfem::quadratic_gradient( tF, tX, k ), 5.0, tTol );
    }
}

TEST( QuadGrad, QuadraticGradientMinimumThreePoints )
{
    // Precondition: n >= 3. Test with exactly 3 points — all branches exercised.
    // f(x) = 2x² + x + 1 → f'(x) = 4x + 1
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0 };
    belfem::Vector< belfem::real > tF = { 1.0, 4.0, 11.0 };

    // index 0 → left boundary formula: f'(0) = 4*0 + 1 = 1
    belfem::real tG0 = belfem::quadratic_gradient( tF, tX, 0 );
    EXPECT_NEAR( tG0, 1.0, tTol );

    // index 1 → interior formula: f'(1) = 4*1 + 1 = 5
    belfem::real tG1 = belfem::quadratic_gradient( tF, tX, 1 );
    EXPECT_NEAR( tG1, 5.0, tTol );

    // index 2 → right boundary formula: f'(2) = 4*2 + 1 = 9
    belfem::real tG2 = belfem::quadratic_gradient( tF, tX, 2 );
    EXPECT_NEAR( tG2, 9.0, tTol );
}

// =============================================================================
// §7.2  Quadratic Gradient  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( QuadGradDebug, QuadraticGradientOutOfBoundsThrows )
{
    // (idea from Gemini)
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0 };
    belfem::Vector< belfem::real > tF = { 0.0, 1.0, 4.0 };
    EXPECT_THROW( belfem::quadratic_gradient( tF, tX, 3 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §8  Rotation Matrices  [semantic]
// =============================================================================

// --- §8.1 Axis-Angle ---

TEST( Rotation, RotationMatrixIdentity )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( tAxis, 0.0, tR );

    for( size_t i = 0; i < 3; ++i )
    for( size_t j = 0; j < 3; ++j )
    {
        belfem::real tExpected = ( i == j ) ? 1.0 : 0.0;
        EXPECT_NEAR( tR( i, j ), tExpected, tEps );
    }
}

TEST( Rotation, RotationMatrix90AboutZ )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( tAxis, M_PI / 2.0, tR );

    // R * {1,0,0} ≈ {0,1,0}
    EXPECT_NEAR( tR( 0, 0 ), 0.0, tTol );
    EXPECT_NEAR( tR( 1, 0 ), 1.0, tTol );
    EXPECT_NEAR( tR( 2, 0 ), 0.0, tTol );
}

TEST( Rotation, RotationMatrixOrthogonal )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    // normalize
    belfem::real tNorm = std::sqrt( 14.0 );
    tAxis( 0 ) /= tNorm;
    tAxis( 1 ) /= tNorm;
    tAxis( 2 ) /= tNorm;

    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( tAxis, 0.73, tR );

    // R * Rᵀ ≈ I — via operator()
    for( size_t i = 0; i < 3; ++i )
    for( size_t j = 0; j < 3; ++j )
    {
        belfem::real tDot = 0.0;
        for( size_t k = 0; k < 3; ++k )
            tDot += tR( i, k ) * tR( j, k );
        belfem::real tExpected = ( i == j ) ? 1.0 : 0.0;
        EXPECT_NEAR( tDot, tExpected, tTol );
    }
}

TEST( Rotation, RotationMatrixDetPlusOne )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::real tNorm = std::sqrt( 14.0 );
    tAxis( 0 ) /= tNorm;
    tAxis( 1 ) /= tNorm;
    tAxis( 2 ) /= tNorm;

    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( tAxis, 1.23, tR );

    EXPECT_NEAR( belfem::det( tR ), 1.0, tTol );
}

// --- §8.2 Euler Angles ---

TEST( Rotation, EulerIdentity )
{
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( 0.0, 0.0, 0.0, tR );

    for( size_t i = 0; i < 3; ++i )
    for( size_t j = 0; j < 3; ++j )
    {
        belfem::real tExpected = ( i == j ) ? 1.0 : 0.0;
        EXPECT_NEAR( tR( i, j ), tExpected, tEps );
    }
}

TEST( Rotation, EulerPureRoll )
{
    // roll = π/2, rest 0 → known matrix (idea from ChatGPT)
    // DIN 9300: roll is first rotation (about x in body frame)
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( 0.0, 0.0, M_PI / 2.0, tR );

    // verify via operator() to catch BUG-M4
    EXPECT_NEAR( tR( 0, 0 ),  0.0, tEps );
    EXPECT_NEAR( tR( 0, 1 ), -1.0, tEps );
    EXPECT_NEAR( tR( 1, 0 ),  1.0, tEps );
    EXPECT_NEAR( tR( 1, 1 ),  0.0, tEps );
    EXPECT_NEAR( tR( 2, 2 ),  1.0, tEps );
}

TEST( Rotation, EulerPureYaw )
{
    // yaw = π/2, rest 0 (idea from ChatGPT)
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( M_PI / 2.0, 0.0, 0.0, tR );

    EXPECT_NEAR( tR( 0, 0 ),  1.0, tEps );
    EXPECT_NEAR( tR( 1, 1 ),  0.0, tEps );
    EXPECT_NEAR( tR( 1, 2 ), -1.0, tEps );
    EXPECT_NEAR( tR( 2, 1 ),  1.0, tEps );
    EXPECT_NEAR( tR( 2, 2 ),  0.0, tEps );
}

TEST( Rotation, EulerPurePitch )
{
    // pitch = π/2, rest 0 (idea from ChatGPT)
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( 0.0, M_PI / 2.0, 0.0, tR );

    EXPECT_NEAR( tR( 0, 0 ),  0.0, tEps );
    EXPECT_NEAR( tR( 0, 2 ),  1.0, tEps );
    EXPECT_NEAR( tR( 1, 1 ),  1.0, tEps );
    EXPECT_NEAR( tR( 2, 0 ), -1.0, tEps );
    EXPECT_NEAR( tR( 2, 2 ),  0.0, tEps );
}

TEST( Rotation, EulerOrthogonal )
{
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( 0.3, 0.5, 0.7, tR );

    // R * Rᵀ ≈ I — via operator() to catch BUG-M4 Blaze padding
    for( size_t i = 0; i < 3; ++i )
    for( size_t j = 0; j < 3; ++j )
    {
        belfem::real tDot = 0.0;
        for( size_t k = 0; k < 3; ++k )
            tDot += tR( i, k ) * tR( j, k );
        belfem::real tExpected = ( i == j ) ? 1.0 : 0.0;
        EXPECT_NEAR( tDot, tExpected, tTol );
    }
}

TEST( Rotation, EulerDetPlusOne )
{
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( 0.3, 0.5, 0.7, tR );

    EXPECT_NEAR( belfem::det( tR ), 1.0, tTol );
}

TEST( Rotation, EulerRegressionBugM4 )
{
    // BUG-M4: Blaze path writes data[4,5,6,8,9,10] skipping 3,7.
    // Verify via operator() not data() to abstract away padding.
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::rotation_matrix( 0.4, 0.6, 0.8, tR );

    // compute expected values manually
    belfem::real tSa = std::sin( 0.8 ), tCa = std::cos( 0.8 );
    belfem::real tSb = std::sin( 0.6 ), tCb = std::cos( 0.6 );
    belfem::real tSc = std::sin( 0.4 ), tCc = std::cos( 0.4 );

    EXPECT_NEAR( tR( 0, 0 ), tCa * tCb, tTol );
    EXPECT_NEAR( tR( 1, 0 ), tCb * tSa, tTol );
    EXPECT_NEAR( tR( 2, 0 ), -tSb,      tTol );
    EXPECT_NEAR( tR( 0, 1 ), tCa * tSb * tSc - tCc * tSa, tTol );
    EXPECT_NEAR( tR( 2, 2 ), tCb * tCc, tTol );
}

// --- §8.3 Strip ---

TEST( Rotation, StripMatchesFullMatrixFirstTwoRows )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::real tNorm = std::sqrt( 14.0 );
    tAxis( 0 ) /= tNorm;
    tAxis( 1 ) /= tNorm;
    tAxis( 2 ) /= tNorm;

    belfem::Matrix< belfem::real > tFull( 3, 3 );
    belfem::rotation_matrix( tAxis, 0.5, tFull );

    belfem::Matrix< belfem::real > tStrip( 2, 3 );
    belfem::rotation_matrix_strip( tAxis, 0.5, tStrip );

    for( size_t i = 0; i < 2; ++i )
    for( size_t j = 0; j < 3; ++j )
    {
        EXPECT_NEAR( tStrip( i, j ), tFull( i, j ), tEps );
    }
}

// --- §8.4 debug ---

#ifndef NDEBUG

TEST( RotationDebug, AxisAngleWrongAxisLengthThrows )
{
    belfem::Vector< belfem::real > tAxis( 2, 1.0 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    EXPECT_THROW( belfem::rotation_matrix( tAxis, 1.0, tR ), std::runtime_error );
}

TEST( RotationDebug, AxisAngleWrongMatrixSizeThrows )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Matrix< belfem::real > tR( 2, 3 );
    EXPECT_THROW( belfem::rotation_matrix( tAxis, 1.0, tR ), std::runtime_error );
}

TEST( RotationDebug, EulerWrongMatrixSizeThrows )
{
    belfem::Matrix< belfem::real > tR( 2, 3 );
    EXPECT_THROW( belfem::rotation_matrix( 0.1, 0.2, 0.3, tR ), std::runtime_error );
}

TEST( RotationDebug, StripWrongMatrixSizeThrows )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Matrix< belfem::real > tR( 3, 3 );
    EXPECT_THROW( belfem::rotation_matrix_strip( tAxis, 1.0, tR ),
                  std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §9  Symratiospace  [semantic]
// =============================================================================

TEST( Symratio, SymratioEndpoints )
{
    belfem::Vector< belfem::real > tX;
    belfem::symratiospace( 0.0, 10.0, 1.2, 11, tX );

    EXPECT_NEAR( tX( 0 ), 0.0,  tEps );
    EXPECT_NEAR( tX( 10 ), 10.0, tEps );
}

TEST( Symratio, SymratioUniformWhenRatioOne )
{
    belfem::Vector< belfem::real > tX;
    belfem::symratiospace( 0.0, 10.0, 1.0, 11, tX );

    belfem::real tExpectedStep = 1.0;
    for( size_t i = 1; i < 11; ++i )
    {
        EXPECT_NEAR( tX( i ) - tX( i - 1 ), tExpectedStep, tTol );
    }
}

TEST( Symratio, SymratioMonotone )
{
    belfem::Vector< belfem::real > tX;
    belfem::symratiospace( 0.0, 10.0, 1.5, 11, tX );

    for( size_t i = 1; i < 11; ++i )
    {
        EXPECT_GT( tX( i ), tX( i - 1 ) );
    }
}

TEST( Symratio, SymratioSymmetric )
{
    belfem::Vector< belfem::real > tX;
    belfem::index_t tN = 11;
    belfem::symratiospace( 0.0, 10.0, 1.3, tN, tX );

    for( size_t k = 1; k < tN / 2; ++k )
    {
        belfem::real tLeft  = tX( k ) - tX( k - 1 );
        belfem::real tRight = tX( tN - k ) - tX( tN - k - 1 );
        EXPECT_NEAR( tLeft, tRight, tTol );
    }
}

TEST( Symratio, SymratioSymmetricAboutCenter )
{
    // tX(k) + tX(N-1-k) ≈ xmin + xmax (idea from ChatGPT)
    belfem::Vector< belfem::real > tX;
    belfem::symratiospace( 0.0, 10.0, 2.0, 7, tX );

    for( size_t k = 0; k < 4; ++k )
    {
        EXPECT_NEAR( tX( k ) + tX( 6 - k ), 10.0, tTol );
    }
}

TEST( Symratio, SymratioRegressionBugM1 )
{
    // BUG-M1: midpoint formula uses 0.5*(aXmax - aXmin) instead of
    // aXmin + 0.5*(aXmax - aXmin). For aXmin=2, aXmax=10:
    // correct midpoint is 6.0, buggy gives 4.0.
    belfem::Vector< belfem::real > tX;
    belfem::symratiospace( 2.0, 10.0, 1.0, 5, tX );

    EXPECT_NEAR( tX( 0 ), 2.0, tEps );
    EXPECT_NEAR( tX( 4 ), 10.0, tEps );
    // midpoint: tX(2) should be 6.0
    EXPECT_NEAR( tX( 2 ), 6.0, tTol );
}

TEST( Symratio, SymratioRegressionBugM1Symmetric )
{
    // same setup — verify symmetry around midpoint
    belfem::Vector< belfem::real > tX;
    belfem::symratiospace( 2.0, 10.0, 1.2, 5, tX );

    belfem::real tMid = 0.5 * ( 2.0 + 10.0 );
    // tX(1) - tX(0) should equal tX(4) - tX(3)
    EXPECT_NEAR( tX( 1 ) - tX( 0 ), tX( 4 ) - tX( 3 ), tTol );
    // midpoint should be center
    EXPECT_NEAR( tX( 2 ), tMid, tTol );
}

// --- §9.2 debug ---

#ifndef NDEBUG

TEST( SymratioDebug, SymratioEvenNThrows )
{
    belfem::Vector< belfem::real > tX;
    EXPECT_THROW( belfem::symratiospace( 0.0, 10.0, 1.0, 4, tX ),
                  std::runtime_error );
}

#endif // NDEBUG
