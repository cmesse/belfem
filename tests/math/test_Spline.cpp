/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the spline module: helper polynomials, help matrix,
 * Spline construction, evaluation, boundary conditions, integral, entropy.
 * See: tests_09_spline.md §1–§7
 *
 * The construction sections are gated on BELFEM_SUPERLU, the solver
 * Spline::update_data actually calls; §1, §2, §7 and §8 need no solver
 * at all and always compile.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Spline.hpp"
#include "fn_Create_Truss_Poly.hpp"
#include "fn_Create_Glue_Poly.hpp"

namespace
{
    const belfem::real tEps = 1e-12;
    const belfem::real tTol = 1e-9;

    // Horner evaluation of cubic: a*x³ + b*x² + c*x + d
    belfem::real eval_cubic(
        const belfem::Vector< belfem::real > & aC,
        belfem::real aX )
    {
        return ( ( aC( 0 ) * aX + aC( 1 ) ) * aX + aC( 2 ) ) * aX + aC( 3 );
    }

    // Derivative of cubic: 3a*x² + 2b*x + c
    belfem::real deval_cubic(
        const belfem::Vector< belfem::real > & aC,
        belfem::real aX )
    {
        return ( 3.0 * aC( 0 ) * aX + 2.0 * aC( 1 ) ) * aX + aC( 2 );
    }

    // Horner evaluation of quartic: a*x⁴ + b*x³ + c*x² + d*x + e
    belfem::real eval_quartic(
        const belfem::Vector< belfem::real > & aC,
        belfem::real aX )
    {
        return ( ( ( aC( 0 ) * aX + aC( 1 ) ) * aX + aC( 2 ) ) * aX + aC( 3 ) ) * aX + aC( 4 );
    }

    // Derivative of quartic
    belfem::real deval_quartic(
        const belfem::Vector< belfem::real > & aC,
        belfem::real aX )
    {
        return ( ( 4.0 * aC( 0 ) * aX + 3.0 * aC( 1 ) ) * aX + 2.0 * aC( 2 ) ) * aX + aC( 3 );
    }

#ifdef BELFEM_SUPERLU
    // Generate equidistant x-vector ( used only by the SuperLU-gated
    // spline sections below )
    belfem::Vector< belfem::real > make_x(
        belfem::real aXmin, belfem::real aXmax, belfem::uint aN )
    {
        belfem::Vector< belfem::real > tX( aN );
        belfem::real tDx = ( aXmax - aXmin ) / ( aN - 1 );
        for( belfem::uint k = 0; k < aN; ++k )
        {
            tX( k ) = aXmin + k * tDx;
        }
        return tX;
    }
#endif // BELFEM_SUPERLU
}

// =============================================================================
// §1.1  create_truss_poly  [semantic]
// =============================================================================

TEST( TrussPoly, TrussPolyInterpolatesValues )
{
    // f(x) = x³ on [0, 1]: f(0)=0, f'(0)=0, f(1)=1, f'(1)=3
    belfem::Vector< belfem::real > tX = { 0.0, 1.0 };
    belfem::Vector< belfem::real > tF = { 0.0, 0.0, 1.0, 3.0 };
    belfem::Vector< belfem::real > tC;

    belfem::create_truss_poly( tX, tF, tC );

    EXPECT_NEAR( eval_cubic( tC, 0.0 ), 0.0, tEps );
    EXPECT_NEAR( eval_cubic( tC, 1.0 ), 1.0, tEps );
}

TEST( TrussPoly, TrussPolyInterpolatesDerivatives )
{
    belfem::Vector< belfem::real > tX = { 0.0, 1.0 };
    belfem::Vector< belfem::real > tF = { 0.0, 0.0, 1.0, 3.0 };
    belfem::Vector< belfem::real > tC;

    belfem::create_truss_poly( tX, tF, tC );

    EXPECT_NEAR( deval_cubic( tC, 0.0 ), 0.0, tEps );
    EXPECT_NEAR( deval_cubic( tC, 1.0 ), 3.0, tEps );
}

TEST( TrussPoly, TrussPolyCubicExact )
{
    // f(x) = x³ → coefficients should be {1, 0, 0, 0}
    belfem::Vector< belfem::real > tX = { 0.0, 1.0 };
    belfem::Vector< belfem::real > tF = { 0.0, 0.0, 1.0, 3.0 };
    belfem::Vector< belfem::real > tC;

    belfem::create_truss_poly( tX, tF, tC );

    EXPECT_NEAR( tC( 0 ), 1.0, tEps );
    EXPECT_NEAR( tC( 1 ), 0.0, tEps );
    EXPECT_NEAR( tC( 2 ), 0.0, tEps );
    EXPECT_NEAR( tC( 3 ), 0.0, tEps );
}

TEST( TrussPoly, TrussPolyLinearFunction )
{
    // f(x) = 3x + 1: f(0)=1, f'(0)=3, f(2)=7, f'(2)=3
    belfem::Vector< belfem::real > tX = { 0.0, 2.0 };
    belfem::Vector< belfem::real > tF = { 1.0, 3.0, 7.0, 3.0 };
    belfem::Vector< belfem::real > tC;

    belfem::create_truss_poly( tX, tF, tC );

    EXPECT_NEAR( tC( 0 ), 0.0, tEps );   // cubic = 0
    EXPECT_NEAR( tC( 1 ), 0.0, tEps );   // quadratic = 0
    EXPECT_NEAR( tC( 2 ), 3.0, tEps );   // linear = 3
    EXPECT_NEAR( tC( 3 ), 1.0, tEps );   // constant = 1
}

TEST( TrussPoly, TrussPolyNonOriginInterval )
{
    // f(x) = x² on [2, 5]: f(2)=4, f'(2)=4, f(5)=25, f'(5)=10
    belfem::Vector< belfem::real > tX = { 2.0, 5.0 };
    belfem::Vector< belfem::real > tF = { 4.0, 4.0, 25.0, 10.0 };
    belfem::Vector< belfem::real > tC;

    belfem::create_truss_poly( tX, tF, tC );

    EXPECT_NEAR( eval_cubic( tC, 2.0 ), 4.0, tEps );
    EXPECT_NEAR( eval_cubic( tC, 5.0 ), 25.0, tEps );
    EXPECT_NEAR( deval_cubic( tC, 2.0 ), 4.0, tEps );
    EXPECT_NEAR( deval_cubic( tC, 5.0 ), 10.0, tEps );
}

// =============================================================================
// §1.2  create_glue_poly  [semantic]
// =============================================================================

TEST( GluePoly, GluePolyInterpolatesAllFiveConditions )
{
    // f(x) = x² on [1, 2, 3]:
    // f(1)=1, f'(1)=2, f(2)=4, f(3)=9, f'(3)=6
    belfem::real tX  = 2.0;
    belfem::real tDX = 1.0;
    belfem::Vector< belfem::real > tF = { 1.0, 2.0, 4.0, 9.0, 6.0 };
    belfem::Vector< belfem::real > tC;

    belfem::create_glue_poly( tX, tDX, tF, tC );

    EXPECT_NEAR( eval_quartic( tC, 1.0 ), 1.0, tTol );
    EXPECT_NEAR( deval_quartic( tC, 1.0 ), 2.0, tTol );
    EXPECT_NEAR( eval_quartic( tC, 2.0 ), 4.0, tTol );
    EXPECT_NEAR( eval_quartic( tC, 3.0 ), 9.0, tTol );
    EXPECT_NEAR( deval_quartic( tC, 3.0 ), 6.0, tTol );
}

// =============================================================================
// §2  Help Matrix  [semantic]
// =============================================================================

TEST( HelpMatrix, HelpMatrixSizeCorrect )
{
    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( 10, 1.0, tA );

    EXPECT_EQ( tA.n_rows(), 10u );
    EXPECT_EQ( tA.n_cols(), 10u );
}

TEST( HelpMatrix, HelpMatrixTridiagonalInterior )
{
    belfem::SpMatrix tA;
    belfem::real tDX = 0.5;
    belfem::spline::create_helpmatrix( 10, tDX, tA );

    tA.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::real tC = 1.0 / tDX;
    const belfem::SpMatrix & tConst = tA;

    // interior row k=5: diagonal = 4/dX, off-diag = 1/dX
    EXPECT_NEAR( tConst( 5, 5 ), 4.0 * tC, tEps );
    EXPECT_NEAR( tConst( 4, 5 ), tC, tEps );
    EXPECT_NEAR( tConst( 6, 5 ), tC, tEps );
}

TEST( HelpMatrix, HelpMatrixNoCurvatureBC )
{
    belfem::SpMatrix tA;
    belfem::real tDX = 1.0;
    belfem::spline::create_helpmatrix( 10, tDX, tA,
        belfem::spline::SplineBC::NoCurvature,
        belfem::spline::SplineBC::NoCurvature );

    tA.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    const belfem::SpMatrix & tConst = tA;

    // first row: diag = 2/dX, off-diag = 1/dX
    EXPECT_NEAR( tConst( 0, 0 ), 2.0, tEps );
    EXPECT_NEAR( tConst( 1, 0 ), 1.0, tEps );
}

TEST( HelpMatrix, HelpMatrixTangentBC )
{
    belfem::SpMatrix tA;
    belfem::real tDX = 1.0;
    belfem::spline::create_helpmatrix( 10, tDX, tA,
        belfem::spline::SplineBC::Tangent,
        belfem::spline::SplineBC::NoCurvature );

    tA.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    const belfem::SpMatrix & tConst = tA;

    // tangent BC: first row diag = 1, (0,1) = 0
    EXPECT_NEAR( tConst( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tConst( 0, 1 ), 0.0, tEps );
}

// --- §2.2 debug ---

#ifndef NDEBUG

TEST( HelpMatrixDebug, HelpMatrixTooFewPointsThrows )
{
    belfem::SpMatrix tA;
    EXPECT_THROW(
        belfem::spline::create_helpmatrix( 3, 1.0, tA ),
        std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §3–§7  Spline tests — require the solver Spline::update_data actually uses
//
// This was gated on BELFEM_SUITESPARSE until 2026-08-10, which compiled every
// spline construction/evaluation test OUT of the default build: SuiteSparse is
// off by default ( legacy, not BSD-3 clean ), so spline coverage in the fast
// gate was zero. The gate was also wrong on its own terms — production
// Spline::update_data solves with SUPERLU, not UMFPACK ( cl_Spline.cpp:427-429 ),
// so these tests were gated on a dependency the code under test does not use.
// =============================================================================

#ifdef BELFEM_SUPERLU

// =============================================================================
// §3  Spline Construction and Evaluation  [semantic]
// =============================================================================

TEST( Spline, SplineConstantFunction )
{
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN, 7.0 );

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    // evaluate at midpoints
    for( belfem::uint k = 0; k < tN - 1; ++k )
    {
        belfem::real tMid = 0.5 * ( tX( k ) + tX( k + 1 ) );
        EXPECT_NEAR( tSpline.eval( tMid ), 7.0, tTol );
    }
}

TEST( Spline, SplineLinearFunction )
{
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = 3.0 * tX( k ) + 1.0;
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    belfem::real tMid = 5.0;
    EXPECT_NEAR( tSpline.eval( tMid ), 3.0 * tMid + 1.0, tTol );
    EXPECT_NEAR( tSpline.deval( tMid ), 3.0, tTol );
}

TEST( Spline, SplineCubicClamped )
{
    // f(x) = x³ with clamped BC: f'(x) = 3x² at endpoints
    // cubic spline should recover x³ exactly with clamped BC
    belfem::uint tN = 20;
    belfem::real tXmin = 1.0, tXmax = 5.0;
    belfem::Vector< belfem::real > tX = make_x( tXmin, tXmax, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = std::pow( tX( k ), 3 );
    }

    belfem::real tDx = tX( 1 ) - tX( 0 );
    belfem::real tDydx0 = 3.0 * tXmin * tXmin;   // f'(xmin)
    belfem::real tDydx1 = 3.0 * tXmax * tXmax;   // f'(xmax)

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tDx, tA,
        belfem::spline::SplineBC::Tangent,
        belfem::spline::SplineBC::Tangent );

    belfem::Spline tSpline( tX, tY, tA,
        belfem::spline::SplineBC::Tangent,
        belfem::spline::SplineBC::Tangent,
        tDydx0, tDydx1 );

    // test at interior points
    belfem::real tTestX = 3.0;
    EXPECT_NEAR( tSpline.eval( tTestX ), 27.0, tTol );
    EXPECT_NEAR( tSpline.deval( tTestX ), 27.0, tTol );
    EXPECT_NEAR( tSpline.ddeval( tTestX ), 18.0, tTol );
}

// =============================================================================
// §3.2  Knot Interpolation  [semantic]
// =============================================================================

TEST( Spline, SplineInterpolatesKnots )
{
    belfem::uint tN = 15;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 7.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = std::sin( tX( k ) );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    for( belfem::uint k = 0; k < tN; ++k )
    {
        EXPECT_NEAR( tSpline.eval( tX( k ) ), tY( k ), tTol );
    }
}

TEST( Spline, SplineContinuousAcrossKnots )
{
    belfem::uint tN = 15;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 7.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = std::sin( tX( k ) );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    belfem::real tEpsX = 1e-8;
    for( belfem::uint k = 1; k < tN - 1; ++k )
    {
        belfem::real tLeft  = tSpline.eval( tX( k ) - tEpsX );
        belfem::real tRight = tSpline.eval( tX( k ) + tEpsX );
        EXPECT_NEAR( tLeft, tRight, 1e-6 );
    }
}

// =============================================================================
// §3.3  Boundary Conditions  [semantic]
// =============================================================================

TEST( Spline, NoCurvatureBC )
{
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = std::sin( tX( k ) );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    // natural spline: second derivative ≈ 0 at endpoints
    EXPECT_NEAR( tSpline.ddeval( tX( 0 ) ), 0.0, tTol );
    EXPECT_NEAR( tSpline.ddeval( tX( tN - 1 ) ), 0.0, tTol );
}

// =============================================================================
// §3.4  Extrapolation (NOTE-S1)  [semantic]
// =============================================================================

// Design note: out-of-range evaluation is intentional.
//
// Spline::find_col() selects the boundary interval for x outside [x_min, x_max],
// but eval()/deval()/ddeval() still use the original x value in the polynomial.
// This means the class EXTRAPOLATES with the first/last interval polynomial
// instead of throwing or hard-clamping to the endpoint value.

TEST( Spline, EvalBelowRangeUsesFirstIntervalExtrapolation )
{
    // (regression pattern from ChatGPT)
    // eval(x) with x < x_min should give the same result as
    // eval(x, 0) — explicit first-interval polynomial evaluation.
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 1.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = tX( k ) * tX( k );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    belfem::real tXout = tSpline.x_min() - 0.25 * tSpline.delta_x();
    belfem::real tExpected = tSpline.eval( tXout, 0 );   // explicit first-interval poly
    EXPECT_NEAR( tSpline.eval( tXout ), tExpected, tEps );
}

TEST( Spline, EvalAboveRangeUsesLastIntervalExtrapolation )
{
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 1.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = tX( k ) * tX( k );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    belfem::real tXout = tSpline.x_max() + 0.25 * tSpline.delta_x();
    // last interval is column n-2 (last column n-1 duplicates n-2)
    belfem::real tExpected = tSpline.eval( tXout, tSpline.n() - 2 );
    EXPECT_NEAR( tSpline.eval( tXout ), tExpected, tEps );
}

// =============================================================================
// §3.5  Accessors  [semantic]
// =============================================================================

TEST( Spline, SplineAccessors )
{
    belfem::uint tN = 20;
    belfem::real tXmin = 1.0, tXmax = 10.0;
    belfem::Vector< belfem::real > tX = make_x( tXmin, tXmax, tN );
    belfem::Vector< belfem::real > tY( tN, 1.0 );

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    EXPECT_EQ( tSpline.n(), tN );
    EXPECT_NEAR( tSpline.x_min(), tXmin, tEps );
    EXPECT_NEAR( tSpline.x_max(), tXmax, tEps );
    EXPECT_NEAR( tSpline.delta_x(), ( tXmax - tXmin ) / ( tN - 1 ), tEps );
    EXPECT_EQ( tSpline.coefficients().n_rows(), 4u );
    EXPECT_EQ( tSpline.coefficients().n_cols(), tN );
}

// =============================================================================
// §4  Integral Mode  [semantic]
// =============================================================================

TEST( SplineIntegral, IntegralOfLinear )
{
    // f(x) = 2x + 1, integral from 0 to L = L² + L
    belfem::uint tN = 20;
    belfem::real tL = 5.0;
    belfem::Vector< belfem::real > tX = make_x( 0.0, tL, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = 2.0 * tX( k ) + 1.0;
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );
    tSpline.create_integral();

    belfem::real tExpected = tL * tL + tL;   // L² + L
    EXPECT_NEAR( tSpline.integrate( 0.0, tL ), tExpected, tTol );
}

TEST( SplineIntegral, IntegralContinuousAcrossKnots )
{
    // (idea from Gemini) — verify integration continuity at interior knots
    belfem::uint tN = 15;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 7.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = std::sin( tX( k ) );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );
    tSpline.create_integral();

    for( belfem::uint k = 2; k < tN - 2; ++k )
    {
        belfem::real tLeft  = tSpline.integrate( tX( k ) - 1e-8 );
        belfem::real tRight = tSpline.integrate( tX( k ) + 1e-8 );
        EXPECT_NEAR( tLeft, tRight, 1e-5 );
    }
}

// --- §4.2 debug ---

#ifndef NDEBUG

TEST( SplineIntegralDebug, IntegrateWithoutCreateThrows )
{
    belfem::uint tN = 10;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 5.0, tN );
    belfem::Vector< belfem::real > tY( tN, 1.0 );

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    // integrate without create_integral → assertion
    EXPECT_THROW( tSpline.integrate( 1.0 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §5  Entropy Mode  [semantic]
// =============================================================================

TEST( SplineEntropy, EntropyModeConstruction )
{
    // Use a simple function — entropy mode requires aXref > 0
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 1.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = tX( k ) * tX( k );   // cp-like function
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::real tXref = 5.0;
    belfem::real tSref = 0.0;

    belfem::Spline tSpline( tX, tY, tA, tXref, tSref );

    // entropy should be callable without assertion
    belfem::real tS = tSpline.entropy( 5.0 );
    EXPECT_NEAR( tS, tSref, tTol );

    // dentropy should be callable
    belfem::real tDS = tSpline.dentropy( 5.0 );
    ( void ) tDS;
}

TEST( SplineEntropy, DentropyConsistentWithFiniteDifference )
{
    // (idea from Gemini) — verify dentropy matches finite difference of entropy
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 1.0, 10.0, tN );
    belfem::Vector< belfem::real > tY( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY( k ) = tX( k ) * tX( k );
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 5.0, 0.0 );

    belfem::real tTestX = 4.0;
    belfem::real tH = 1e-5;
    belfem::real tFD = ( tSpline.entropy( tTestX + tH )
                       - tSpline.entropy( tTestX - tH ) ) / ( 2.0 * tH );

    EXPECT_NEAR( tSpline.dentropy( tTestX ), tFD, 1e-4 );
}

// --- §5.2 debug ---

#ifndef NDEBUG

TEST( SplineEntropyDebug, EntropyWithoutCreateThrows )
{
    belfem::uint tN = 10;
    belfem::Vector< belfem::real > tX = make_x( 1.0, 5.0, tN );
    belfem::Vector< belfem::real > tY( tN, 1.0 );

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    // plain spline — no entropy (aXref=0 → entropy not created)
    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    EXPECT_THROW( tSpline.entropy( 2.0 ), std::runtime_error );
}

TEST( SplineEntropyDebug, DentropyWithoutCreateThrows )
{
    belfem::uint tN = 10;
    belfem::Vector< belfem::real > tX = make_x( 1.0, 5.0, tN );
    belfem::Vector< belfem::real > tY( tN, 1.0 );

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    EXPECT_THROW( tSpline.dentropy( 2.0 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §6  update_data  [semantic]
// =============================================================================

TEST( SplineUpdate, UpdateDataChangesCoefficients )
{
    // Codex finding: use exactly recoverable target (constant) for strong check
    belfem::uint tN = 20;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 10.0, tN );
    belfem::Vector< belfem::real > tY1( tN );
    belfem::Vector< belfem::real > tY2( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tY1( k ) = tX( k );   // linear: f(x) = x
        tY2( k ) = 7.0;       // constant: f(x) = 7 (exactly recoverable)
    }

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY1, tA, 0.0 );

    // verify linear
    EXPECT_NEAR( tSpline.eval( 5.0 ), 5.0, tTol );

    // update with constant data
    tSpline.update_data( tA, tY2 );

    // constant is exactly recoverable with any BCs
    EXPECT_NEAR( tSpline.eval( 5.0 ), 7.0, tTol );
    EXPECT_NEAR( tSpline.eval( 2.5 ), 7.0, tTol );
}

TEST( SplineUpdate, UpdateDataPreservesGrid )
{
    belfem::uint tN = 20;
    belfem::real tXmin = 1.0, tXmax = 10.0;
    belfem::Vector< belfem::real > tX = make_x( tXmin, tXmax, tN );
    belfem::Vector< belfem::real > tY( tN, 1.0 );

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tA );

    belfem::Spline tSpline( tX, tY, tA, 0.0 );

    belfem::Vector< belfem::real > tY2( tN, 2.0 );
    tSpline.update_data( tA, tY2 );

    EXPECT_NEAR( tSpline.x_min(), tXmin, tEps );
    EXPECT_NEAR( tSpline.x_max(), tXmax, tEps );
    EXPECT_NEAR( tSpline.delta_x(), ( tXmax - tXmin ) / ( tN - 1 ), tEps );
}

// --- §7.2 debug (CheckInputTooFewPointsThrows is already tested
//     as HelpMatrixTooFewPointsThrows above, outside SUITESPARSE guard) ---

#endif // BELFEM_SUPERLU

// =============================================================================
// §7  Construction Variants — no solver needed
// =============================================================================

TEST( SplineVariant, EmptyContainerConstructor )
{
    belfem::Spline tSpline( 50, 0.0, 10.0 );

    EXPECT_EQ( tSpline.n(), 50u );
    EXPECT_NEAR( tSpline.x_min(), 0.0, tEps );
    EXPECT_NEAR( tSpline.x_max(), 10.0, tEps );
}

// =============================================================================
// §8  Explicit column overloads — no solver needed
// =============================================================================

// eval, deval, ddeval, entropy and dentropy each have a second form taking
// the coefficient column explicitly. That form skips find_col() and must
// therefore use the column it is given and no other — also for an x that
// lies outside that interval, which is exactly what makes it usable as the
// reference in §3.4 and as the fast path for a caller that has already
// located its interval.
//
// The contract has two halves, and only the first is testable without
// trusting the class:
//
//   1. eval( x, k ) is the polynomial of column k evaluated at x. Checked
//      here against coefficients written by hand through matrix_data(),
//      so the expected value is known in closed form.
//   2. eval( x ) is eval( x, find_col( x ) ). Checked as an exact identity,
//      since both sides are the same arithmetic in the same order.
//
// Writing the table directly keeps the solver out of the loop, so these
// tests stay outside the SuperLU gate above and run in every configuration.

namespace
{
    const belfem::index_t tColN = 11;         // 11 points -> 10 intervals
    const belfem::real    tColXmin = 1.0;     // > 0, entropy takes log( x )
    const belfem::real    tColXmax = 11.0;

    // coefficients of column k, all five rows distinct per column so that
    // reading the neighboring column cannot pass by coincidence
    belfem::real col_a( const belfem::index_t aK ) { return  0.5  + 0.5  * aK; }
    belfem::real col_b( const belfem::index_t aK ) { return -1.0  - 0.1  * aK; }
    belfem::real col_c( const belfem::index_t aK ) { return  2.0  + 0.25 * aK; }
    belfem::real col_d( const belfem::index_t aK ) { return -3.0  + 1.0  * aK; }
    belfem::real col_e( const belfem::index_t aK ) { return  0.75 + 0.5  * aK; }

    // builds a spline with a hand-written coefficient table; five rows, so
    // the entropy row is present and aEntropy can switch the mode on
    void make_column_spline( belfem::Spline & aSpline, const bool aEntropy )
    {
        belfem::Matrix< belfem::real > & tData = aSpline.matrix_data();
        tData.set_size( 5, tColN, 0.0 );

        for( belfem::index_t k = 0; k < tColN; ++k )
        {
            tData( 0, k ) = col_a( k );
            tData( 1, k ) = col_b( k );
            tData( 2, k ) = col_c( k );
            tData( 3, k ) = col_d( k );
            tData( 4, k ) = col_e( k );
        }

        if( aEntropy )
        {
            aSpline.set_extra_mode( belfem::spline::ExtraMode::Entropy );
        }
    }

    belfem::real ref_eval( const belfem::index_t aK, const belfem::real aX )
    {
        return ( ( col_a( aK ) * aX + col_b( aK ) ) * aX + col_c( aK ) ) * aX
               + col_d( aK );
    }

    belfem::real ref_deval( const belfem::index_t aK, const belfem::real aX )
    {
        return ( 3.0 * col_a( aK ) * aX + 2.0 * col_b( aK ) ) * aX + col_c( aK );
    }

    belfem::real ref_ddeval( const belfem::index_t aK, const belfem::real aX )
    {
        return 6.0 * col_a( aK ) * aX + 2.0 * col_b( aK );
    }

    belfem::real ref_entropy( const belfem::index_t aK, const belfem::real aX )
    {
        return ( 1.5 * col_a( aK ) * aX + 2.0 * col_b( aK ) ) * aX
               + col_c( aK ) * std::log( aX ) + col_e( aK );
    }

    belfem::real ref_dentropy( const belfem::index_t aK, const belfem::real aX )
    {
        return 3.0 * col_a( aK ) * aX + 2.0 * col_b( aK ) + col_c( aK ) / aX;
    }
}

TEST( SplineColumn, EvalWithColumnUsesGivenColumn )
{
    belfem::Spline tSpline( tColN, tColXmin, tColXmax );
    make_column_spline( tSpline, false );

    // every valid column, evaluated at three x: the left knot of the
    // interval, its midpoint, and a point far outside it. The last one is
    // the case that separates "uses column k" from "clamps like find_col".
    for( belfem::index_t k = 0; k < tColN - 1; ++k )
    {
        const belfem::real tXk = tColXmin + k;

        const belfem::real tX[ 3 ] = { tXk, tXk + 0.5, tColXmax + 4.0 };

        for( belfem::uint j = 0; j < 3; ++j )
        {
            EXPECT_NEAR( tSpline.eval(   tX[ j ], k ), ref_eval(   k, tX[ j ] ), tTol );
            EXPECT_NEAR( tSpline.deval(  tX[ j ], k ), ref_deval(  k, tX[ j ] ), tTol );
            EXPECT_NEAR( tSpline.ddeval( tX[ j ], k ), ref_ddeval( k, tX[ j ] ), tTol );
        }
    }
}

TEST( SplineColumn, EntropyWithColumnUsesGivenColumn )
{
    belfem::Spline tSpline( tColN, tColXmin, tColXmax );
    make_column_spline( tSpline, true );

    for( belfem::index_t k = 0; k < tColN - 1; ++k )
    {
        const belfem::real tXk = tColXmin + k;

        const belfem::real tX[ 3 ] = { tXk, tXk + 0.5, tColXmax + 4.0 };

        for( belfem::uint j = 0; j < 3; ++j )
        {
            EXPECT_NEAR( tSpline.entropy(  tX[ j ], k ), ref_entropy(  k, tX[ j ] ), tTol );
            EXPECT_NEAR( tSpline.dentropy( tX[ j ], k ), ref_dentropy( k, tX[ j ] ), tTol );
        }
    }
}

TEST( SplineColumn, ColumnOverloadsAgreeWithFindCol )
{
    belfem::Spline tSpline( tColN, tColXmin, tColXmax );
    make_column_spline( tSpline, true );

    // the column-free forms are defined as find_col() followed by the
    // column form, so the two must agree bit for bit, not just to a
    // tolerance — including below x_min and above x_max, where find_col
    // clamps to the boundary interval but the polynomial still sees the
    // original x ( the extrapolation of §3.4 )
    // every x stays positive — entropy takes log( x ), so the below-range
    // sample is 0.25 rather than a negative abscissa
    const belfem::real tX[ 7 ] =
        { 0.25, tColXmin, tColXmin + 0.25, 5.5,
          tColXmax - 1e-9, tColXmax, tColXmax + 2.0 };

    for( belfem::uint j = 0; j < 7; ++j )
    {
        const belfem::index_t tCol = tSpline.find_col( tX[ j ] );

        EXPECT_LT( tCol, tSpline.n() - 1 );

        EXPECT_DOUBLE_EQ( tSpline.eval(     tX[ j ] ), tSpline.eval(     tX[ j ], tCol ) );
        EXPECT_DOUBLE_EQ( tSpline.deval(    tX[ j ] ), tSpline.deval(    tX[ j ], tCol ) );
        EXPECT_DOUBLE_EQ( tSpline.ddeval(   tX[ j ] ), tSpline.ddeval(   tX[ j ], tCol ) );
        EXPECT_DOUBLE_EQ( tSpline.entropy(  tX[ j ] ), tSpline.entropy(  tX[ j ], tCol ) );
        EXPECT_DOUBLE_EQ( tSpline.dentropy( tX[ j ] ), tSpline.dentropy( tX[ j ], tCol ) );
    }
}

// --- §8.1 debug ---

#ifndef NDEBUG

TEST( SplineColumnDebug, ColumnOutOfRangeThrows )
{
    belfem::Spline tSpline( tColN, tColXmin, tColXmax );
    make_column_spline( tSpline, true );

    // the coefficient table has n() columns but only n()-1 intervals; the
    // last column is a duplicate that find_col never returns, so the valid
    // range ends one short of the table width
    const belfem::index_t tBad = tSpline.n() - 1;

    EXPECT_EQ( tSpline.matrix_data().n_cols(), tColN );

    EXPECT_THROW( tSpline.eval(     2.0, tBad ), std::runtime_error );
    EXPECT_THROW( tSpline.deval(    2.0, tBad ), std::runtime_error );
    EXPECT_THROW( tSpline.ddeval(   2.0, tBad ), std::runtime_error );
    EXPECT_THROW( tSpline.entropy(  2.0, tBad ), std::runtime_error );
    EXPECT_THROW( tSpline.dentropy( 2.0, tBad ), std::runtime_error );
}

TEST( SplineColumnDebug, EntropyColumnWithoutEntropyModeThrows )
{
    belfem::Spline tSpline( tColN, tColXmin, tColXmax );
    make_column_spline( tSpline, false );   // extra mode stays None

    EXPECT_THROW( tSpline.entropy(  2.0, 0 ), std::runtime_error );
    EXPECT_THROW( tSpline.dentropy( 2.0, 0 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// Input Validation  [debug]
// =============================================================================

#ifdef BELFEM_SUPERLU
#ifndef NDEBUG

TEST( SplineDebug, TooFewPointsThrows )
{
    // check_input asserts aX.length() > 3 (need at least 4 points)
    // create_helpmatrix also asserts this, so test it directly
    belfem::SpMatrix tA;
    EXPECT_THROW( belfem::spline::create_helpmatrix( 3, 1.0, tA ), std::runtime_error );
}

TEST( SplineDebug, LengthMismatchThrows )
{
    belfem::uint tN = 10;
    belfem::Vector< belfem::real > tX = make_x( 0.0, 9.0, tN );
    belfem::Vector< belfem::real > tY( tN + 1, 1.0 );   // wrong length

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( tN, 1.0, tA );

    EXPECT_THROW( belfem::Spline( tX, tY, tA, 0.0 ), std::runtime_error );
}

TEST( SplineDebug, NonEquidistantThrows )
{
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 3.0, 4.0, 5.0 };  // gap at index 2
    belfem::Vector< belfem::real > tY = { 0.0, 1.0, 9.0, 16.0, 25.0 };

    belfem::SpMatrix tA;
    belfem::spline::create_helpmatrix( 5, 1.0, tA );

    EXPECT_THROW( belfem::Spline( tX, tY, tA, 0.0 ), std::runtime_error );
}

#endif // NDEBUG
#endif // BELFEM_SUPERLU
