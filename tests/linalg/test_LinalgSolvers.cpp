/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for linear algebra solvers, polynomial functions, and
 * statistical functions.
 * See: tests_02_linalg.md §4.4 (crossmat), §5 (gesv/posv), §6 (polyfit/polyval/r2/eigen)
 *
 * Solver tests verify via residual norm — NOT exact solution comparison.
 * gesv and posv mutate their inputs; copies of A and b are saved before solving.
 *
 * Includes BUG-L2 regression test for r2() matrix version.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"
#include "fn_posv.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "fn_r2.hpp"
#include "fn_eigen.hpp"
#include "fn_norm.hpp"
#include "fn_crossmat.hpp"
#include "op_MatrixTimes.hpp"
#include "op_VectorMinus.hpp"

namespace
{
    const belfem::real tEps = 1e-12;  // exact-in-theory results
    const belfem::real tTol = 1e-9;   // solver residuals, polyfit
}

// =============================================================================
// §5.1  gesv (Vector RHS)  [semantic]
// =============================================================================

TEST( LinalgSolver, GesvIdentitySystem )
{
    // I * x = b  →  x == b
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 1.0 } };
    belfem::Vector< belfem::real > tX = { 3.0, 5.0, 7.0 };
    belfem::Vector< belfem::int_t >          tP( 3 );

    // save originals for residual check
    belfem::Vector< belfem::real > tB( tX );   // b = original RHS
    belfem::Matrix< belfem::real > tAorig( tA );

    belfem::gesv( tA, tX, tP );

    // for identity system, solution equals RHS exactly
    EXPECT_NEAR( tX( 0 ), 3.0, tEps );
    EXPECT_NEAR( tX( 1 ), 5.0, tEps );
    EXPECT_NEAR( tX( 2 ), 7.0, tEps );
}

TEST( LinalgSolver, GesvKnown2x2 )
{
    // {{2,1},{1,3}} * x = {5,7}  →  x = {1.6, 1.8}
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 },
                                           { 1.0, 3.0 } };
    belfem::Vector< belfem::real > tX = { 5.0, 7.0 };
    belfem::Vector< belfem::int_t >          tP( 2 );

    // save copies
    belfem::Matrix< belfem::real > tAorig( tA );
    belfem::Vector< belfem::real > tB( tX );

    belfem::gesv( tA, tX, tP );

    // verify via residual: ||A_orig * x - b|| < tTol
    belfem::Vector< belfem::real > tResidual( tAorig * tX - tB );
    EXPECT_NEAR( belfem::norm( tResidual ), 0.0, tTol );
}

TEST( LinalgSolver, GesvKnown3x3 )
{
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0, -1.0 },
                                           { -3.0, -1.0, 2.0 },
                                           { -2.0, 1.0, 2.0 } };
    belfem::Vector< belfem::real > tX = { 8.0, -11.0, -3.0 };
    belfem::Vector< belfem::int_t >          tP( 3 );

    belfem::Matrix< belfem::real > tAorig( tA );
    belfem::Vector< belfem::real > tB( tX );

    belfem::gesv( tA, tX, tP );

    belfem::Vector< belfem::real > tResidual( tAorig * tX - tB );
    EXPECT_NEAR( belfem::norm( tResidual ), 0.0, tTol );
}

TEST( LinalgSolver, GesvMutatesMatrix )
{
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 }, { 1.0, 3.0 } };
    belfem::Matrix< belfem::real > tAorig( tA );
    belfem::Vector< belfem::real > tX = { 5.0, 7.0 };
    belfem::Vector< belfem::int_t >          tP( 2 );

    belfem::gesv( tA, tX, tP );

    // A should be overwritten (LU factors) — at least one entry differs
    bool tChanged = false;
    for( size_t i = 0; i < 2; ++i )
    {
        for( size_t j = 0; j < 2; ++j )
        {
            if( std::abs( tA( i, j ) - tAorig( i, j ) ) > tEps )
            {
                tChanged = true;
            }
        }
    }
    EXPECT_TRUE( tChanged );
}

TEST( LinalgSolver, GesvMutatesRhs )
{
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 }, { 1.0, 3.0 } };
    belfem::Vector< belfem::real > tX = { 5.0, 7.0 };
    belfem::Vector< belfem::real > tBoriginal( tX );
    belfem::Vector< belfem::int_t >          tP( 2 );

    belfem::gesv( tA, tX, tP );

    // tX now contains the solution, not the original RHS
    bool tChanged = false;
    for( size_t i = 0; i < 2; ++i )
    {
        if( std::abs( tX( i ) - tBoriginal( i ) ) > tEps )
        {
            tChanged = true;
        }
    }
    EXPECT_TRUE( tChanged );
}

TEST( LinalgSolver, GesvPivotVectorPopulated )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0, 3.0 },
                                           { 4.0, 5.0, 6.0 },
                                           { 7.0, 8.0, 10.0 } };
    belfem::Vector< belfem::real > tX = { 1.0, 1.0, 1.0 };
    belfem::Vector< belfem::int_t >          tP( 3, 0 );

    belfem::gesv( tA, tX, tP );

    // pivot vector should have been written — at least one entry != 0
    bool tPopulated = false;
    for( size_t i = 0; i < 3; ++i )
    {
        if( tP( i ) != 0 )
        {
            tPopulated = true;
        }
    }
    EXPECT_TRUE( tPopulated );
}

// =============================================================================
// §5.2  gesv (Matrix RHS)  [semantic]
// =============================================================================

TEST( LinalgSolver, GesvMultipleRhs )
{
    // Solve A * X = B where B is 3×2 (two right-hand sides)
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0, 0.0 },
                                           { 1.0, 3.0, 1.0 },
                                           { 0.0, 1.0, 2.0 } };
    belfem::Matrix< belfem::real > tX = { { 1.0, 4.0 },
                                           { 2.0, 5.0 },
                                           { 3.0, 6.0 } };
    belfem::Vector< belfem::int_t >          tP( 3 );

    belfem::Matrix< belfem::real > tAorig( tA );
    belfem::Matrix< belfem::real > tB( tX );

    belfem::gesv( tA, tX, tP );

    // verify each column of the solution via residual
    for( size_t j = 0; j < 2; ++j )
    {
        belfem::Vector< belfem::real > tXcol( tX.col( j ) );
        belfem::Vector< belfem::real > tBcol( tB.col( j ) );
        belfem::Vector< belfem::real > tRes( tAorig * tXcol - tBcol );
        EXPECT_NEAR( belfem::norm( tRes ), 0.0, tTol );
    }
}

// =============================================================================
// §5.3  gesv  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( LinalgSolverDebug, GesvRowMismatchThrows )
{
    belfem::Matrix< belfem::real > tA( 3, 3, 1.0 );
    belfem::Vector< belfem::real > tX( 4, 1.0 );   // wrong: 4 != 3
    belfem::Vector< belfem::int_t >          tP( 4 );

    EXPECT_THROW( belfem::gesv( tA, tX, tP ), std::runtime_error );
}

TEST( LinalgSolverDebug, GesvColMismatchThrows )
{
    // non-square A → n_cols != n_rows → second BELFEM_ASSERT fails
    belfem::Matrix< belfem::real > tA( 3, 4, 1.0 );
    belfem::Vector< belfem::real > tX( 3, 1.0 );
    belfem::Vector< belfem::int_t >          tP( 3 );

    EXPECT_THROW( belfem::gesv( tA, tX, tP ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §5.4  posv  [semantic]
// =============================================================================

TEST( LinalgSolver, PosvSpdSystem )
{
    // symmetric positive definite: A = {{4,2},{2,3}}
    belfem::Matrix< belfem::real > tA = { { 4.0, 2.0 },
                                           { 2.0, 3.0 } };
    belfem::Vector< belfem::real > tX = { 1.0, 2.0 };

    belfem::Matrix< belfem::real > tAorig( tA );
    belfem::Vector< belfem::real > tB( tX );

    belfem::posv( tA, tX );

    // verify via residual
    belfem::Vector< belfem::real > tResidual( tAorig * tX - tB );
    EXPECT_NEAR( belfem::norm( tResidual ), 0.0, tTol );
}

// =============================================================================
// §4.4  Crossmat  [semantic]
// =============================================================================

TEST( Crossmat, Crossmat2DBasic )
{
    // 2D: n × A, where n is length-2 normal, A is 2×N matrix
    belfem::Vector< belfem::real > tN = { 1.0, 0.0 };
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0 },
                                           { 0.0, 1.0 } };

    belfem::Vector< belfem::real > tResult( 2, 0.0 );
    belfem::crossmat( tN, tA, tResult );

    // n(0)*A(1,k) - n(1)*A(0,k) = 1*{0,1} - 0*{1,0} = {0,1}
    EXPECT_NEAR( tResult( 0 ), 0.0, tEps );
    EXPECT_NEAR( tResult( 1 ), 1.0, tEps );
}

TEST( Crossmat, Crossmat2DWithScale )
{
    belfem::Vector< belfem::real > tN = { 1.0, 0.0 };
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0 },
                                           { 0.0, 1.0 } };

    belfem::Vector< belfem::real > tResult( 2, 0.0 );
    belfem::crossmat( tN, tA, 2.0, tResult );

    // scale factor 2.0 applied, accumulated onto initial 0.0
    EXPECT_NEAR( tResult( 0 ), 0.0, tEps );
    EXPECT_NEAR( tResult( 1 ), 2.0, tEps );
}

TEST( Crossmat, Crossmat3DBasic )
{
    // 3D: n × A, where n is length-3 normal, A is 3×N matrix
    // use n = e_x = {1,0,0}, A = I_3
    belfem::Vector< belfem::real > tN = { 1.0, 0.0, 0.0 };
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 1.0 } };

    belfem::Matrix< belfem::real > tResult( 3, 3, 0.0 );
    belfem::crossmat( tN, tA, tResult );

    // e_x × e_x = 0, e_x × e_y = e_z, e_x × e_z = -e_y
    // column 0: n × A(:,0) = {1,0,0} × {1,0,0} = {0,0,0}
    EXPECT_NEAR( tResult( 0, 0 ), 0.0, tEps );
    EXPECT_NEAR( tResult( 1, 0 ), 0.0, tEps );
    EXPECT_NEAR( tResult( 2, 0 ), 0.0, tEps );

    // column 1: n × A(:,1) = {1,0,0} × {0,1,0} = {0,0,1}
    EXPECT_NEAR( tResult( 0, 1 ), 0.0, tEps );
    EXPECT_NEAR( tResult( 1, 1 ), 0.0, tEps );
    EXPECT_NEAR( tResult( 2, 1 ), 1.0, tEps );

    // column 2: n × A(:,2) = {1,0,0} × {0,0,1} = {0,-1,0}
    EXPECT_NEAR( tResult( 0, 2 ), 0.0,  tEps );
    EXPECT_NEAR( tResult( 1, 2 ), -1.0, tEps );
    EXPECT_NEAR( tResult( 2, 2 ), 0.0,  tEps );
}

TEST( Crossmat, Crossmat3DWithScale )
{
    belfem::Vector< belfem::real > tN = { 1.0, 0.0, 0.0 };
    belfem::Matrix< belfem::real > tA = { { 0.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 0.0 } };

    // pre-fill with non-zero to verify accumulation (Codex finding:
    // the scaled overload is additive, so result = existing + scale*(n×A))
    belfem::Matrix< belfem::real > tResult( 3, 3, 1.0 );
    belfem::crossmat( tN, tA, 3.0, tResult );

    // column 1: n × {0,1,0} = {0,0,1}, scaled by 3 → {0,0,3}, accumulated onto 1.0
    EXPECT_NEAR( tResult( 2, 1 ), 4.0, tEps );   // 1.0 + 3.0
    // zero cross-product entries still accumulated onto 1.0
    EXPECT_NEAR( tResult( 0, 0 ), 1.0, tEps );    // 1.0 + 0.0
    EXPECT_NEAR( tResult( 1, 1 ), 1.0, tEps );    // 1.0 + 0.0
}

TEST( Crossmat, Crossmat2DDustRemoval )
{
    // near-zero entries relative to norm should be zeroed.
    // Use a 2-column matrix where one result dominates so that
    // the tiny entry is small relative to the result norm.
    belfem::Vector< belfem::real > tN = { 1.0, 0.0 };
    belfem::Matrix< belfem::real > tA = { { 0.0,   1.0 },       // row 0
                                           { 1.0,   1e-18 } };   // row 1

    // result(k) = n(0)*A(1,k) - n(1)*A(0,k)
    // result(0) = 1*1.0 - 0*0.0 = 1.0  (dominant)
    // result(1) = 1*1e-18 - 0*1.0 = 1e-18  (dust)
    // norm ≈ 1.0, so |1e-18 / 1.0| < BELFEM_EPSILON → zeroed
    belfem::Vector< belfem::real > tResult( 2, 0.0 );
    belfem::crossmat( tN, tA, tResult );

    EXPECT_NEAR( tResult( 0 ), 1.0, tEps );   // dominant entry preserved
    EXPECT_NEAR( tResult( 1 ), 0.0, tEps );   // dust removed
}

// =============================================================================
// §6.1  Polyval  [semantic]
// =============================================================================

TEST( Polynomial, PolyvalConstant )
{
    // coeffs {5.0} → p(x) = 5 for any x
    belfem::Vector< belfem::real > tCoeffs = { 5.0 };
    EXPECT_NEAR( belfem::polyval( tCoeffs, 0.0 ),   5.0, tEps );
    EXPECT_NEAR( belfem::polyval( tCoeffs, 100.0 ), 5.0, tEps );
    EXPECT_NEAR( belfem::polyval( tCoeffs, -7.0 ),  5.0, tEps );
}

TEST( Polynomial, PolyvalLinear )
{
    // coeffs {2.0, 3.0} → p(x) = 2x + 3
    belfem::Vector< belfem::real > tCoeffs = { 2.0, 3.0 };
    EXPECT_NEAR( belfem::polyval( tCoeffs, 4.0 ), 11.0, tEps );
    EXPECT_NEAR( belfem::polyval( tCoeffs, 0.0 ), 3.0,  tEps );
}

TEST( Polynomial, PolyvalQuadratic )
{
    // coeffs {1.0, 0.0, -1.0} → p(x) = x² - 1
    belfem::Vector< belfem::real > tCoeffs = { 1.0, 0.0, -1.0 };
    EXPECT_NEAR( belfem::polyval( tCoeffs, 3.0 ), 8.0,  tEps );
    EXPECT_NEAR( belfem::polyval( tCoeffs, 0.0 ), -1.0, tEps );
    EXPECT_NEAR( belfem::polyval( tCoeffs, 1.0 ), 0.0,  tEps );
}

TEST( Polynomial, PolyvalVectorized )
{
    // vectorized overload: polyval(coeffs, xVec, yVec)
    belfem::Vector< belfem::real > tCoeffs = { 1.0, -2.0, 1.0 };  // x² - 2x + 1
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tY;

    belfem::polyval( tCoeffs, tX, tY );

    EXPECT_EQ( tY.length(), 4u );
    for( size_t i = 0; i < tX.length(); ++i )
    {
        belfem::real tExpected = belfem::polyval( tCoeffs, tX( i ) );
        EXPECT_NEAR( tY( i ), tExpected, tEps );
    }
}

// =============================================================================
// §6.2  Dpolyval (first derivative)  [semantic]
// =============================================================================

TEST( Polynomial, DpolyvalLinear )
{
    // d/dx(2x + 3) = 2 for any x
    belfem::Vector< belfem::real > tCoeffs = { 2.0, 3.0 };
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, 0.0 ), 2.0, tEps );
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, 5.0 ), 2.0, tEps );
}

TEST( Polynomial, DpolyvalQuadratic )
{
    // d/dx(x² - 1) = 2x  →  at x=3: 6.0
    belfem::Vector< belfem::real > tCoeffs = { 1.0, 0.0, -1.0 };
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, 3.0 ), 6.0, tEps );
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, 0.0 ), 0.0, tEps );
}

TEST( Polynomial, DpolyvalConstantIsZero )
{
    // d/dx(5) = 0
    belfem::Vector< belfem::real > tCoeffs = { 5.0 };
    EXPECT_NEAR( belfem::dpolyval( tCoeffs, 42.0 ), 0.0, tEps );
}

// =============================================================================
// §6.3  DDpolyval (second derivative)  [semantic]
// =============================================================================

TEST( Polynomial, DDpolyvalQuadratic )
{
    // d²/dx²(x² - 1) = 2
    belfem::Vector< belfem::real > tCoeffs = { 1.0, 0.0, -1.0 };
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, 0.0 ),  2.0, tEps );
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, 99.0 ), 2.0, tEps );
}

TEST( Polynomial, DDpolyvalLinearIsZero )
{
    // d²/dx²(2x + 3) = 0
    belfem::Vector< belfem::real > tCoeffs = { 2.0, 3.0 };
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, 5.0 ), 0.0, tEps );
}

TEST( Polynomial, DDpolyvalConstantIsZero )
{
    belfem::Vector< belfem::real > tCoeffs = { 7.0 };
    EXPECT_NEAR( belfem::ddpolyval( tCoeffs, 5.0 ), 0.0, tEps );
}

// =============================================================================
// §6.4  Polyfit  [semantic]
// =============================================================================

TEST( Polynomial, PolyfitLinearRoundTrip )
{
    // fit y = 2x + 3 to degree-1 polynomial
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY( 5 );
    for( size_t i = 0; i < 5; ++i )
    {
        tY( i ) = 2.0 * tX( i ) + 3.0;
    }

    belfem::Vector< belfem::real > tCoeffs;
    belfem::polyfit( tX, tY, 1, tCoeffs );

    // coefficients should be approximately {2.0, 3.0}
    EXPECT_NEAR( tCoeffs( 0 ), 2.0, tTol );
    EXPECT_NEAR( tCoeffs( 1 ), 3.0, tTol );
}

TEST( Polynomial, PolyfitQuadraticRoundTrip )
{
    // fit y = x² to degree-2
    belfem::Vector< belfem::real > tX = { -2.0, -1.0, 0.0, 1.0, 2.0 };
    belfem::Vector< belfem::real > tY( 5 );
    for( size_t i = 0; i < 5; ++i )
    {
        tY( i ) = tX( i ) * tX( i );
    }

    belfem::Vector< belfem::real > tCoeffs;
    belfem::polyfit( tX, tY, 2, tCoeffs );

    // evaluate fitted polynomial at training points
    for( size_t i = 0; i < 5; ++i )
    {
        EXPECT_NEAR( belfem::polyval( tCoeffs, tX( i ) ), tY( i ), tTol );
    }
}

TEST( Polynomial, PolyfitReproducesData )
{
    // generic round-trip: polyval(polyfit(x, y, n), x) ≈ y
    belfem::Vector< belfem::real > tX = { 0.0, 0.5, 1.0, 1.5, 2.0 };
    belfem::Vector< belfem::real > tY = { 1.0, 1.25, 2.0, 3.25, 5.0 };

    belfem::Vector< belfem::real > tCoeffs;
    belfem::polyfit( tX, tY, 2, tCoeffs );

    for( size_t i = 0; i < tX.length(); ++i )
    {
        EXPECT_NEAR( belfem::polyval( tCoeffs, tX( i ) ), tY( i ), tTol );
    }
}

TEST( Polynomial, PolyfitDegreeEightRecoversCoefficients )
{
    // Both backends scale x by ( n - 1 ) / ( x_n - x_0 ), which maps [ 2, 5 ]
    // to unit spacing on [ 39.3, 98.3 ], and solve the Vandermonde system by
    // QR least squares: Blaze through gels, Armadillo inside arma::polyfit.
    // The coefficient error is dominated by the low powers. The worst
    // coefficient measured on this fixture is 1.2e-5 ( Armadillo, MKL ),
    // 5.9e-6 ( Blaze, MKL ), 6.8e-6 ( reference LAPACK ) and 3.6e-6
    // ( OpenBLAS ); the bound sits a factor of eight above the worst of them.
    // A normal-equations solve fails every coefficient of the same fixture,
    // the worst by 8e+2, and the off-grid value too.
    const belfem::real tCoeffTol = 1e-4;

    const belfem::Vector< belfem::real > tExact =
            { 1.0, -2.5, 3.25, -0.75, 1.5, -4.0, 2.0, -0.5, 3.0 };

    // coefficients are stored highest power first, so Horner runs forward
    auto tHorner = []( const belfem::Vector< belfem::real > & aC,
                       const belfem::real aX )
    {
        belfem::real tValue = aC( 0 );
        for( size_t i = 1; i < aC.length(); ++i )
        {
            tValue = tValue * aX + aC( i );
        }
        return tValue;
    };

    // 60 samples on [ 2, 5 ]
    const size_t tNumSamples = 60;
    belfem::Vector< belfem::real > tX( tNumSamples );
    belfem::Vector< belfem::real > tY( tNumSamples );

    for( size_t k = 0; k < tNumSamples; ++k )
    {
        tX( k ) = 2.0 + 3.0 * ( belfem::real ) k / 59.0;

        // the samples come from Horner, not from polyval, so the fixture does
        // not depend on the function the neighboring cases pin
        tY( k ) = tHorner( tExact, tX( k ) );
    }

    belfem::Vector< belfem::real > tCoeffs;
    belfem::polyfit( tX, tY, 8, tCoeffs );

    ASSERT_EQ( tCoeffs.length(), tExact.length() );

    for( size_t i = 0; i < tExact.length(); ++i )
    {
        EXPECT_LT( std::abs( tCoeffs( i ) - tExact( i ) ) / std::abs( tExact( i ) ),
                   tCoeffTol ) << "coefficient of power " << ( 8 - i );
    }

    const belfem::real tProbe = 3.37;
    const belfem::real tRef = tHorner( tExact, tProbe );

    // the fitted value is far better conditioned than the coefficients that
    // produce it; measured at 4e-16 or better
    EXPECT_LT( std::abs( tHorner( tCoeffs, tProbe ) - tRef ) / std::abs( tRef ), 1e-12 );
}

// =============================================================================
// §6.5  Polyfit  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( PolynomialDebug, PolyfitLengthMismatchThrows )
{
    belfem::Vector< belfem::real > tX( 5, 1.0 );
    belfem::Vector< belfem::real > tY( 4, 1.0 );   // different length
    belfem::Vector< belfem::real > tCoeffs;

    EXPECT_THROW( belfem::polyfit( tX, tY, 2, tCoeffs ), std::runtime_error );
}

TEST( PolynomialDebug, PolyfitInsufficientSamplesThrows )
{
    // degree 3 needs at least 4 samples, but we only have 3
    belfem::Vector< belfem::real > tX = { 0.0, 1.0, 2.0 };
    belfem::Vector< belfem::real > tY = { 0.0, 1.0, 4.0 };
    belfem::Vector< belfem::real > tCoeffs;

    EXPECT_THROW( belfem::polyfit( tX, tY, 3, tCoeffs ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §6.6  R² (coefficient of determination)  [semantic]
// =============================================================================

TEST( R2, R2PerfectFit )
{
    belfem::Vector< belfem::real > tY = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    EXPECT_NEAR( belfem::r2( tY, tY ), 1.0, tEps );
}

TEST( R2, R2KnownValue )
{
    // exact: {1, 2, 3}, approximated: {1.1, 2.0, 2.9}
    belfem::Vector< belfem::real > tExact  = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tApprox = { 1.1, 2.0, 2.9 };

    belfem::real tResult = belfem::r2( tApprox, tExact );

    // manually: SSres = 0.01 + 0 + 0.01 = 0.02
    //           mean = 2.0
    //           SStot = 1 + 0 + 1 = 2
    //           R² = 1 - 0.02/2 = 0.99
    EXPECT_NEAR( tResult, 0.99, tTol );
}

TEST( R2, R2ConstantExact )
{
    // all exact values constant, approximation matches → R² = 1.0
    // (SStot < BELFEM_EPSILON triggers the early-return branch)
    belfem::Vector< belfem::real > tExact  = { 5.0, 5.0, 5.0 };
    belfem::Vector< belfem::real > tApprox = { 5.0, 5.0, 5.0 };

    EXPECT_NEAR( belfem::r2( tApprox, tExact ), 1.0, tEps );
}

TEST( R2, R2MatrixVersion )
{
    // BUG-L2 regression test: use non-square matrix (3 rows × 5 cols)
    // to catch the inner loop iterating j<n (cols) instead of j<m (rows)
    belfem::Matrix< belfem::real > tExact( 3, 5 );
    belfem::Matrix< belfem::real > tApprox( 3, 5 );

    // fill with known values
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 5; ++j )
        {
            belfem::real tVal = static_cast< belfem::real >( i * 5 + j + 1 );
            tExact( i, j )  = tVal;
            tApprox( i, j ) = tVal;   // perfect approximation
        }
    }

    // perfect fit → R² should be 1.0
    EXPECT_NEAR( belfem::r2( tApprox, tExact ), 1.0, tEps );

    // now introduce known error in the approximation
    tApprox( 0, 0 ) += 0.1;
    belfem::real tResult = belfem::r2( tApprox, tExact );
    EXPECT_GT( tResult, 0.99 );
    EXPECT_LT( tResult, 1.0 );
}

// =============================================================================
// §6.7  Eigen  [semantic]
// =============================================================================

TEST( Eigen, EigenDiagonalMatrix )
{
    // diag(1, 2, 3) → eigenvalues {1, 2, 3} in some order
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0, 0.0 },
                                           { 0.0, 2.0, 0.0 },
                                           { 0.0, 0.0, 3.0 } };
    belfem::Vector< belfem::real > tValues;
    belfem::eigen( tA, tValues );

    EXPECT_EQ( tValues.length(), 3u );

    // eigenvalues may be in any order — check sorted
    belfem::Vector< belfem::real > tSorted( tValues );
    // simple bubble for 3 elements
    if( tSorted( 0 ) > tSorted( 1 ) ) std::swap( tSorted( 0 ), tSorted( 1 ) );
    if( tSorted( 1 ) > tSorted( 2 ) ) std::swap( tSorted( 1 ), tSorted( 2 ) );
    if( tSorted( 0 ) > tSorted( 1 ) ) std::swap( tSorted( 0 ), tSorted( 1 ) );

    EXPECT_NEAR( tSorted( 0 ), 1.0, tTol );
    EXPECT_NEAR( tSorted( 1 ), 2.0, tTol );
    EXPECT_NEAR( tSorted( 2 ), 3.0, tTol );
}

TEST( Eigen, EigenSymmetricMatrix )
{
    // symmetric → all real eigenvalues, no NaN
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 },
                                           { 1.0, 3.0 } };
    belfem::Vector< belfem::real > tValues;
    belfem::eigen( tA, tValues );

    EXPECT_EQ( tValues.length(), 2u );
    // both must be real (no NaN)
    EXPECT_FALSE( std::isnan( tValues( 0 ) ) );
    EXPECT_FALSE( std::isnan( tValues( 1 ) ) );

    // eigenvalues of {{2,1},{1,3}} are (5±√5)/2 ≈ 1.382, 3.618
    belfem::real tLambda1 = ( 5.0 - std::sqrt( 5.0 ) ) / 2.0;
    belfem::real tLambda2 = ( 5.0 + std::sqrt( 5.0 ) ) / 2.0;

    belfem::real tMin = std::min( tValues( 0 ), tValues( 1 ) );
    belfem::real tMax = std::max( tValues( 0 ), tValues( 1 ) );
    EXPECT_NEAR( tMin, tLambda1, tTol );
    EXPECT_NEAR( tMax, tLambda2, tTol );
}

TEST( Eigen, EigenComplexEigenvaluesAbortByDefault )
{
    // rotation matrix: eigenvalues are complex (cos(θ) ± i·sin(θ)), which a real
    // Vector cannot hold. Since 2026-08-18 that is an error by default rather than
    // a silent NaN.
    belfem::real tTheta = 1.0;   // ~57 degrees
    belfem::Matrix< belfem::real > tA = {
        { std::cos( tTheta ), -std::sin( tTheta ) },
        { std::sin( tTheta ),  std::cos( tTheta ) } };

    belfem::Vector< belfem::real > tValues;

    EXPECT_THROW( belfem::eigen( tA, tValues ), std::runtime_error );
}

TEST( Eigen, EigenComplexEigenvaluesReturnNaNWhenNotAborting )
{
    // the same matrix, with the abort switched off: the complex entries come back
    // as NaN and the return value counts them
    belfem::real tTheta = 1.0;
    belfem::Matrix< belfem::real > tA = {
        { std::cos( tTheta ), -std::sin( tTheta ) },
        { std::sin( tTheta ),  std::cos( tTheta ) } };

    belfem::Vector< belfem::real > tValues;
    belfem::int_t tNumComplex = belfem::eigen( tA, tValues, false );

    EXPECT_EQ( tValues.length(), 2u );
    EXPECT_EQ( tNumComplex, 2 );
    EXPECT_TRUE( std::isnan( tValues( 0 ) ) );
    EXPECT_TRUE( std::isnan( tValues( 1 ) ) );
}

TEST( Eigen, EigenRealEigenvaluesReturnZeroCount )
{
    // a matrix with real eigenvalues reports none complex
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 },
                                           { 1.0, 3.0 } };
    belfem::Vector< belfem::real > tValues;

    EXPECT_EQ( belfem::eigen( tA, tValues, false ), 0 );
}

// -----------------------------------------------------------------------------
// eigen_sym: symmetric matrices have real eigenvalues by construction
// -----------------------------------------------------------------------------

TEST( EigenSym, EigenSymReturnsAscendingRealValues )
{
    // eigenvalues of {{2,1},{1,3}} are (5±√5)/2 ≈ 1.381966, 3.618034
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 },
                                           { 1.0, 3.0 } };
    belfem::Vector< belfem::real > tValues;
    belfem::eigen_sym( tA, tValues );

    belfem::real tTol     = 1e-9;
    belfem::real tLambda1 = 0.5 * ( 5.0 - std::sqrt( 5.0 ) );
    belfem::real tLambda2 = 0.5 * ( 5.0 + std::sqrt( 5.0 ) );

    ASSERT_EQ( tValues.length(), 2u );
    EXPECT_NEAR( tValues( 0 ), tLambda1, tTol );   // ascending order
    EXPECT_NEAR( tValues( 1 ), tLambda2, tTol );
}

TEST( EigenSym, EigenSymDiagonal )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0, 0.0 },
                                           { 0.0, 3.0, 0.0 },
                                           { 0.0, 0.0, 2.0 } };
    belfem::Vector< belfem::real > tValues;
    belfem::eigen_sym( tA, tValues );

    belfem::real tTol = 1e-9;
    ASSERT_EQ( tValues.length(), 3u );
    EXPECT_NEAR( tValues( 0 ), 1.0, tTol );
    EXPECT_NEAR( tValues( 1 ), 2.0, tTol );
    EXPECT_NEAR( tValues( 2 ), 3.0, tTol );
}

TEST( EigenSym, EigenSymAgreesWithEigenOnASymmetricMatrix )
{
    // the two entry points must not disagree where both are valid
    belfem::Matrix< belfem::real > tA = { { 4.0, 1.0 },
                                           { 1.0, 4.0 } };

    belfem::Vector< belfem::real > tSym;
    belfem::eigen_sym( tA, tSym );

    belfem::Vector< belfem::real > tGen;
    belfem::eigen( tA, tGen );

    belfem::real tTol = 1e-9;
    ASSERT_EQ( tSym.length(), tGen.length() );

    belfem::real tGenMin = std::min( tGen( 0 ), tGen( 1 ) );
    belfem::real tGenMax = std::max( tGen( 0 ), tGen( 1 ) );

    EXPECT_NEAR( tSym( 0 ), tGenMin, tTol );
    EXPECT_NEAR( tSym( 1 ), tGenMax, tTol );
}

// -----------------------------------------------------------------------------
// eigen_sym: the spectrum is known in advance, only the upper triangle is read,
// and the input is not written
//
// S = Q diag( 1, 2, 3, 7 ) Q^T with the rational orthogonal
//     Q = 1/2 * { { 1, 1, 1, 1 }, { 1, 1, -1, -1 },
//                 { 1, -1, 1, -1 }, { 1, -1, -1, 1 } } .
// Every entry of S is a multiple of 1/4, so the construction below is exact
// in floating point.
// -----------------------------------------------------------------------------

namespace
{
    void
    build_spectral_matrix( belfem::Matrix< belfem::real > & aS )
    {
        const belfem::Matrix< belfem::real > tQ = { {  0.5,  0.5,  0.5,  0.5 },
                                                    {  0.5,  0.5, -0.5, -0.5 },
                                                    {  0.5, -0.5,  0.5, -0.5 },
                                                    {  0.5, -0.5, -0.5,  0.5 } };

        const belfem::Vector< belfem::real > tLambda = { 1.0, 2.0, 3.0, 7.0 };

        aS.set_size( 4, 4, 0.0 );

        for ( belfem::uint j = 0; j < 4; ++j )
            for ( belfem::uint i = 0; i < 4; ++i )
            {
                belfem::real tSum = 0.0;
                for ( belfem::uint k = 0; k < 4; ++k )
                {
                    tSum += tQ( i, k ) * tLambda( k ) * tQ( j, k );
                }
                aS( i, j ) = tSum;
            }
    }

    // the six strictly-lower entries, none of them the mirror of its upper twin
    void
    scribble_lower_triangle( belfem::Matrix< belfem::real > & aA )
    {
        aA( 1, 0 ) = -7.5;
        aA( 2, 0 ) =  9.25;
        aA( 2, 1 ) = 11.0;
        aA( 3, 0 ) = -0.125;
        aA( 3, 1 ) =  3.75;
        aA( 3, 2 ) = -6.5;
    }
}

TEST( EigenSym, EigenSymRecoversAConstructedSpectrum )
{
    belfem::Matrix< belfem::real > tS;
    build_spectral_matrix( tS );

    belfem::Vector< belfem::real > tValues;
    belfem::eigen_sym( tS, tValues );

    ASSERT_EQ( tValues.length(), 4u );
    EXPECT_NEAR( tValues( 0 ), 1.0, tTol );
    EXPECT_NEAR( tValues( 1 ), 2.0, tTol );
    EXPECT_NEAR( tValues( 2 ), 3.0, tTol );
    EXPECT_NEAR( tValues( 3 ), 7.0, tTol );
}

TEST( EigenSym, EigenSymReadsOnlyTheUpperTriangle )
{
    // the whole lower triangle is garbage, the spectrum must not move.
    // Armadillo prints one
    //     warning: eig_sym(): given matrix is not symmetric
    // line for this call in any build -- its conformance checks are on
    // regardless of NDEBUG and sample two of the entries scribbled here. It
    // is a warning, not a failure, and Blaze has no such check.
    belfem::Matrix< belfem::real > tA;
    build_spectral_matrix( tA );
    scribble_lower_triangle( tA );

    belfem::Vector< belfem::real > tValues;
    belfem::eigen_sym( tA, tValues );

    ASSERT_EQ( tValues.length(), 4u );
    EXPECT_NEAR( tValues( 0 ), 1.0, tTol );
    EXPECT_NEAR( tValues( 1 ), 2.0, tTol );
    EXPECT_NEAR( tValues( 2 ), 3.0, tTol );
    EXPECT_NEAR( tValues( 3 ), 7.0, tTol );
}

TEST( EigenSym, EigenSymLeavesTheInputIntact )
{
    // nothing may write the caller's matrix: the Blaze wrapper works on a copy
    // because syev overwrites its argument, Armadillo copies inside eig_sym.
    // Entries are compared through the accessor with an exact equality --
    // data() is not the same layout under both backends, and a decomposition
    // that touched the matrix at all is a defect. The scribbled lower
    // triangle draws the same Armadillo warning as the previous test.
    belfem::Matrix< belfem::real > tA;
    build_spectral_matrix( tA );
    scribble_lower_triangle( tA );

    const belfem::Matrix< belfem::real > tA0( tA );

    belfem::Vector< belfem::real > tValues;
    belfem::eigen_sym( tA, tValues );

    for ( belfem::uint j = 0; j < 4; ++j )
        for ( belfem::uint i = 0; i < 4; ++i )
        {
            EXPECT_EQ( tA( i, j ), tA0( i, j ) );
        }
}

TEST( Eigen, EigenResultLength )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 0.0, 0.0 },
                                           { 0.0, 2.0, 0.0 },
                                           { 0.0, 0.0, 3.0 } };
    belfem::Vector< belfem::real > tValues;
    belfem::eigen( tA, tValues );

    EXPECT_EQ( tValues.length(), tA.n_cols() );
}

// =============================================================================
// §6.8  Eigen  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( EigenDebug, EigenNonSquareThrows )
{
    belfem::Matrix< belfem::real > tA( 3, 4, 1.0 );
    belfem::Vector< belfem::real > tValues;

    EXPECT_THROW( belfem::eigen( tA, tValues ), std::runtime_error );
}

#endif // NDEBUG
