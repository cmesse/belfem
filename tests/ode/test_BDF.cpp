/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for BDF (Backward Differentiation Formula) integrator.
 * See: tests_12_ode.md §1–§2
 *
 * Deferred: stiff systems, implicit BDF, BDF order > 6, event handling.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_ShiftRegister.hpp"
#include "cl_BDF.hpp"

namespace
{
    const belfem::real tEps = 1e-12;
}

// =============================================================================
// §1  BDF Coefficient Verification  [semantic]
// =============================================================================

TEST( BDFCoefficients, Order1 )
{
    // BDF-1 (Backward Euler): α = [1, -1]
    belfem::ShiftRegister< belfem::real > tH( 1 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tSolver.compute_coefficients();

    const belfem::Vector< belfem::real > & tC = tSolver.coefficients();
    ASSERT_EQ( tC.length(), 2u );
    EXPECT_NEAR( tC( 0 ),  1.0, tEps );
    EXPECT_NEAR( tC( 1 ), -1.0, tEps );
}

TEST( BDFCoefficients, Order2 )
{
    // BDF-2: α = [3/2, -2, 1/2]
    belfem::ShiftRegister< belfem::real > tH( 2 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 1.0 );
    tSolver.compute_coefficients();

    const belfem::Vector< belfem::real > & tC = tSolver.coefficients();
    ASSERT_EQ( tC.length(), 3u );
    EXPECT_NEAR( tC( 0 ),  1.5,  tEps );
    EXPECT_NEAR( tC( 1 ), -2.0,  tEps );
    EXPECT_NEAR( tC( 2 ),  0.5,  tEps );
}

TEST( BDFCoefficients, Order3 )
{
    // BDF-3: α = [11/6, -3, 3/2, -1/3]
    belfem::ShiftRegister< belfem::real > tH( 3 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 1.0 );
    tH.push( 1.0 );
    tSolver.compute_coefficients();

    const belfem::Vector< belfem::real > & tC = tSolver.coefficients();
    ASSERT_EQ( tC.length(), 4u );
    EXPECT_NEAR( tC( 0 ),  11.0 / 6.0, tEps );
    EXPECT_NEAR( tC( 1 ), -3.0,        tEps );
    EXPECT_NEAR( tC( 2 ),  3.0 / 2.0,  tEps );
    EXPECT_NEAR( tC( 3 ), -1.0 / 3.0,  tEps );
}

TEST( BDFCoefficients, Order4 )
{
    // BDF-4: α = [25/12, -4, 3, -4/3, 1/4]
    belfem::ShiftRegister< belfem::real > tH( 4 );
    belfem::ode::BDF tSolver( tH );

    for( int k = 0; k < 4; ++k )
    {
        tH.push( 1.0 );
    }
    tSolver.compute_coefficients();

    const belfem::Vector< belfem::real > & tC = tSolver.coefficients();
    ASSERT_EQ( tC.length(), 5u );
    EXPECT_NEAR( tC( 0 ),  25.0 / 12.0, tEps );
    EXPECT_NEAR( tC( 1 ), -4.0,         tEps );
    EXPECT_NEAR( tC( 2 ),  3.0,         tEps );
    EXPECT_NEAR( tC( 3 ), -4.0 / 3.0,   tEps );
    EXPECT_NEAR( tC( 4 ),  1.0 / 4.0,   tEps );
}

TEST( BDFCoefficients, VariableStepConsistency )
{
    // For ANY step sizes, sum of coefficients must be zero
    // (constant function has zero derivative)
    belfem::ShiftRegister< belfem::real > tH( 4 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 0.8 );
    tH.push( 1.2 );
    tH.push( 0.5 );
    tSolver.compute_coefficients();

    const belfem::Vector< belfem::real > & tC = tSolver.coefficients();

    belfem::real tSum = 0.0;
    for( belfem::index_t k = 0; k < tC.length(); ++k )
    {
        tSum += tC( k );
    }
    EXPECT_NEAR( tSum, 0.0, tEps );
}

TEST( BDFCoefficients, OrderRampsWithHistory )
{
    // (idea from ChatGPT) — coefficient vector length grows as history fills
    belfem::ShiftRegister< belfem::real > tH( 5 );
    belfem::ode::BDF tSolver( tH );

    for( belfem::uint k = 1; k <= 5; ++k )
    {
        tH.push( 1.0 );
        tSolver.compute_coefficients();
        EXPECT_EQ( tSolver.coefficients().length(), k + 1 );
    }
}

TEST( BDFCoefficients, NoUpdateFlagReusesCoefficients )
{
    // Verify that the false flag prevents coefficient recomputation.
    // Note: deval always divides by current mH(0), so we check the
    // coefficient vector itself rather than the deval result.
    belfem::ShiftRegister< belfem::real > tH( 3 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 1.0 );
    tH.push( 1.0 );
    tSolver.compute_coefficients();

    // record coefficients before mutation
    const belfem::Vector< belfem::real > & tCoeffs = tSolver.coefficients();
    belfem::Vector< belfem::real > tCoeffsBefore( tCoeffs.length() );
    for( belfem::index_t k = 0; k < tCoeffs.length(); ++k )
    {
        tCoeffsBefore( k ) = tCoeffs( k );
    }

    // mutate step history
    tH.push( 2.0 );
    tH.push( 0.5 );

    // deval with false — coefficients must remain unchanged
    belfem::ShiftRegister< belfem::real > tY( 4 );
    tY.push( 0.0 );
    tY.push( 1.0 );
    tY.push( 2.0 );
    tY.push( 3.0 );

    tSolver.deval( tY, false );
    for( belfem::index_t k = 0; k < tCoeffsBefore.length(); ++k )
    {
        EXPECT_NEAR( tCoeffs( k ), tCoeffsBefore( k ), tEps );
    }

    // deval with true — coefficients must change
    tSolver.deval( tY, true );
    bool tChanged = false;
    for( belfem::index_t k = 0; k < tCoeffsBefore.length(); ++k )
    {
        if( std::abs( tCoeffs( k ) - tCoeffsBefore( k ) ) > 1e-6 )
        {
            tChanged = true;
        }
    }
    EXPECT_TRUE( tChanged );

    // Also verify eval(..., false) preserves coefficients (Codex finding #1)
    belfem::Vector< belfem::real > tCoeffsAfterUpdate( tCoeffs.length() );
    for( belfem::index_t k = 0; k < tCoeffs.length(); ++k )
    {
        tCoeffsAfterUpdate( k ) = tCoeffs( k );
    }

    tH.push( 0.3 );  // mutate again

    belfem::ShiftRegister< belfem::real > tYe( 4 );
    tYe.push( 0.0 );
    tYe.push( 1.0 );
    tYe.push( 2.0 );
    tYe.push( 0.0 );

    tSolver.eval( tYe, 1.0, false );
    for( belfem::index_t k = 0; k < tCoeffsAfterUpdate.length(); ++k )
    {
        EXPECT_NEAR( tCoeffs( k ), tCoeffsAfterUpdate( k ), tEps );
    }
}

// =============================================================================
// §2  BDF Eval and Deval  [semantic]
// =============================================================================

TEST( BDFDeval, ConstantSignalReturnsZero )
{
    // (idea from ChatGPT) — derivative of constant is zero
    belfem::ShiftRegister< belfem::real > tH( 3 );
    belfem::ShiftRegister< belfem::real > tY( 4 );
    belfem::ode::BDF tSolver( tH );

    // push uniform steps and constant y values
    for( int k = 0; k < 3; ++k )
    {
        tH.push( 1.0 );
    }
    for( int k = 0; k < 4; ++k )
    {
        tY.push( 5.0 );  // constant value
    }

    belfem::real tDeriv = tSolver.deval( tY );
    EXPECT_NEAR( tDeriv, 0.0, tEps );
}

TEST( BDFDeval, ScalarLinear )
{
    // y(t) = t → dy/dt = 1
    // Use BDF-3 with uniform h=1: t = {3, 2, 1, 0}, y = {3, 2, 1, 0}
    belfem::ShiftRegister< belfem::real > tH( 3 );
    belfem::ShiftRegister< belfem::real > tY( 4 );
    belfem::ode::BDF tSolver( tH );

    for( int k = 0; k < 3; ++k )
    {
        tH.push( 1.0 );
    }

    // push in chronological order (oldest first → newest at index 0)
    tY.push( 0.0 );  // y(t=0)
    tY.push( 1.0 );  // y(t=1)
    tY.push( 2.0 );  // y(t=2)
    tY.push( 3.0 );  // y(t=3)

    belfem::real tDeriv = tSolver.deval( tY );
    EXPECT_NEAR( tDeriv, 1.0, tEps );
}

TEST( BDFDeval, ScalarQuadratic )
{
    // y(t) = t² → dy/dt = 2t
    // BDF-3 with h=1: times {3,2,1,0}, values {9,4,1,0}
    // deval at t=3 should approximate dy/dt(3) = 6
    belfem::ShiftRegister< belfem::real > tH( 3 );
    belfem::ShiftRegister< belfem::real > tY( 4 );
    belfem::ode::BDF tSolver( tH );

    for( int k = 0; k < 3; ++k )
    {
        tH.push( 1.0 );
    }

    tY.push( 0.0 );  // y(t=0) = 0
    tY.push( 1.0 );  // y(t=1) = 1
    tY.push( 4.0 );  // y(t=2) = 4
    tY.push( 9.0 );  // y(t=3) = 9

    belfem::real tDeriv = tSolver.deval( tY );

    // BDF-3 is exact for polynomials up to degree 3, so this should be exact
    EXPECT_NEAR( tDeriv, 6.0, tEps );
}

TEST( BDFEval, ScalarLinear )
{
    // y(t) = t, dy/dt = 1
    // Given history y(t=1)=1, y(t=0)=0 and f=1, eval should return y(t=2) = 2
    belfem::ShiftRegister< belfem::real > tH( 2 );
    belfem::ShiftRegister< belfem::real > tY( 3 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tY.push( 0.0 );  // y(t=0) — oldest
    tY.push( 1.0 );  // y(t=1)

    // push new step and placeholder (newest position)
    tH.push( 1.0 );
    tY.push( 0.0 );  // placeholder for y(t=2), will be filled by eval

    belfem::real tResult = tSolver.eval( tY, 1.0 );  // f = dy/dt = 1
    EXPECT_NEAR( tResult, 2.0, tEps );
}

TEST( BDFEval, ScalarQuadratic )
{
    // y(t) = t², dy/dt = 2t
    // History: y(2)=4, y(1)=1, y(0)=0, with h=1
    // At t=3: f = 2*3 = 6, expected y(3) = 9
    belfem::ShiftRegister< belfem::real > tH( 3 );
    belfem::ShiftRegister< belfem::real > tY( 4 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 1.0 );
    tY.push( 0.0 );  // y(t=0) — oldest
    tY.push( 1.0 );  // y(t=1)
    tY.push( 4.0 );  // y(t=2)

    // new step and placeholder
    tH.push( 1.0 );
    tY.push( 0.0 );  // placeholder for y(t=3)

    belfem::real tResult = tSolver.eval( tY, 6.0 );  // f = 2*3 = 6
    EXPECT_NEAR( tResult, 9.0, 1e-9 );
}

TEST( BDFDeval, VectorMatchesScalar )
{
    // Vector deval should produce same results as scalar on each component
    belfem::ShiftRegister< belfem::real > tH( 2 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 1.0 );

    // Scalar: y = t, push oldest first
    belfem::ShiftRegister< belfem::real > tYs( 3 );
    tYs.push( 0.0 );
    tYs.push( 1.0 );
    tYs.push( 2.0 );

    belfem::real tDerivScalar = tSolver.deval( tYs );

    // Vector: same values in a 1-component vector
    belfem::ShiftRegister< belfem::Vector< belfem::real > > tYv( 3 );
    belfem::Vector< belfem::real > tV0( 1 ), tV1( 1 ), tV2( 1 );
    tV0( 0 ) = 0.0;
    tV1( 0 ) = 1.0;
    tV2( 0 ) = 2.0;
    tYv.push( tV0 );
    tYv.push( tV1 );
    tYv.push( tV2 );

    const belfem::Vector< belfem::real > & tDerivVec = tSolver.deval( tYv );

    EXPECT_NEAR( tDerivVec( 0 ), tDerivScalar, tEps );
}

TEST( BDFEval, VectorMatchesScalar )
{
    // Vector eval should produce same results as scalar on each component
    belfem::ShiftRegister< belfem::real > tH( 2 );
    belfem::ode::BDF tSolver( tH );

    tH.push( 1.0 );
    tH.push( 1.0 );

    // Scalar: y = t, f = 1 — push oldest first, placeholder last
    belfem::ShiftRegister< belfem::real > tYs( 3 );
    tYs.push( 0.0 );   // y(t=0) — oldest
    tYs.push( 1.0 );   // y(t=1)
    tYs.push( 0.0 );   // placeholder for y(t=2)

    belfem::real tResultScalar = tSolver.eval( tYs, 1.0 );

    // Vector: same values
    belfem::ShiftRegister< belfem::Vector< belfem::real > > tYv( 3 );
    belfem::Vector< belfem::real > tV0( 1 ), tV1( 1 ), tV2( 1 );
    tV0( 0 ) = 0.0;   // y(t=0) — oldest
    tV1( 0 ) = 1.0;   // y(t=1)
    tV2( 0 ) = 0.0;   // placeholder
    tYv.push( tV0 );
    tYv.push( tV1 );
    tYv.push( tV2 );

    belfem::Vector< belfem::real > tF( 1 );
    tF( 0 ) = 1.0;

    const belfem::Vector< belfem::real > & tResultVec = tSolver.eval( tYv, tF );

    EXPECT_NEAR( tResultVec( 0 ), tResultScalar, tEps );
}

// =============================================================================
// §2.2  BDF Error Paths
// =============================================================================

TEST( BDFError, HighOrderThrows )
{
    // (idea from Gemini) — BELFEM_ERROR (always active): capacity >= 7
    belfem::ShiftRegister< belfem::real > tH( 7 );
    EXPECT_THROW( belfem::ode::BDF tSolver( tH ), std::runtime_error );
}

#ifndef NDEBUG
TEST( BDFError, EvalCapacityMismatchThrows )
{
    // (idea from Gemini) — BELFEM_ASSERT (debug only)
    belfem::ShiftRegister< belfem::real > tH( 2 );
    belfem::ode::BDF tSolver( tH );

    // tY should have capacity 3 (= tH.capacity() + 1), but we use 2
    belfem::ShiftRegister< belfem::real > tY( 2 );
    tH.push( 1.0 );
    tY.push( 1.0 );
    tY.push( 0.0 );

    EXPECT_THROW( tSolver.eval( tY, 1.0 ), std::runtime_error );
}

TEST( BDFError, DevalCapacityMismatchThrows )
{
    belfem::ShiftRegister< belfem::real > tH( 2 );
    belfem::ode::BDF tSolver( tH );

    belfem::ShiftRegister< belfem::real > tY( 2 );  // wrong: should be 3
    tH.push( 1.0 );
    tY.push( 1.0 );
    tY.push( 0.0 );

    EXPECT_THROW( tSolver.deval( tY ), std::runtime_error );
}

TEST( BDFError, EvalVectorCapacityMismatchThrows )
{
    // (idea from Gemini) — vector overload capacity check
    belfem::ShiftRegister< belfem::real > tH( 1 );
    belfem::ode::BDF tSolver( tH );

    // capacity should be tH.capacity()+1 = 2, but we use 3
    belfem::ShiftRegister< belfem::Vector< belfem::real > > tY( 3 );
    belfem::Vector< belfem::real > tF = { 0.0 };

    EXPECT_THROW( tSolver.eval( tY, tF ), std::runtime_error );
}
#endif
