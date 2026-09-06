/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_ShiftRegister.hpp"
#include "fn_FEM_anderson_mixing.hpp"

using namespace belfem;

namespace
{
    constexpr index_t tN     = 4 ;
    constexpr uint    tDepth = 3 ;

    //! bundles the scratch the mixing step needs, presized per contract
    struct MixingScratch
    {
        Matrix< real > mDeltaR ;
        Vector< real > mRhs ;
        Vector< real > mWork ;
        Vector< real > mGamma ;
        Vector< real > mColNorm ;
        Vector< real > mXNew ;

        MixingScratch()
        {
            mDeltaR.set_size( tN, tDepth );
            mRhs.set_size( tN );
            mGamma.set_size( tDepth );
            mColNorm.set_size( tDepth );
            mXNew.set_size( tN );
        }
    };

    //! the affine test map G( x ) = a * x + c with uniform contraction a
    void
    eval_map( const real aA, const Vector< real > & aC,
              const Vector< real > & aX, Vector< real > & aR )
    {
        for ( index_t k = 0; k < tN; ++k )
        {
            aR( k ) = aA * aX( k ) + aC( k ) - aX( k );
        }
    }
}

//------------------------------------------------------------------------------

TEST( AndersonMixing, empty_history_is_plain_relaxed_step )
{
    MixingScratch tScratch ;
    ShiftRegister< Vector< real > > tXHist( tDepth );
    ShiftRegister< Vector< real > > tRHist( tDepth );

    Vector< real > tX( tN );
    Vector< real > tR( tN );
    for ( index_t k = 0; k < tN; ++k )
    {
        tX( k ) = 1.0 + k ;
        tR( k ) = 0.5 - 0.1 * k ;
    }

    const real tBeta = 0.7 ;
    uint tUsed = fem::anderson_mixing_step( tX, tR, tXHist, tRHist, tBeta,
        tScratch.mDeltaR, tScratch.mRhs, tScratch.mWork,
        tScratch.mGamma, tScratch.mColNorm, tScratch.mXNew );

    EXPECT_EQ( tUsed, 0u );
    for ( index_t k = 0; k < tN; ++k )
    {
        EXPECT_NEAR( tScratch.mXNew( k ), tX( k ) + tBeta * tR( k ), 1e-14 );
    }
}

//------------------------------------------------------------------------------

TEST( AndersonMixing, affine_map_converges_in_one_mixed_step )
{
    // uniform contraction: one committed pair makes the secant exact,
    // so the first mixed step must land on the fixed point.
    // NOTE when retuning ( tA, tBeta ): the mixing coefficient is
    // gamma = ( 1 + beta ( a - 1 ) ) / ( beta ( a - 1 ) ) — here -19; a
    // choice with |gamma| > gAndersonGammaMax would trip the guard and
    // fail this test for the wrong reason
    const real tA = 0.9 ;
    Vector< real > tC( tN );
    for ( index_t k = 0; k < tN; ++k )
    {
        tC( k ) = 0.1 * ( k + 1 );
    }
    // fixed point x* = c / ( 1 - a )
    Vector< real > tXStar( tN );
    for ( index_t k = 0; k < tN; ++k )
    {
        tXStar( k ) = tC( k ) / ( 1.0 - tA );
    }

    MixingScratch tScratch ;
    ShiftRegister< Vector< real > > tXHist( tDepth );
    ShiftRegister< Vector< real > > tRHist( tDepth );

    Vector< real > tX( tN, 0.0 );
    Vector< real > tR( tN );
    const real tBeta = 0.5 ;

    // iteration 0: empty history -> plain step, then commit
    eval_map( tA, tC, tX, tR );
    uint tUsed = fem::anderson_mixing_step( tX, tR, tXHist, tRHist, tBeta,
        tScratch.mDeltaR, tScratch.mRhs, tScratch.mWork,
        tScratch.mGamma, tScratch.mColNorm, tScratch.mXNew );
    EXPECT_EQ( tUsed, 0u );
    tXHist.push( tX );
    tRHist.push( tR );
    tX = tScratch.mXNew ;

    // iteration 1: one column -> exact for the affine uniform map
    eval_map( tA, tC, tX, tR );
    tUsed = fem::anderson_mixing_step( tX, tR, tXHist, tRHist, tBeta,
        tScratch.mDeltaR, tScratch.mRhs, tScratch.mWork,
        tScratch.mGamma, tScratch.mColNorm, tScratch.mXNew );
    EXPECT_EQ( tUsed, 1u );

    for ( index_t k = 0; k < tN; ++k )
    {
        EXPECT_NEAR( tScratch.mXNew( k ), tXStar( k ), 1e-10 );
    }
}

//------------------------------------------------------------------------------

TEST( AndersonMixing, degenerate_column_shrinks_window )
{
    // two identical committed pairs: the oldest difference column vanishes,
    // the solver must drop it and succeed with the newest column only
    MixingScratch tScratch ;
    ShiftRegister< Vector< real > > tXHist( tDepth );
    ShiftRegister< Vector< real > > tRHist( tDepth );

    Vector< real > tXOld( tN ), tROld( tN );
    for ( index_t k = 0; k < tN; ++k )
    {
        tXOld( k ) = 1.0 ;
        tROld( k ) = 0.25 * ( k + 1 );
    }
    tXHist.push( tXOld );   tRHist.push( tROld );
    tXHist.push( tXOld );   tRHist.push( tROld );   // duplicate

    Vector< real > tX( tN ), tR( tN );
    for ( index_t k = 0; k < tN; ++k )
    {
        tX( k ) = 2.0 ;
        tR( k ) = 0.125 * ( k + 1 );
    }

    uint tUsed = fem::anderson_mixing_step( tX, tR, tXHist, tRHist, 1.0,
        tScratch.mDeltaR, tScratch.mRhs, tScratch.mWork,
        tScratch.mGamma, tScratch.mColNorm, tScratch.mXNew );

    EXPECT_EQ( tUsed, 1u );
}

//------------------------------------------------------------------------------

TEST( AndersonMixing, window_is_clamped_to_system_size )
{
    // more history columns than dofs: the window must clamp to n ( the
    // release-safe guard — an underdetermined solve would overrun the rhs )
    const index_t tSmallN   = 2 ;
    const uint    tBigDepth = 3 ;

    ShiftRegister< Vector< real > > tXHist( tBigDepth );
    ShiftRegister< Vector< real > > tRHist( tBigDepth );

    Matrix< real > tDeltaR( tSmallN, tBigDepth );
    Vector< real > tRhs( tSmallN );
    Vector< real > tWork ;
    Vector< real > tGamma( tBigDepth );
    Vector< real > tColNorm( tBigDepth );
    Vector< real > tXNew( tSmallN );

    // three distinct commits in R^2 — at most two columns can be used
    Vector< real > tX( tSmallN ), tR( tSmallN );
    tX( 0 ) = 0.0 ; tX( 1 ) = 0.0 ; tR( 0 ) = 4.0 ; tR( 1 ) = 0.0 ;
    tXHist.push( tX ); tRHist.push( tR );
    tX( 0 ) = 1.0 ; tX( 1 ) = 0.0 ; tR( 0 ) = 0.0 ; tR( 1 ) = 2.0 ;
    tXHist.push( tX ); tRHist.push( tR );
    tX( 0 ) = 0.0 ; tX( 1 ) = 1.0 ; tR( 0 ) = 1.0 ; tR( 1 ) = 1.0 ;
    tXHist.push( tX ); tRHist.push( tR );

    tX( 0 ) = 1.0 ; tX( 1 ) = 1.0 ; tR( 0 ) = 2.0 ; tR( 1 ) = 1.0 ;

    uint tUsed = fem::anderson_mixing_step( tX, tR, tXHist, tRHist, 1.0,
        tDeltaR, tRhs, tWork, tGamma, tColNorm, tXNew );

    EXPECT_EQ( tUsed, ( uint ) tSmallN );
}

//------------------------------------------------------------------------------

TEST( AndersonMixing, unusable_history_falls_back_to_plain_step )
{
    // current residual equals the newest committed one: every attempt sees a
    // vanishing newest column, so mixing must fall back to the plain step
    MixingScratch tScratch ;
    ShiftRegister< Vector< real > > tXHist( tDepth );
    ShiftRegister< Vector< real > > tRHist( tDepth );

    Vector< real > tX( tN ), tR( tN );
    for ( index_t k = 0; k < tN; ++k )
    {
        tX( k ) = 3.0 ;
        tR( k ) = 1.0 + k ;
    }
    tXHist.push( tX );
    tRHist.push( tR );

    const real tBeta = 0.25 ;
    uint tUsed = fem::anderson_mixing_step( tX, tR, tXHist, tRHist, tBeta,
        tScratch.mDeltaR, tScratch.mRhs, tScratch.mWork,
        tScratch.mGamma, tScratch.mColNorm, tScratch.mXNew );

    EXPECT_EQ( tUsed, 0u );
    for ( index_t k = 0; k < tN; ++k )
    {
        EXPECT_NEAR( tScratch.mXNew( k ), tX( k ) + tBeta * tR( k ), 1e-14 );
    }
}
