/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * Unit tests for IWG_Timestep's scheme selection and the BDF startup ramp.
 *
 * The defect this guards is the one the BDF Jacobian work was opened on:
 * compute_bdf_coefficients() was never called, so mAlpha/mBeta stayed NaN and
 * BDF2-5 assembled an all-NaN system.
 *
 * §1-§3 use the public interface only. §2 is not cosmetic: at step zero mH is
 * all zeros, so a ramp that failed to hold the first step at order 1 would
 * walk straight into compute_bdf_coefficients()'s own "Invalid BDF history
 * step size" assert.
 *
 * §4 covers the variable-step coefficient formulas themselves, and needs the
 * BdfCoefficientProbe seam declared in cl_IWG_Timestep.hpp: mAlpha, mBeta, mH
 * and mStepCount are private. The production writers of the step-size history
 * are shift_fields / reset_fields / restore_savepoint ( all needing an
 * initialized DofManager ) plus the validated memdump seam
 * restore_history_state ( 2026-08-15, needs none — tested in §5 below ). Nothing in the tree exercises BDF2-5 end to end —
 * it is opt-in, and the deck that would have done it no longer exists
 * (todo/bdf_nonlinear_mass_verification.md, V1) — so this is the only place
 * those formulas are checked at all.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <algorithm>
#include <limits>

#include "typedefs.hpp"
#include "en_SolverEnums.hpp"
#include "en_IWGs.hpp"
#include "cl_IWG_Timestep.hpp"

using namespace belfem;
using namespace belfem::fem;

namespace
{
    // IWG_Timestep is concrete and its constructor allocates nothing beyond
    // its own coefficient buffers, so it stands up without a mesh, a kernel
    // or a dof manager. IwgType::UNDEFINED keeps the physics out of it —
    // nothing below assembles anything.
    IWG_Timestep *
    make_iwg()
    {
        return new IWG_Timestep(
            IwgType::UNDEFINED,
            ModelDimensionality::ThreeD );
    }

    const EulerMethod tBdfMethods[ 4 ] =
    {
        EulerMethod::BackwardDifference2,
        EulerMethod::BackwardDifference3,
        EulerMethod::BackwardDifference4,
        EulerMethod::BackwardDifference5
    };
}

// -----------------------------------------------------------------------------
// the coefficient seam
// -----------------------------------------------------------------------------

namespace belfem
{
    namespace fem
    {
        /**
         * Poses a step-size history on an IWG_Timestep and reads back the
         * coefficients it produces. Declared a friend in cl_IWG_Timestep.hpp;
         * this is its only definition and nothing in src/ uses it.
         */
        class BdfCoefficientProbe
        {
            IWG_Timestep & mIWG ;

        public:

            BdfCoefficientProbe( IWG_Timestep & aIWG ) : mIWG( aIWG ) {}

            /**
             * Pose the state a BDF-p step would see well past the startup
             * ramp: aStepCount steps already taken, and the previous step
             * sizes h1..h4 in aH. aDeltaTime is the UPCOMING step.
             */
            void
            pose(
                const real aDeltaTime,
                const real aH[ 4 ],
                const uint aStepCount = 16 )
            {
                mIWG.delta_time() = aDeltaTime ;
                mIWG.mStepCount   = aStepCount ;

                // the coefficient formulas read mH( 0 ) .. mH( 3 ); mH( 4 )
                // is only ever a shift target and stays at its ctor zero
                for( uint k = 0; k < 4; ++k )
                {
                    mIWG.mH( k ) = aH[ k ] ;
                }
            }

            real alpha() const                { return mIWG.mAlpha ; }
            real beta( const uint aK ) const  { return mIWG.mBeta( aK ) ; }

            real h( const uint aK ) const     { return mIWG.mH( aK ) ; }
            uint step_count() const           { return mIWG.mStepCount ; }
            real delta_time() const           { return mIWG.mDeltaTime ; }

            /**
             * The five-line rotation of shift_fields, WITHOUT the field
             * machinery ( which needs a DofManager ). The warm-start
             * round-trip below must exercise exactly this rotation: a
             * restore-then-compute test with no shift would go green while
             * the restart still injected the fresh-object mDeltaTime = 1.0
             * into the history ( the P0 of the 2026-08-15 review ).
             */
            void
            shift()
            {
                mIWG.mHDropped = mIWG.mH( 3 );
                mIWG.mH( 3 ) = mIWG.mH( 2 );
                mIWG.mH( 2 ) = mIWG.mH( 1 );
                mIWG.mH( 1 ) = mIWG.mH( 0 );
                mIWG.mH( 0 ) = mIWG.mDeltaTime ;
                ++mIWG.mStepCount ;
                mIWG.mCoeffsDirty = true ;
            }
        };
    }
}

namespace
{
    /**
     * The quantity every BDF scheme exists to produce: the approximation of
     * dq/dt at t_n, assembled the way the solver assembles it —
     *
     *     ( alpha*q( t_n ) - qhist ) / h,
     *     qhist = beta0*q0 - beta1*q1 + beta2*q2 - beta3*q3 + beta4*q4
     *
     * with q0 the previous step. The alternating signs live in collect_qhist,
     * not in mBeta, so they are applied here too.
     *
     * Sample abscissae: t_n, then t_n - h, t_n - h - h1, and so on, since
     * mDeltaTime is the upcoming step and mH(0) the one before it.
     */
    template< typename F >
    belfem::real
    bdf_derivative(
        const belfem::fem::BdfCoefficientProbe & aProbe,
        const uint    aOrder,
        const real    aTn,
        const real    aDeltaTime,
        const real    aH[ 4 ],
        F             aFunction )
    {
        real tSum = aProbe.alpha() * aFunction( aTn );

        real tT = aTn - aDeltaTime ;

        for( uint i = 0; i < aOrder; ++i )
        {
            const real tSign = ( i % 2 == 0 ) ? -1.0 : 1.0 ;

            tSum += tSign * aProbe.beta( i ) * aFunction( tT );

            // aH holds h1..h4, so the last term has no gap left to walk
            if( i + 1 < aOrder )
            {
                tT -= aH[ i ];
            }
        }

        return tSum / aDeltaTime ;
    }
}

// =============================================================================
// §1  scheme selection
// =============================================================================

TEST( BdfTimestep, DefaultMethodIsBdf1 )
{
    // the constructor configures BDF1 itself, so the dispatch pointer is
    // never null even if a driver forgets to choose a scheme
    IWG_Timestep * tIWG = make_iwg();

    EXPECT_EQ( tIWG->method(), EulerMethod::BackwardDifference1 );
    EXPECT_EQ( tIWG->order_active(), 1u );

    delete tIWG ;
}

TEST( BdfTimestep, MethodRoundTrips )
{
    IWG_Timestep * tIWG = make_iwg();

    tIWG->set_timestepping_method( EulerMethod::BackwardDifference1 );
    EXPECT_EQ( tIWG->method(), EulerMethod::BackwardDifference1 );

    for( uint k = 0; k < 4; ++k )
    {
        tIWG->set_timestepping_method( tBdfMethods[ k ] );
        EXPECT_EQ( tIWG->method(), tBdfMethods[ k ] );
    }

    delete tIWG ;
}

// =============================================================================
// §2  startup ramp
// =============================================================================

TEST( BdfTimestep, RampHoldsFirstStepAtOrderOne )
{
    // Before any step has been shifted there is no history: mH is all zeros
    // and only the current state exists. Every BDF-p run must therefore
    // execute its first step at order 1, whatever p is. This is the
    // "never read an unpopulated mH slot" guarantee — without it,
    // compute_bdf_coefficients() reaches its own history-step-size assert.
    IWG_Timestep * tIWG = make_iwg();

    tIWG->delta_time() = 1.0e-3 ;

    for( uint k = 0; k < 4; ++k )
    {
        tIWG->set_timestepping_method( tBdfMethods[ k ] );
        tIWG->compute_bdf_coefficients();

        EXPECT_EQ( tIWG->order_active(), 1u )
            << "BDF" << ( k + 2 ) << " did not start at order 1";
    }

    delete tIWG ;
}

TEST( BdfTimestep, Bdf1IsNotRamped )
{
    // BDF1 has no ramp to run — it is order 1 from the first step and stays
    // there, so its behaviour must not depend on the step counter
    IWG_Timestep * tIWG = make_iwg();

    tIWG->delta_time() = 1.0e-3 ;
    tIWG->set_timestepping_method( EulerMethod::BackwardDifference1 );

    tIWG->compute_bdf_coefficients();
    EXPECT_EQ( tIWG->order_active(), 1u );

    // idempotent: recomputing without a shift changes nothing
    tIWG->compute_bdf_coefficients();
    EXPECT_EQ( tIWG->order_active(), 1u );

    delete tIWG ;
}

TEST( BdfTimestep, NonBdfSchemesSkipCoefficientWork )
{
    // compute_bdf_coefficients() gates on the method before touching the
    // ramp. That gate is what protects the temporary MassOnly/StiffnessOnly
    // switch in the eigenvalue path from inheriting stale BDF order state —
    // these schemes scale the mass matrix by one and must run exactly as
    // configured. Note delta_time is deliberately left at its default here:
    // these branches must not depend on it.
    IWG_Timestep * tIWG = make_iwg();

    tIWG->set_timestepping_method( EulerMethod::BackwardDifference5 );

    tIWG->set_timestepping_method( EulerMethod::MassOnly );
    tIWG->compute_bdf_coefficients();
    EXPECT_EQ( tIWG->method(), EulerMethod::MassOnly );
    EXPECT_EQ( tIWG->order_active(), 1u );

    tIWG->set_timestepping_method( EulerMethod::StiffnessOnly );
    tIWG->compute_bdf_coefficients();
    EXPECT_EQ( tIWG->method(), EulerMethod::StiffnessOnly );
    EXPECT_EQ( tIWG->order_active(), 1u );

    delete tIWG ;
}

// =============================================================================
// §3  schemes that are parseable but not usable  [debug]
// =============================================================================

// The rejection below is a BELFEM_ERROR, which is always active but does not
// always THROW: it throws under debug and calls error_abort() under release
// (assert.hpp). Catching it is therefore a debug-build test — under NDEBUG it
// would take the whole test binary down with it.
#ifndef NDEBUG

TEST( BdfTimestepDebug, CrankNicolsonAndGalerkinAreRejected )
{
    // Both are accepted by the input parser and then hard-error here, before
    // any dispatch pointer or tangent is built. That is deliberate — the
    // Newton tangent for them is inconsistent at variable dt — and it is why
    // the CN/Galerkin clause of the timestep debt row was struck rather than
    // fixed. This test is what keeps the rejection from being quietly lost.
    IWG_Timestep * tIWG = make_iwg();

    EXPECT_THROW( tIWG->set_timestepping_method( EulerMethod::CrankNicolson ),
                  std::runtime_error );
    EXPECT_THROW( tIWG->set_timestepping_method( EulerMethod::Galerkin ),
                  std::runtime_error );

    delete tIWG ;
}

#endif // NDEBUG

// =============================================================================
// §4  the variable-step coefficient formulas
// =============================================================================

TEST( BdfCoefficients, ConstantStepReducesToTextbook )
{
    // At uniform step size the variable-step Lagrange formulas must collapse
    // onto the classical BDF coefficients. Written out with the sign
    // convention of collect_qhist ( beta stored as magnitudes, alternating
    // signs applied there ), the schemes are
    //
    //   BDF2   3/2    q_n - 2 q_1 + 1/2 q_2
    //   BDF3  11/6    q_n - 3 q_1 + 3/2 q_2 - 1/3 q_3
    //   BDF4  25/12   q_n - 4 q_1 + 3   q_2 - 4/3 q_3 + 1/4 q_4
    //   BDF5 137/60   q_n - 5 q_1 + 5   q_2 - 10/3 q_3 + 5/4 q_4 - 1/5 q_5
    //
    // Independent reference, not a restatement of the implementation: these
    // are the published constants ( Hairer & Wanner 1996, II.4 ).
    const real tH = 2.5e-4 ;
    const real tHist[ 4 ] = { tH, tH, tH, tH };

    const real tAlpha[ 4 ] = { 3.0 / 2.0, 11.0 / 6.0, 25.0 / 12.0, 137.0 / 60.0 };

    const real tBeta[ 4 ][ 5 ] =
    {
        { 2.0, 1.0 / 2.0, 0.0,       0.0,       0.0       },
        { 3.0, 3.0 / 2.0, 1.0 / 3.0, 0.0,       0.0       },
        { 4.0, 3.0,       4.0 / 3.0, 1.0 / 4.0, 0.0       },
        { 5.0, 5.0,      10.0 / 3.0, 5.0 / 4.0, 1.0 / 5.0 }
    };

    for( uint k = 0; k < 4; ++k )
    {
        IWG_Timestep * tIWG = make_iwg();
        BdfCoefficientProbe tProbe( *tIWG );

        tIWG->set_timestepping_method( tBdfMethods[ k ] );
        tProbe.pose( tH, tHist );
        tIWG->compute_bdf_coefficients();

        const uint tOrder = k + 2 ;

        // past the ramp, the configured order is the one that runs
        EXPECT_EQ( tIWG->order_active(), tOrder );

        EXPECT_NEAR( tProbe.alpha(), tAlpha[ k ], 1e-13 )
            << "BDF" << tOrder << " alpha";

        for( uint i = 0; i < tOrder; ++i )
        {
            EXPECT_NEAR( tProbe.beta( i ), tBeta[ k ][ i ], 1e-13 )
                << "BDF" << tOrder << " beta" << i ;
        }

        delete tIWG ;
    }
}

TEST( BdfCoefficients, VariableStepIsExactOnPolynomials )
{
    // The property that actually defines a BDF-p scheme, and the one the
    // constant-step check above cannot reach: with NON-uniform step sizes the
    // formula must still differentiate every polynomial of degree <= p
    // exactly. That is the whole content of "Lagrange derivative formula" in
    // the implementation comments, and it is what a mistyped cumulative sum
    // would break while leaving the uniform case intact.
    //
    // Deliberately ragged history — growing, shrinking and a sharp cut, the
    // shapes an adaptive controller actually produces.
    //
    // WHY THE STEPS ARE O( 0.1 ) AND NOT A REALISTIC 1e-4. Exactness is scale
    // free ( exact is exact, down to roundoff ), but the FAILURE this test
    // must be able to see is not: a scheme of the wrong order misses by
    // O( h^p ), which at h = 3e-4 is around 1e-13 — indistinguishable from
    // the roundoff of a correct answer. Measured: with a 3e-4 history, BDF5
    // reproduces the derivative of t^6 to 2e-13, so the test would have
    // passed for a scheme it is supposed to reject. At these step sizes the
    // separation is fifteen orders of magnitude ( ~1e-15 exact against
    // 16-73 % wrong ), which is what the assertions below rely on.
    const real tDeltaTime = 0.30 ;
    const real tHist[ 4 ] = { 0.50, 0.125, 0.80, 0.20 };

    const real tTn = 0.7 ;

    for( uint k = 0; k < 4; ++k )
    {
        IWG_Timestep * tIWG = make_iwg();
        BdfCoefficientProbe tProbe( *tIWG );

        tIWG->set_timestepping_method( tBdfMethods[ k ] );
        tProbe.pose( tDeltaTime, tHist );
        tIWG->compute_bdf_coefficients();

        const uint tOrder = k + 2 ;
        ASSERT_EQ( tIWG->order_active(), tOrder );

        for( uint d = 0; d <= tOrder; ++d )
        {
            // q( t ) = t^d, so dq/dt at t_n is d * t_n^( d-1 )
            const real tExact = ( d == 0 )
                ? 0.0
                : ( real ) d * std::pow( tTn, ( real ) ( d - 1 ) );

            const real tNumeric = bdf_derivative(
                tProbe, tOrder, tTn, tDeltaTime, tHist,
                [ d ]( const real aT ) { return std::pow( aT, ( real ) d ); } );

            // relative: the sum cancels terms of order q/h, so the absolute
            // error scales with 1/h and a fixed tolerance would be meaningless
            EXPECT_NEAR( tNumeric, tExact, 1e-12 * std::max( 1.0, std::abs( tExact ) ) )
                << "BDF" << tOrder << " is not exact on t^" << d ;
        }

        // ...and it must NOT be exact one degree higher. Without this, the
        // loop above would pass just as happily for a scheme of too HIGH an
        // order — the property "exact up to degree p" is only a fingerprint
        // of BDF-p together with its own negation at p+1.
        {
            const uint d = tOrder + 1 ;

            const real tExact = ( real ) d * std::pow( tTn, ( real ) ( d - 1 ) );

            const real tNumeric = bdf_derivative(
                tProbe, tOrder, tTn, tDeltaTime, tHist,
                [ d ]( const real aT ) { return std::pow( aT, ( real ) d ); } );

            EXPECT_GT( std::abs( tNumeric - tExact ) / std::abs( tExact ), 1e-3 )
                << "BDF" << tOrder << " is exact on t^" << d
                << ", which no scheme of that order can be" ;
        }

        delete tIWG ;
    }
}

TEST( BdfCoefficients, RampUsesReducedOrderCoefficients )
{
    // Mid-ramp: a BDF5 run that has taken only two steps must produce BDF2's
    // coefficients, not BDF5's. This is the pairing collect_qhist depends on
    // — it truncates the history to mOrderActive, so an mAlpha computed at
    // the configured order would be contracted against a shorter history and
    // the Jacobian would silently stop being the tangent of its residual.
    const real tH = 1.0e-3 ;
    const real tHist[ 4 ] = { tH, tH, tH, tH };

    IWG_Timestep * tIWG = make_iwg();
    BdfCoefficientProbe tProbe( *tIWG );

    tIWG->set_timestepping_method( EulerMethod::BackwardDifference5 );

    // two steps taken: the ramp runs BDF2
    tProbe.pose( tH, tHist, 2 );
    tIWG->compute_bdf_coefficients();

    EXPECT_EQ( tIWG->order_active(), 2u );
    EXPECT_NEAR( tProbe.alpha(), 3.0 / 2.0, 1e-13 );
    EXPECT_NEAR( tProbe.beta( 0 ), 2.0, 1e-13 );
    EXPECT_NEAR( tProbe.beta( 1 ), 0.5, 1e-13 );

    // four steps taken: BDF4
    tProbe.pose( tH, tHist, 4 );
    tIWG->compute_bdf_coefficients();

    EXPECT_EQ( tIWG->order_active(), 4u );
    EXPECT_NEAR( tProbe.alpha(), 25.0 / 12.0, 1e-13 );

    // five steps taken: the configured order is finally reached
    tProbe.pose( tH, tHist, 5 );
    tIWG->compute_bdf_coefficients();

    EXPECT_EQ( tIWG->order_active(), 5u );
    EXPECT_NEAR( tProbe.alpha(), 137.0 / 60.0, 1e-13 );

    delete tIWG ;
}

//------------------------------------------------------------------------------
// §5  Memdump warm start: the ( mH, mStepCount, mDeltaTime ) round-trip
//------------------------------------------------------------------------------

// The resumed process must reach the SAME coefficients as the uninterrupted
// one. The uninterrupted twin runs one more step from a posed full-order
// state; the resumed twin saves that state, restores it into a FRESH object
// ( whose mDeltaTime would otherwise be the ctor 1.0 — the restart poison ),
// applies the shift rotation, and computes. Bitwise agreement required.
TEST( BdfWarmStart, RoundTripMatchesUninterruptedRun )
{
    const real tHist[ 4 ] = { 3.1e-4, 2.7e-4, 3.4e-4, 2.9e-4 };
    const real tLastDt = 2.55e-4 ;   // the just-completed step h_n
    const real tNextDt = 3.05e-4 ;   // the upcoming step h_n+1

    // --- uninterrupted twin -------------------------------------------------
    IWG_Timestep * tA = make_iwg();
    tA->set_timestepping_method( EulerMethod::BackwardDifference5 );
    BdfCoefficientProbe tProbeA( *tA );
    // state right after accepting step n: history h_n-1..h_n-4, the IWG's
    // delta_time still holds h_n ( the controller assigns h_n+1 only after
    // the next shift )
    tProbeA.pose( tLastDt, tHist, 16 );
    tProbeA.shift();                       // start step n+1
    tA->delta_time() = tNextDt ;
    tA->compute_bdf_coefficients();

    // --- resumed twin -------------------------------------------------------
    IWG_Timestep * tB = make_iwg();
    tB->set_timestepping_method( EulerMethod::BackwardDifference5 );
    BdfCoefficientProbe tProbeB( *tB );

    // save from a third object posed like the run at dump time
    IWG_Timestep * tDump = make_iwg();
    tDump->set_timestepping_method( EulerMethod::BackwardDifference5 );
    BdfCoefficientProbe tProbeD( *tDump );
    tProbeD.pose( tLastDt, tHist, 16 );

    Vector< real > tH ;
    uint tCount ;
    real tDt ;
    tDump->save_history_state( tH, tCount, tDt );

    EXPECT_EQ( tDt, tLastDt );             // the P0 quantity travels

    EXPECT_TRUE( tB->restore_history_state( tH, tCount, tDt ) );

    tProbeB.shift();                       // the first post-restore shift
    tB->delta_time() = tNextDt ;
    tB->compute_bdf_coefficients();

    EXPECT_EQ( tProbeB.alpha(), tProbeA.alpha() );
    for ( uint k = 0; k < 5; ++k )
    {
        EXPECT_EQ( tProbeB.beta( k ), tProbeA.beta( k ) );
        EXPECT_EQ( tProbeB.h( k ), tProbeA.h( k ) );
    }
    EXPECT_EQ( tB->order_active(), 5u );

    delete tA ;
    delete tB ;
    delete tDump ;
}

// The validation matrix: everything here must REJECT and leave the
// cold-start ramp untouched — except the spare-slot case, which must pass.
TEST( BdfWarmStart, RestoreValidation )
{
    const real tGood[ 5 ] = { 3.1e-4, 2.7e-4, 3.4e-4, 2.9e-4, 0.0 };

    Vector< real > tH( 5 );
    for ( uint k = 0; k < 5; ++k ) tH( k ) = tGood[ k ] ;

    IWG_Timestep * tIWG = make_iwg();
    tIWG->set_timestepping_method( EulerMethod::BackwardDifference5 );
    BdfCoefficientProbe tProbe( *tIWG );

    // aH( 4 ) == 0 is the NORMAL spare slot: accept
    EXPECT_TRUE( tIWG->restore_history_state( tH, 16, 2.5e-4 ) );

    // step count zero would silently ramp to BDF1 behind a full history
    EXPECT_FALSE( tIWG->restore_history_state( tH, 0, 2.5e-4 ) );

    // the last completed step must be positive and finite
    EXPECT_FALSE( tIWG->restore_history_state( tH, 16, 0.0 ) );
    EXPECT_FALSE( tIWG->restore_history_state( tH, 16, -1.0e-4 ) );
    EXPECT_FALSE( tIWG->restore_history_state( tH, 16,
        std::numeric_limits< real >::quiet_NaN() ) );

    // wrong length
    Vector< real > tShort( 3, 1.0e-4 );
    EXPECT_FALSE( tIWG->restore_history_state( tShort, 16, 2.5e-4 ) );

    // NaN in the history
    Vector< real > tNaN = tH ;
    tNaN( 1 ) = std::numeric_limits< real >::quiet_NaN();
    EXPECT_FALSE( tIWG->restore_history_state( tNaN, 16, 2.5e-4 ) );

    // a REQUIRED slot at zero: BDF5 with step count 16 needs aH( 0..3 ) > 0
    Vector< real > tHole = tH ;
    tHole( 2 ) = 0.0 ;
    EXPECT_FALSE( tIWG->restore_history_state( tHole, 16, 2.5e-4 ) );

    // negative history entry
    Vector< real > tNeg = tH ;
    tNeg( 1 ) = -1.0e-4 ;
    EXPECT_FALSE( tIWG->restore_history_state( tNeg, 16, 2.5e-4 ) );

    // +Inf in the history and in the last step
    Vector< real > tInf = tH ;
    tInf( 0 ) = std::numeric_limits< real >::infinity();
    EXPECT_FALSE( tIWG->restore_history_state( tInf, 16, 2.5e-4 ) );
    EXPECT_FALSE( tIWG->restore_history_state( tH, 16,
        std::numeric_limits< real >::infinity() ) );

    // rejects leave the ramp state untouched, checked for EVERY gate on a
    // fresh object ( step count, last dt, length, NaN entry, required hole )
    IWG_Timestep * tFresh = make_iwg();
    tFresh->set_timestepping_method( EulerMethod::BackwardDifference5 );
    BdfCoefficientProbe tProbeF( *tFresh );
    EXPECT_FALSE( tFresh->restore_history_state( tH, 0, 2.5e-4 ) );
    EXPECT_FALSE( tFresh->restore_history_state( tH, 16, 0.0 ) );
    EXPECT_FALSE( tFresh->restore_history_state( tShort, 16, 2.5e-4 ) );
    EXPECT_FALSE( tFresh->restore_history_state( tNaN, 16, 2.5e-4 ) );
    EXPECT_FALSE( tFresh->restore_history_state( tHole, 16, 2.5e-4 ) );
    EXPECT_EQ( tProbeF.step_count(), 0u );
    EXPECT_EQ( tProbeF.delta_time(), 1.0 );   // the ctor default survives
    for ( uint k = 0; k < 5; ++k )
    {
        EXPECT_EQ( tProbeF.h( k ), 0.0 );
    }

    // a successful restore caps the counter at the configured order ( the
    // ramp saturates there; a huge dumped counter carries no information )
    IWG_Timestep * tCap = make_iwg();
    tCap->set_timestepping_method( EulerMethod::BackwardDifference5 );
    BdfCoefficientProbe tProbeC( *tCap );
    EXPECT_TRUE( tCap->restore_history_state( tH, 100000u, 2.5e-4 ) );
    EXPECT_EQ( tProbeC.step_count(), 5u );
    EXPECT_EQ( tProbeC.delta_time(), 2.5e-4 );

    // non-BDF2-5 methods accept the state harmlessly ( nothing reads mH )
    IWG_Timestep * tEuler = make_iwg();
    tEuler->set_timestepping_method( EulerMethod::BackwardDifference1 );
    EXPECT_TRUE( tEuler->restore_history_state( tH, 16, 2.5e-4 ) );
    tEuler->delta_time() = 3.0e-4 ;
    tEuler->compute_bdf_coefficients();
    BdfCoefficientProbe tProbeE( *tEuler );
    EXPECT_EQ( tProbeE.alpha(), 1.0 );        // BDF1 alpha untouched

    delete tIWG ;
    delete tFresh ;
    delete tCap ;
    delete tEuler ;
}
