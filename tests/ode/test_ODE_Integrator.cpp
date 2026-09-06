/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for ODE base class, Integrator wrapper, RK45, and DOP853.
 * See: tests_12_ode.md §3–§6
 *
 * Deferred: event handling (check_events commented out in RK45),
 *           DOP853 via Integrator (not wired up), stiff systems, RK78.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_ODE.hpp"
#include "cl_ODE_Integrator.hpp"
#include "en_ODE_Type.hpp"
#include "en_ODE_Status.hpp"
#include "fn_ODE_RK45.hpp"
#include "fn_ODE_DOP853.hpp"

namespace
{
    const belfem::real tTol = 1e-6;      // tolerance for RK45 integration
    const belfem::real tTolDOP = 1e-10;  // tolerance for DOP853 integration

    // =========================================================================
    // Concrete ODE subclasses with known analytic solutions
    // =========================================================================

    // y' = -y,  y(0) = 1  →  y(t) = e^(-t)
    class ExponentialDecayODE : public belfem::ode::ODE
    {
    public:
        ExponentialDecayODE() : ODE( 1 ) {}

        void
        compute( const belfem::real           & aT,
                 const belfem::Vector< belfem::real > & aY,
                       belfem::Vector< belfem::real > & adYdT ) override
        {
            adYdT( 0 ) = -aY( 0 );
        }
    };

    // y' = cos(t),  y(0) = 0  →  y(t) = sin(t)
    class SinCosODE : public belfem::ode::ODE
    {
    public:
        SinCosODE() : ODE( 1 ) {}

        void
        compute( const belfem::real           & aT,
                 const belfem::Vector< belfem::real > & aY,
                       belfem::Vector< belfem::real > & adYdT ) override
        {
            adYdT( 0 ) = std::cos( aT );
        }
    };

    // y1' = y2,  y2' = -y1  (simple harmonic oscillator)
    // y1(0) = 1, y2(0) = 0  →  y1(t) = cos(t), y2(t) = -sin(t)
    class HarmonicODE : public belfem::ode::ODE
    {
    public:
        HarmonicODE() : ODE( 2 ) {}

        void
        compute( const belfem::real           & aT,
                 const belfem::Vector< belfem::real > & aY,
                       belfem::Vector< belfem::real > & adYdT ) override
        {
            adYdT( 0 ) =  aY( 1 );
            adYdT( 1 ) = -aY( 0 );
        }
    };

    // y' = c (constant), y(0) = 0 → y(t) = c*t
    // (idea from ChatGPT) — simplest possible ODE for smoke tests
    class ConstantODE : public belfem::ode::ODE
    {
        belfem::real mC;
    public:
        explicit ConstantODE( belfem::real aC ) : ODE( 1 ), mC( aC ) {}

        void
        compute( const belfem::real           & aT,
                 const belfem::Vector< belfem::real > & aY,
                       belfem::Vector< belfem::real > & adYdT ) override
        {
            adYdT( 0 ) = mC;
        }
    };

    // Helper: integrate to target time using Integrator, return final status
    belfem::ode::Status
    integrate_to( belfem::ode::Integrator & aIntegrator,
                  belfem::real            & aT,
                  belfem::Vector< belfem::real > & aY,
                  const belfem::real        aTargetTime )
    {
        aIntegrator.maxtime() = aTargetTime;
        belfem::ode::Status tStatus = belfem::ode::Status::OK;

        while( aT < aTargetTime && tStatus != belfem::ode::Status::MAXIT )
        {
            tStatus = aIntegrator.step( aT, aY );
            if( tStatus == belfem::ode::Status::TRAPPED )
            {
                break;
            }
        }
        return tStatus;
    }
}

// =============================================================================
// §3  ODE Base Class  [semantic]
// =============================================================================

TEST( ODEBase, DimensionAccessor )
{
    ExponentialDecayODE tODE1;
    EXPECT_EQ( tODE1.dimension(), 1u );

    HarmonicODE tODE2;
    EXPECT_EQ( tODE2.dimension(), 2u );
}

TEST( ODEBase, CheckEventsDefaultReturnsZero )
{
    ExponentialDecayODE tODE;

    belfem::real tT = 1.0;
    belfem::Vector< belfem::real > tY = { 0.5 };

    EXPECT_NEAR( tODE.check_events( tT, tY ), 0.0, 1e-15 );
}

// =============================================================================
// §4  RK45 Integration (via Integrator)  [semantic]
// =============================================================================

TEST( RK45, ExponentialDecay )
{
    // y' = -y, y(0) = 1 → y(1) = e^(-1)
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-8;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    integrate_to( tIntegrator, tT, tY, 1.0 );

    EXPECT_NEAR( tY( 0 ), std::exp( -1.0 ), tTol );
}

TEST( RK45, SinCosine )
{
    // y' = cos(t), y(0) = 0 → y(π/2) = 1
    SinCosODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-8;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 0.0 };

    belfem::real tTarget = 0.5 * std::acos( -1.0 );  // π/2
    integrate_to( tIntegrator, tT, tY, tTarget );

    EXPECT_NEAR( tY( 0 ), 1.0, tTol );
}

TEST( RK45, TwoDimensionalHarmonic )
{
    // y1' = y2, y2' = -y1, y(0) = (1, 0)
    // After t = 2π: y1 ≈ 1, y2 ≈ 0
    HarmonicODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-8;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0, 0.0 };

    belfem::real tTarget = 2.0 * std::acos( -1.0 );  // 2π
    integrate_to( tIntegrator, tT, tY, tTarget );

    EXPECT_NEAR( tY( 0 ),  1.0, tTol );
    EXPECT_NEAR( tY( 1 ),  0.0, tTol );
}

TEST( RK45, StatusOK )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.maxtime()  = 10.0;  // far away — won't trap

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    belfem::ode::Status tStatus = tIntegrator.step( tT, tY );
    EXPECT_EQ( tStatus, belfem::ode::Status::OK );
    EXPECT_GT( tT, 0.0 );
}

TEST( RK45, StatusTrapped )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 1.0;   // large step
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.maxtime()  = 0.05;  // small maxtime — will trap immediately

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    belfem::ode::Status tStatus = tIntegrator.step( tT, tY );
    EXPECT_EQ( tStatus, belfem::ode::Status::TRAPPED );
}

TEST( RK45, TrappedPreservesStoredTimestep )
{
    // (idea from ChatGPT) — when TRAPPED, aStep is NOT overwritten
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.2;
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.maxtime()  = 0.05;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    belfem::ode::Status tStatus = tIntegrator.step( tT, tY );
    EXPECT_EQ( tStatus, belfem::ode::Status::TRAPPED );

    // stored timestep should remain at original value
    // (line 1321: aStep = h only when Status::OK)
    EXPECT_NEAR( tIntegrator.timestep(), 0.2, 1e-15 );
}

TEST( RK45, FixedTimestep )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.01;
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.set_auto_timestep( false );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    integrate_to( tIntegrator, tT, tY, 1.0 );

    EXPECT_NEAR( tY( 0 ), std::exp( -1.0 ), tTol );
}

TEST( RK45, AutoTimestepOffLeavesStepUnchanged )
{
    // (idea from ChatGPT) — with auto off, timestep() stays constant
    // Step must be small enough that the error stays below epsilon
    // without adaptive refinement.
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.01;
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.maxtime()  = 10.0;
    tIntegrator.set_auto_timestep( false );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    ASSERT_EQ( tIntegrator.step( tT, tY ), belfem::ode::Status::OK );
    EXPECT_NEAR( tIntegrator.timestep(), 0.01, 1e-15 );
}

TEST( RK45, AutoTimestepOnCanAdaptStep )
{
    // (idea from ChatGPT) — with auto on, timestep() changes after a step
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.5;
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.maxtime()  = 10.0;
    tIntegrator.set_auto_timestep( true );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    ASSERT_EQ( tIntegrator.step( tT, tY ), belfem::ode::Status::OK );
    // timestep should have been adapted (changed from 0.5)
    EXPECT_GT( std::abs( tIntegrator.timestep() - 0.5 ), 1e-12 );
}

TEST( RK45, MaxIterationsZeroReturnsMaxit )
{
    // (idea from ChatGPT) — zero iterations → immediate MAXIT, state unchanged
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-15;
    tIntegrator.max_num_iterations() = 0;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    belfem::ode::Status tStatus = tIntegrator.step( tT, tY );
    EXPECT_EQ( tStatus, belfem::ode::Status::MAXIT );
    EXPECT_NEAR( tT, 0.0, 1e-15 );
    EXPECT_NEAR( tY( 0 ), 1.0, 1e-15 );
}

// =============================================================================
// §5  DOP853 Integration (Direct Call)  [semantic]
// =============================================================================

TEST( DOP853, ExponentialDecay )
{
    // y' = -y, y(0) = 1 → y(1) = e^(-1)
    ExponentialDecayODE tODE;

    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    belfem::ode::DOP853_init( tODE, tWork );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };
    belfem::real tStep = 0.1;
    belfem::real tMaxTime = 1.0;

    belfem::ode::Status tStatus = belfem::ode::Status::OK;
    while( tT < tMaxTime && tStatus != belfem::ode::Status::MAXIT )
    {
        tStatus = belfem::ode::DOP853(
            tODE, tT, tY, tStep, tWork,
            1e-12, 1000, tMaxTime, true );

        if( tStatus == belfem::ode::Status::TRAPPED )
        {
            break;
        }
    }

    EXPECT_NEAR( tY( 0 ), std::exp( -1.0 ), tTolDOP );
}

TEST( DOP853, SinCosine )
{
    // y' = cos(t), y(0) = 0 → y(π/2) = 1
    SinCosODE tODE;

    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    belfem::ode::DOP853_init( tODE, tWork );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 0.0 };
    belfem::real tStep = 0.1;
    belfem::real tTarget = 0.5 * std::acos( -1.0 );

    belfem::ode::Status tStatus = belfem::ode::Status::OK;
    while( tT < tTarget && tStatus != belfem::ode::Status::MAXIT )
    {
        tStatus = belfem::ode::DOP853(
            tODE, tT, tY, tStep, tWork,
            1e-12, 1000, tTarget, true );

        if( tStatus == belfem::ode::Status::TRAPPED )
        {
            break;
        }
    }

    EXPECT_NEAR( tY( 0 ), 1.0, tTolDOP );
}

TEST( DOP853, HigherOrderAccuracy )
{
    // DOP853 should achieve tighter accuracy than RK45 on exponential decay
    ExponentialDecayODE tODE;

    // RK45 via Integrator
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );
    tIntegrator.timestep() = 0.5;  // deliberately large initial step
    tIntegrator.epsilon()  = 1e-8;

    belfem::real tT_rk = 0.0;
    belfem::Vector< belfem::real > tY_rk = { 1.0 };
    integrate_to( tIntegrator, tT_rk, tY_rk, 1.0 );
    belfem::real tError_rk = std::abs( tY_rk( 0 ) - std::exp( -1.0 ) );

    // DOP853 direct
    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    belfem::ode::DOP853_init( tODE, tWork );

    belfem::real tT_dop = 0.0;
    belfem::Vector< belfem::real > tY_dop = { 1.0 };
    belfem::real tStep = 0.5;

    belfem::ode::Status tStatus = belfem::ode::Status::OK;
    while( tT_dop < 1.0 && tStatus != belfem::ode::Status::MAXIT )
    {
        tStatus = belfem::ode::DOP853(
            tODE, tT_dop, tY_dop, tStep, tWork,
            1e-12, 1000, 1.0, true );
        if( tStatus == belfem::ode::Status::TRAPPED ) break;
    }
    belfem::real tError_dop = std::abs( tY_dop( 0 ) - std::exp( -1.0 ) );

    // DOP853 error should be smaller than RK45 error
    EXPECT_LT( tError_dop, tError_rk );
}

TEST( DOP853, StatusTrapped )
{
    ExponentialDecayODE tODE;

    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    belfem::ode::DOP853_init( tODE, tWork );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };
    belfem::real tStep = 1.0;
    belfem::real tMaxTime = 0.05;

    belfem::ode::Status tStatus = belfem::ode::DOP853(
        tODE, tT, tY, tStep, tWork,
        1e-8, 1000, tMaxTime, true );

    EXPECT_EQ( tStatus, belfem::ode::Status::TRAPPED );
}

TEST( DOP853, InitAllocatesFourteenWorkVectors )
{
    // (idea from ChatGPT) — verify work cell size and vector dimensions
    HarmonicODE tODE;

    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    belfem::ode::DOP853_init( tODE, tWork );

    ASSERT_EQ( tWork.size(), 14u );
    for( belfem::index_t k = 0; k < tWork.size(); ++k )
    {
        EXPECT_EQ( tWork( k ).length(), 2u );  // HarmonicODE dimension = 2
    }
}

TEST( DOP853, MaxIterationsZeroReturnsMaxit )
{
    // (idea from ChatGPT) — direct DOP853 with maxIter=0
    ExponentialDecayODE tODE;

    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    belfem::ode::DOP853_init( tODE, tWork );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };
    belfem::real tStep = 0.1;

    belfem::ode::Status tStatus = belfem::ode::DOP853(
        tODE, tT, tY, tStep, tWork,
        1e-12, 0, 10.0, true );

    EXPECT_EQ( tStatus, belfem::ode::Status::MAXIT );
    EXPECT_NEAR( tT, 0.0, 1e-15 );
    EXPECT_NEAR( tY( 0 ), 1.0, 1e-15 );
}

#ifndef NDEBUG
TEST( DOP853, UndersizedWorkCellThrows )
{
    // (idea from ChatGPT) — work cell with 13 entries triggers BELFEM_ASSERT
    ExponentialDecayODE tODE;

    belfem::Cell< belfem::Vector< belfem::real > > tWork;
    tWork.set_size( 13, belfem::Vector< belfem::real >( 1, 0.0 ) );

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };
    belfem::real tStep = 0.1;

    EXPECT_THROW(
        belfem::ode::DOP853( tODE, tT, tY, tStep, tWork ),
        std::runtime_error );
}
#endif

// =============================================================================
// §5b  DOP853 via Integrator  [semantic]
// =============================================================================

TEST( DOP853, ViaIntegratorExponentialDecay )
{
    // y' = -y, y(0) = 1 → y(1) = e^(-1)
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::DOP853 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-12;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    integrate_to( tIntegrator, tT, tY, 1.0 );

    EXPECT_NEAR( tY( 0 ), std::exp( -1.0 ), tTolDOP );
}

TEST( DOP853, ViaIntegratorHarmonic )
{
    // y1' = y2, y2' = -y1, y(0) = (1, 0) → after 2π: (1, 0)
    HarmonicODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::DOP853 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-12;

    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0, 0.0 };

    belfem::real tTarget = 2.0 * std::acos( -1.0 );
    integrate_to( tIntegrator, tT, tY, tTarget );

    EXPECT_NEAR( tY( 0 ),  1.0, tTolDOP );
    EXPECT_NEAR( tY( 1 ),  0.0, tTolDOP );
}

// =============================================================================
// §6  Integrator Configuration  [semantic]
// =============================================================================

TEST( Integrator, TypeAccessorRK45 )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    EXPECT_EQ( tIntegrator.type(), belfem::ode::Type::RK45 );
}

TEST( Integrator, TypeAccessorDOP853 )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::DOP853 );

    EXPECT_EQ( tIntegrator.type(), belfem::ode::Type::DOP853 );
}

TEST( Integrator, EpsilonMutableRef )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.epsilon() = 1e-10;
    EXPECT_NEAR( tIntegrator.epsilon(), 1e-10, 1e-15 );
}

TEST( Integrator, TimestepMutableRef )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.01;
    EXPECT_NEAR( tIntegrator.timestep(), 0.01, 1e-15 );
}

TEST( Integrator, MaxtimeMutableRef )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.maxtime() = 100.0;
    EXPECT_NEAR( tIntegrator.maxtime(), 100.0, 1e-15 );
}

TEST( Integrator, MaxIterationsMutableRef )
{
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.max_num_iterations() = 500;
    EXPECT_EQ( tIntegrator.max_num_iterations(), 500u );
}

TEST( Integrator, TimeDoesNotDriveStep )
{
    // time() does not drive integration — step() uses the caller's aT.
    // But time() is synced to aT after the step returns.
    ExponentialDecayODE tODE;
    belfem::ode::Integrator tIntegrator( tODE, belfem::ode::Type::RK45 );

    tIntegrator.timestep() = 0.1;
    tIntegrator.epsilon()  = 1e-8;
    tIntegrator.maxtime()  = 10.0;

    // set time() to something unrelated
    tIntegrator.time() = 999.0;

    // pass aT = 0.0 into step — the step should use aT, not time()
    belfem::real tT = 0.0;
    belfem::Vector< belfem::real > tY = { 1.0 };

    belfem::ode::Status tStatus = tIntegrator.step( tT, tY );

    // step should have advanced tT from 0, not from 999
    EXPECT_GT( tT, 0.0 );
    EXPECT_LT( tT, 2.0 );
    EXPECT_EQ( tStatus, belfem::ode::Status::OK );

    // time() should now reflect the actual integration time
    EXPECT_NEAR( tIntegrator.time(), tT, 1e-15 );
}

// =============================================================================
// §6.2  Integrator Error Paths
// =============================================================================

TEST( IntegratorError, UnknownTypeThrows )
{
    // (idea from Gemini) — BELFEM_ERROR (always active)
    ExponentialDecayODE tODE;

    EXPECT_THROW(
        belfem::ode::Integrator tIntegrator(
            tODE,
            static_cast< belfem::ode::Type >( 99 ) ),
        std::runtime_error );
}
