/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * Phase handling of the periodic source functions.
 *
 * set_periodic() used to follow its own correct set_phase() call with
 * set_time_offset( (aPhase/2*pi)*aPeriod ). Left-to-right evaluation makes
 * that (phase/2)*pi*T instead of phase*T/(2*pi), and because
 * set_time_offset() back-computes the phase from the offset, it also
 * overwrote the phase with phase*pi^2. Both slots were wrong by the same
 * factor pi^2 ~ 9.8696, and the corrupted phase reached the sine too --
 * function_sine() reads phase(), not time_offset().
 *
 * The first group of tests pins the two slots and the sine sample.
 *
 * The second group covers the fmod by-catch that the phase work made
 * observable:
 * std::fmod keeps the sign of its dividend, so before the fix every waveform
 * that wraps with it read a NEGATIVE tau at t < time_offset() instead of the
 * wrapped positive one. Square was therefore stuck at +A for all t below the
 * offset -- not a delayed square, a constant -- sawtooth ran negative, and
 * triangle escaped its own amplitude bound once the offset passed 3T/4.
 *
 * The third group covers the phase-sign convention. BELFEM follows
 * the SPICE definition: a POSITIVE phase ADVANCES every periodic waveform,
 * so the shape is evaluated at t + time_offset. Before that convention
 * landed, the sine
 * advanced ( it reads phase() through sin( omega*t + phase ) ) while the
 * other three delayed ( they read t - time_offset ) -- the same deck key
 * moving two ways at once.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_SourceFunction.hpp"

using belfem::real;
using belfem::uint;
using belfem::SourceFunction;
using belfem::SourceFunctionType;

//-----------------------------------------------------------------------------

// A phase must land as the matching fraction of the period: phase/(2*pi)*T.
// With T = 1 s and phase = pi/2 that is 0.25 s; the defect produced
// (pi/2)/2*pi*1 = pi^2/4 = 2.4674 s and a phase of pi^3/2 = 15.503 rad.
TEST( SourceFunction, PhaseToTimeOffsetSquare )
{
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Square,
                       1.0,
                       1.0,
                       belfem::constant::pi * 0.5 );

    EXPECT_NEAR( tFun.time_offset(), 0.25, 1e-12 );

    // the offset write-back must not corrupt the stored phase
    EXPECT_NEAR( tFun.phase(), belfem::constant::pi * 0.5, 1e-12 );
}

//-----------------------------------------------------------------------------

// The sine reads phase() rather than time_offset(), so the corrupted phase
// reached it as well: sin( pi^3/2 ) = 0.2034 instead of sin( pi/2 ) = 1.
// This is the assertion that refutes the "sine is unaffected" reading.
TEST( SourceFunction, PhaseReachesSine )
{
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Sine,
                       1.0,
                       1.0,
                       belfem::constant::pi * 0.5 );

    EXPECT_NEAR( tFun.phase(), belfem::constant::pi * 0.5, 1e-12 );

    // a quarter-period shift puts the sine at its maximum at t = 0
    EXPECT_NEAR( tFun.compute( 0.0 ), 1.0, 1e-12 );
}

//-----------------------------------------------------------------------------

// One conversion serves all four periodic waveforms, so pin it on all four.
// T = 4 s with phase = pi rad is half a period, i.e. 2 s.
TEST( SourceFunction, PhaseToTimeOffsetAllWaveforms )
{
    const SourceFunctionType tTypes[ 4 ] = { SourceFunctionType::Sine,
                                             SourceFunctionType::Square,
                                             SourceFunctionType::Triangle,
                                             SourceFunctionType::Sawtooth } ;

    for( uint k = 0 ; k < 4 ; ++k )
    {
        SourceFunction tFun ;
        tFun.set_periodic( tTypes[ k ], 1.0, 4.0, belfem::constant::pi );

        EXPECT_NEAR( tFun.time_offset(), 2.0, 1e-12 );
        EXPECT_NEAR( tFun.phase(), belfem::constant::pi, 1e-12 );
    }
}

//-----------------------------------------------------------------------------

// Regression guard: the zero-phase path was always correct and must stay so.
// Closed-form samples, not a comparison against the previous behaviour.
TEST( SourceFunction, ZeroPhaseUnchanged )
{
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Sine, 2.0, 1.0, 0.0 );

    EXPECT_NEAR( tFun.time_offset(), 0.0, 1e-15 );
    EXPECT_NEAR( tFun.phase(), 0.0, 1e-15 );

    // sin( 0 ) = 0 and sin( pi/2 ) = 1 at a quarter period
    EXPECT_NEAR( tFun.compute( 0.0 ), 0.0, 1e-12 );
    EXPECT_NEAR( tFun.compute( 0.25 ), 2.0, 1e-12 );
}

//-----------------------------------------------------------------------------

// A zero period drives set_period() into its omega = INFINITY branch, where
// set_phase() parks the offset at zero on purpose. The removed line called
// set_time_offset( 0 ), whose phase write-back then evaluated 0 * INFINITY
// and poisoned the phase with NaN. Recomputing the offset here -- as the
// originally proposed one-line repair would have -- keeps that NaN.
TEST( SourceFunction, ZeroPeriodLeavesNoNaN )
{
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Square,
                       1.0,
                       0.0,
                       belfem::constant::pi * 0.5 );

    EXPECT_FALSE( std::isnan( tFun.phase() ) );
    EXPECT_FALSE( std::isnan( tFun.time_offset() ) );

    EXPECT_NEAR( tFun.time_offset(), 0.0, 1e-15 );
    EXPECT_NEAR( tFun.phase(), belfem::constant::pi * 0.5, 1e-12 );
}

//-----------------------------------------------------------------------------
//  fmod sign -- the periodic extension before the offset
//-----------------------------------------------------------------------------

// DISCRIMINATOR. Before the fix fmod( -0.25, 1 ) = -0.25, which is < 0.5 and
// took the +A branch; the wrapped tau is 0.75 and must take the -A branch.
TEST( SourceFunction, SquareWrapsBelowZero )
{
    // a NEGATIVE phase is what drives the shifted time below zero once the
    // waveform advances, which is where fmod hands back a negative remainder
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Square, 1.0, 1.0,
                       -belfem::constant::pi * 0.5 );

    EXPECT_NEAR( tFun.compute( 0.0 ), -1.0, 1e-12 );
}

//-----------------------------------------------------------------------------

// DISCRIMINATOR. Before the fix this returned -0.25, outside the sawtooth's
// own [0,A) range.
TEST( SourceFunction, SawtoothWrapsBelowZero )
{
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Sawtooth, 1.0, 1.0,
                       -belfem::constant::pi * 0.5 );

    EXPECT_NEAR( tFun.compute( 0.0 ), 0.75, 1e-12 );
}

//-----------------------------------------------------------------------------

// DISCRIMINATOR. The +3T/4 shift inside function_triangle only raises the
// threshold, it does not wrap: once |offset| exceeds 3T/4 the negative tau
// escapes the fabs reduction and the wave leaves its declared amplitude.
// The +1.8 overshoot on an amplitude of 1 was measured before the fmod fix,
// with the mirror-image case of a POSITIVE phase under the old delay
// convention; under the advance convention a negative phase is what reaches
// it. Against an intermediate tree ( fmod wrapped, advance convention not
// yet applied ) this same call returns -0.2, so the assertion discriminates
// the sign flip.
TEST( SourceFunction, TriangleStaysInsideItsAmplitude )
{
    SourceFunction tFun ;
    tFun.set_periodic( SourceFunctionType::Triangle, 1.0, 1.0,
                       -1.9 * belfem::constant::pi );

    EXPECT_NEAR( tFun.time_offset(), -0.95, 1e-12 );
    EXPECT_NEAR( tFun.compute( 0.0 ), 0.2, 1e-12 );
}

//-----------------------------------------------------------------------------

// GUARD. A zero phase is the one case invariant under BOTH the wrap and the
// sign convention, so these closed-form samples must never move.
TEST( SourceFunction, ZeroPhaseWaveformsUnmoved )
{
    SourceFunction tSquare, tSaw, tTri ;
    tSquare.set_periodic( SourceFunctionType::Square,   1.0, 1.0, 0.0 );
    tSaw.set_periodic(    SourceFunctionType::Sawtooth, 1.0, 1.0, 0.0 );
    tTri.set_periodic(    SourceFunctionType::Triangle, 1.0, 1.0, 0.0 );

    EXPECT_NEAR( tSquare.compute( 0.3 ),  1.0,  1e-12 );
    EXPECT_NEAR( tSquare.compute( 0.7 ), -1.0,  1e-12 );
    EXPECT_NEAR( tSaw.compute( 0.25 ),    0.25, 1e-12 );
    EXPECT_NEAR( tSaw.compute( 0.75 ),    0.75, 1e-12 );
    EXPECT_NEAR( tTri.compute( 0.0 ),     0.0,  1e-12 );
    EXPECT_NEAR( tTri.compute( 0.25 ),    1.0,  1e-12 );
}

//-----------------------------------------------------------------------------

// The property the whole row is about: these are periodic functions, so
// f(t) == f(t+T) everywhere, including across the offset seam. Mixed roles
// at this phase -- square is a guard at t=-0.6 and a discriminator at t=0.1,
// triangle the other way round, sawtooth a discriminator at both -- so both
// sample times are needed to cover all three waveforms.
TEST( SourceFunction, PeriodicExtensionHolds )
{
    const real tPeriod = 1.0 ;
    const SourceFunctionType tTypes[ 3 ] = { SourceFunctionType::Square,
                                             SourceFunctionType::Triangle,
                                             SourceFunctionType::Sawtooth } ;
    const real tTimes[ 2 ] = { -0.6, 0.1 } ;

    for( uint k = 0 ; k < 3 ; ++k )
    {
        SourceFunction tFun ;
        tFun.set_periodic( tTypes[ k ], 1.0, tPeriod,
                           belfem::constant::pi * 0.5 );

        for( uint i = 0 ; i < 2 ; ++i )
        {
            EXPECT_NEAR( tFun.compute( tTimes[ i ] ),
                         tFun.compute( tTimes[ i ] + tPeriod ), 1e-12 );
        }
    }
}

//-----------------------------------------------------------------------------
//  one phase-sign convention for the whole periodic family
//-----------------------------------------------------------------------------

// A positive phase ADVANCES ( SPICE ). A quarter-period phase must therefore
// give exactly the zero-phase waveform sampled a quarter period LATER, and it
// must do so for all four members -- that equality is what the sine satisfied
// and the other three did not.
TEST( SourceFunction, PositivePhaseAdvancesEveryWaveform )
{
    const real tPeriod  = 1.0 ;
    const real tQuarter = tPeriod * 0.25 ;
    const SourceFunctionType tTypes[ 4 ] = { SourceFunctionType::Sine,
                                             SourceFunctionType::Square,
                                             SourceFunctionType::Triangle,
                                             SourceFunctionType::Sawtooth } ;
    // sample away from the square's discontinuities at tau in { 0, 0.5 }
    const real tTimes[ 4 ] = { 0.05, 0.1, 0.4, 0.8 } ;

    for( uint k = 0 ; k < 4 ; ++k )
    {
        SourceFunction tShifted, tPlain ;
        tShifted.set_periodic( tTypes[ k ], 1.0, tPeriod,
                               belfem::constant::pi * 0.5 );
        tPlain.set_periodic( tTypes[ k ], 1.0, tPeriod, 0.0 );

        for( uint i = 0 ; i < 4 ; ++i )
        {
            EXPECT_NEAR( tShifted.compute( tTimes[ i ] ),
                         tPlain.compute( tTimes[ i ] + tQuarter ), 1e-12 );
        }
    }
}
