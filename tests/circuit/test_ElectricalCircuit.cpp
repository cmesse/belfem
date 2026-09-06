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
 * Solve-level tests for ElectricalCircuit timestep semantics:
 *   - mixed switch-before-V dof counting
 *   - switch latch revert on a rejected timestep
 *   - warm-restart state surviving the first stamp ( unit-level
 *     discriminator only — the coupled restart gate stays owed )
 *   - order-2 RLC time integration across a forced rejection ( library
 *     semantics; the demo executable itself is not exercised here )
 *   - the demo circuits harvested as fixed-timestep transients with
 *     analytic gates: resistive ( exact, pins the branch-current sign ),
 *     Frederic's switched RLC ( two regimes ), the diode bridge ( first
 *     transient exercise of the diode Newton ), and a netlist twin of
 *     Frederic's circuit tying NgspiceCircuitFactory to the same gate
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <memory>
#include <string>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_ElectricalCircuit.hpp"
#include "cl_Component.hpp"
#include "cl_SourceFunction.hpp"
#include "cl_NgspiceCircuitFactory.hpp"
#include "constants.hpp"
#include "fn_norm.hpp"

#ifdef BELFEM_HDF5
#include "cl_HDF5.hpp"
#include "hdf5_tools.hpp"
#include "filetools.hpp"
#endif

using belfem::real;
using belfem::uint;
using belfem::index_t;
using belfem::Vector;
using belfem::SourceFunction;
using belfem::SourceFunctionType;
using belfem::electronics::ElectricalCircuit;
using belfem::electronics::NgspiceCircuitFactory;

namespace
{
//-----------------------------------------------------------------------------

    // Newton driver mirroring the FEM controller loop: the caller has already
    // run set_timestep() + shift() + compute_MNA_matrix() for this attempt.
    // aNumIterations, when given, receives the iteration count so a test
    // can report Newton health ( no relaxation: omega stays 1 )
    bool
    solve_attempt( ElectricalCircuit & aCircuit,
                   const uint          aMaxIt         = 50,
                   const real          aTolerance     = 1e-9,
                   uint              * aNumIterations = nullptr )
    {
        aCircuit.set_omega( 1.0 );
        real tEpsilon = BELFEM_REAL_MAX;
        uint tIt = 0;
        while ( tIt < aMaxIt && tEpsilon > aTolerance )
        {
            ++tIt;
            aCircuit.compute_jacobian_and_rhs();
            aCircuit.solve();
            tEpsilon = aCircuit.residual();
        }
        if ( aNumIterations != nullptr )
        {
            *aNumIterations = tIt;
        }
        return tEpsilon <= aTolerance;
    }

    // one accepted controller-parity step
    bool
    take_step( ElectricalCircuit & aCircuit,
               const real          aDeltaTime,
               uint              * aNumIterations = nullptr )
    {
        aCircuit.set_timestep( aDeltaTime );
        aCircuit.shift();
        aCircuit.compute_MNA_matrix();
        return solve_attempt( aCircuit, 50, 1e-9, aNumIterations );
    }

    // tolerance with an absolute floor: every trace below passes through
    // zero, where a purely relative comparison divides by roundoff
    real
    mixed_tol( const real aWant, const real aRel )
    {
        return aRel * std::max( 1.0, std::abs( aWant ) );
    }

    // write aContent to an absolute scratch path and return it ( same
    // helper as test_HybridFactory.cpp; not an HDF5 FileGuard, so it
    // compiles on a USE_HDF5=OFF tree )
    std::string
    write_scratch( const std::string & aName, const std::string & aContent )
    {
        const std::string tPath = ::testing::TempDir() + aName;
        std::ofstream tFile( tPath );
        tFile << aContent;
        return tPath;
    }

//-----------------------------------------------------------------------------

    // series RLC at resonance, driven by a 1 V sine ( demo's second circuit ):
    // V( 0-gnd ), L( 0-1 ), C( 1-2 ), R( 2-gnd ), order 2
    std::unique_ptr< ElectricalCircuit >
    create_rlc_circuit( const real aFreq )
    {
        std::unique_ptr< ElectricalCircuit > tCircuit(
                new ElectricalCircuit( 4 ) );
        SourceFunction * tFun = new SourceFunction();
        tFun->set_periodic( SourceFunctionType::Sine, 1.0, 1.0 / aFreq, 0.0 );
        tCircuit->create_voltage_source( tFun, 0, 3, "V1" );
        tCircuit->create_inductor( 15.915e-6, 2, 0, 1, "L1" );
        tCircuit->create_capacitor( 15.915e-6, 2, 1, 2, "C1" );
        tCircuit->create_resistor( 1.0, 2, 3, "R1" );
        return tCircuit;
    }

//-----------------------------------------------------------------------------

    // Frederic's circuit ( demo's first circuit ), 1 V / 50 Hz sine:
    // V( 0-gnd ), L( 0-1 ), R( 1-gnd ), S( 1-2 ), C( 2-gnd ), ground = 3,
    // L and C at BDF order 3.
    //
    // The voltage source is created before the switch; the regression
    // SwitchBeforeVoltageSourceAgrees is what catches a reorder.
    //
    // Values are written in the arithmetic form the SPICE number parser
    // produces ( "200u" is 200 * 1e-6, "SIN(0 1 50)" is 1.0 / 50 ), so the
    // netlist twin in NetlistTwinMatchesFredericCircuit reproduces them bit
    // for bit. Labels are lowercase because the netlist parser folds case.
    //
    // t_switch sits half a step BEFORE the demo's 25 ms: Switch::shift()
    // compares the accumulated end-of-step time against it, and a value on
    // an exact step boundary is a floating-point tie that fires at step 500
    // or 501 depending on the build. The first end-of-step time
    // >= 24.975 ms is still 25.0 ms, the end of step 500 at 50 us.
    const real gFredericFreq = 50.0;

    std::unique_ptr< ElectricalCircuit >
    create_frederic_circuit()
    {
        std::unique_ptr< ElectricalCircuit > tCircuit(
                new ElectricalCircuit( 4 ) );
        SourceFunction * tFun = new SourceFunction();
        tFun->set_periodic( SourceFunctionType::Sine, 1.0, 1.0 / gFredericFreq, 0.0 );
        tCircuit->create_voltage_source( tFun, 0, 3, "v1" );
        tCircuit->create_inductor( 200 * 1e-6, 3, 0, 1, "l1" );
        tCircuit->create_resistor( 1.0, 1, 3, "r1" );
        tCircuit->create_switch( false, 24.975 * 1e-3, 1, 2, "s1" );
        tCircuit->create_capacitor( 3.183 * 1e-3, 3, 2, 3, "c1" );
        return tCircuit;
    }

//-----------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------
// a switch created before a voltage source must not shift the
// V-branch MNA stamps onto the switch's unknown-current dof
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, SwitchBeforeVoltageSourceAgrees )
{
    // divider: V1 ( 5 V ) at node 0, R1 node0-node1, R2 node1-ground,
    // open switch node1-ground. Switch and V-source terminals differ, so
    // a misplaced V stamp lands outside the switch column's sparsity.
    // Expected solution: V( n0 ) = 5, V( n1 ) = 2.5, I_switch = 0.

    // creation order A: switch BEFORE the voltage source ( the failing layout )
    ElectricalCircuit tCircuitA( 3 );
    tCircuitA.create_switch( false, 1.0, 1, 2, "S1" );
    SourceFunction * tFunA = new SourceFunction();
    tFunA->set_constant( 5.0 );
    tCircuitA.create_voltage_source( tFunA, 0, 2, "V1" );
    tCircuitA.create_resistor( 1000.0, 0, 1, "R1" );
    tCircuitA.create_resistor( 1000.0, 1, 2, "R2" );

    // creation order B: voltage source first ( the always-legal layout )
    ElectricalCircuit tCircuitB( 3 );
    SourceFunction * tFunB = new SourceFunction();
    tFunB->set_constant( 5.0 );
    tCircuitB.create_voltage_source( tFunB, 0, 2, "V1" );
    tCircuitB.create_switch( false, 1.0, 1, 2, "S1" );
    tCircuitB.create_resistor( 1000.0, 0, 1, "R1" );
    tCircuitB.create_resistor( 1000.0, 1, 2, "R2" );

    real tDeltaTime = 1e-6;

    // on the pre-fix tree order A throws ( structural zero ) or solves wrong;
    // an uncaught throw fails the test either way
    ASSERT_TRUE( take_step( tCircuitA, tDeltaTime ) );
    ASSERT_TRUE( take_step( tCircuitB, tDeltaTime ) );

    // both orders must produce the same physics
    EXPECT_NEAR( tCircuitA.get_voltage_on_node( 0 ),
                 tCircuitB.get_voltage_on_node( 0 ), 1e-9 );
    EXPECT_NEAR( tCircuitA.get_voltage_on_node( 1 ),
                 tCircuitB.get_voltage_on_node( 1 ), 1e-9 );

    // and the right physics
    EXPECT_NEAR( tCircuitA.get_voltage_on_node( 0 ), 5.0, 1e-6 );
    EXPECT_NEAR( tCircuitA.get_voltage_on_node( 1 ), 2.5, 1e-6 );
}

//-----------------------------------------------------------------------------
// a switch that fired during a rejected timestep attempt must
// revert with shift_back() and fire again on the retry
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, SwitchLatchRevertsOnRejectedStep )
{
    // R and an open switch across node0-ground; t_switch inside the first step
    ElectricalCircuit tCircuit( 2 );
    tCircuit.create_resistor( 1000.0, 0, 1, "R1" );
    tCircuit.create_switch( false, 0.5e-6, 0, 1, "S1" );
    const belfem::electronics::Component * tSwitch = tCircuit.component( 1 );

    EXPECT_FALSE( tSwitch->is_closed() );

    // attempt crossing t_switch: the latch fires
    tCircuit.set_timestep( 1e-6 );
    tCircuit.shift();
    tCircuit.compute_MNA_matrix();
    EXPECT_TRUE( tSwitch->is_closed() );

    // the attempt is rejected: the latch must revert with everything else
    tCircuit.shift_back();
    EXPECT_FALSE( tSwitch->is_closed() );

    // the retry crosses t_switch again: the latch must fire again
    tCircuit.shift();
    EXPECT_TRUE( tSwitch->is_closed() );
}

//-----------------------------------------------------------------------------
// state restored by load_state() must survive the first
// compute_MNA_matrix(). Unit-level wipe discriminator only — NOT the
// coupled restart gate, which stays owed.
//-----------------------------------------------------------------------------

#ifdef BELFEM_HDF5

namespace
{
    void
    remove_file_if_exists( const std::string & aPath )
    {
        if ( belfem::file_exists( aPath ) )
        {
            std::remove( aPath.c_str() );
        }
    }

    // removes its file on scope exit, including gtest fatal-assert unwinding
    struct FileGuard
    {
        std::string mPath;

        explicit FileGuard( const std::string & aPath ) : mPath( aPath )
        {
            remove_file_if_exists( mPath );
        }

        ~FileGuard()
        {
            remove_file_if_exists( mPath );
        }
    };

    // V1 ( 5 V ) with a parallel resistor: one node dof + one branch dof
    std::unique_ptr< ElectricalCircuit >
    create_vr_circuit()
    {
        std::unique_ptr< ElectricalCircuit > tCircuit(
                new ElectricalCircuit( 2 ) );
        SourceFunction * tFun = new SourceFunction();
        tFun->set_constant( 5.0 );
        tCircuit->create_voltage_source( tFun, 0, 1, "V1" );
        tCircuit->create_resistor( 1000.0, 0, 1, "R1" );
        return tCircuit;
    }

    void
    save_circuit( ElectricalCircuit & aCircuit, const std::string & aPath )
    {
        belfem::HDF5 tFile( aPath, belfem::FileMode::NEW );
        tFile.create_group( "circuit" );
        aCircuit.save_state( tFile.active_group() );
    }

    void
    load_circuit( ElectricalCircuit & aCircuit, const std::string & aPath )
    {
        belfem::HDF5 tFile( aPath, belfem::FileMode::OPEN_RDONLY );
        tFile.select_group( "circuit" );
        aCircuit.load_state( tFile.active_group() );
    }

    void
    read_state_vectors( const std::string & aPath,
                        Vector< real >    & aX,
                        Vector< real >    & aPrevX )
    {
        belfem::HDF5 tFile( aPath, belfem::FileMode::OPEN_RDONLY );
        hid_t tGroup = tFile.select_group( "circuit" );
        herr_t tStatus = 0;
        belfem::hdf5::load_vector_from_file( tGroup, "x", aX, tStatus );
        belfem::hdf5::load_vector_from_file( tGroup, "prev_x", aPrevX, tStatus );
    }
}

TEST( ElectricalCircuit, WarmRestartKeepsState )
{
    FileGuard tPathA( "test_circuit_warm_restart_a.hdf5" );
    FileGuard tPathB( "test_circuit_warm_restart_b.hdf5" );

    // run two converged steps so BOTH dumped vectors are nonzero
    // ( after one step prev_x would still be the cold zeros and its
    // comparison below would be vacuous )
    {
        std::unique_ptr< ElectricalCircuit > tSource = create_vr_circuit();
        ASSERT_TRUE( take_step( *tSource, 1e-6 ) );
        ASSERT_TRUE( take_step( *tSource, 1e-6 ) );
        save_circuit( *tSource, tPathA.mPath );
    }

    Vector< real > tX0;
    Vector< real > tPrevX0;
    read_state_vectors( tPathA.mPath, tX0, tPrevX0 );
    ASSERT_GT( belfem::norm( tX0 ), 1e-6 );     // dump really is nonzero
    ASSERT_GT( belfem::norm( tPrevX0 ), 1e-6 ); // and so is the history

    // restart: load, stamp, dump again — the state must survive the stamp
    {
        std::unique_ptr< ElectricalCircuit > tRestart = create_vr_circuit();
        load_circuit( *tRestart, tPathA.mPath );
        tRestart->compute_MNA_matrix();
        save_circuit( *tRestart, tPathB.mPath );
    }

    Vector< real > tX1;
    Vector< real > tPrevX1;
    read_state_vectors( tPathB.mPath, tX1, tPrevX1 );

    ASSERT_EQ( tX1.length(), tX0.length() );
    ASSERT_EQ( tPrevX1.length(), tPrevX0.length() );
    for ( index_t k = 0; k < tX0.length(); ++k )
    {
        // pure serialization roundtrip: bitwise, not approximate
        EXPECT_EQ( tX1( k ), tX0( k ) );
        EXPECT_EQ( tPrevX1( k ), tPrevX0( k ) );
    }
}

TEST( ElectricalCircuit, WarmRestartControllerParity )
{
    // same discriminator through the controller call order:
    // load_state, then set_timestep + shift + compute_MNA_matrix.
    // shift() copies the loaded x into prev_x; the stamp must keep it.
    FileGuard tPathA( "test_circuit_warm_parity_a.hdf5" );
    FileGuard tPathB( "test_circuit_warm_parity_b.hdf5" );

    {
        std::unique_ptr< ElectricalCircuit > tSource = create_vr_circuit();
        ASSERT_TRUE( take_step( *tSource, 1e-6 ) );
        save_circuit( *tSource, tPathA.mPath );
    }

    Vector< real > tX0;
    Vector< real > tPrevX0;
    read_state_vectors( tPathA.mPath, tX0, tPrevX0 );

    {
        std::unique_ptr< ElectricalCircuit > tRestart = create_vr_circuit();
        load_circuit( *tRestart, tPathA.mPath );
        tRestart->set_timestep( 1e-6 );
        tRestart->shift();
        tRestart->compute_MNA_matrix();
        save_circuit( *tRestart, tPathB.mPath );
    }

    Vector< real > tX1;
    Vector< real > tPrevX1;
    read_state_vectors( tPathB.mPath, tX1, tPrevX1 );

    // after the shift, prev_x must carry the state loaded into x
    ASSERT_EQ( tPrevX1.length(), tX0.length() );
    for ( index_t k = 0; k < tX0.length(); ++k )
    {
        EXPECT_NEAR( tPrevX1( k ), tX0( k ), 1e-12 );
    }
}

//-----------------------------------------------------------------------------
// /circuit format v2: a restart must continue the run, not restart the
// physics. The traces below are the unit-level core of the R6 gate ( the
// coupled hphirun A/B remains owed — the controller's restart-entry dt cap
// is outside the circuit's control ).
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, WarmRestartContinuationEquivalence )
{
    const real tFreq = 1e4;
    const real tDeltaTime = 1e-7;
    const uint tNumBefore = 500;           // full order-2 registers at the dump
    const uint tNumAfter  = 500;
    const uint tInductor  = 1;
    const uint tCapacitorC = 2;

    FileGuard tPath( "test_circuit_continuation.hdf5" );

    // twin: uninterrupted, record the continuation window
    std::unique_ptr< ElectricalCircuit > tTwin = create_rlc_circuit( tFreq );
    Vector< real > tCurrentLA( tNumAfter, 0.0 );
    Vector< real > tCurrentCA( tNumAfter, 0.0 );
    Vector< real > tVoltageA( tNumAfter, 0.0 );
    for ( uint k = 0; k < tNumBefore + tNumAfter; ++k )
    {
        ASSERT_TRUE( take_step( *tTwin, tDeltaTime ) );
        if ( k >= tNumBefore )
        {
            tCurrentLA( k - tNumBefore ) = tTwin->get_current_on_component( tInductor );
            tCurrentCA( k - tNumBefore ) = tTwin->get_current_on_component( tCapacitorC );
            tVoltageA( k - tNumBefore ) = tTwin->get_voltage_on_node( 1 );
        }
    }

    // dumped run
    {
        std::unique_ptr< ElectricalCircuit > tSource = create_rlc_circuit( tFreq );
        for ( uint k = 0; k < tNumBefore; ++k )
        {
            ASSERT_TRUE( take_step( *tSource, tDeltaTime ) );
        }
        save_circuit( *tSource, tPath.mPath );
    }

    // restarted run
    std::unique_ptr< ElectricalCircuit > tRestart = create_rlc_circuit( tFreq );
    load_circuit( *tRestart, tPath.mPath );
    for ( uint k = 0; k < tNumAfter; ++k )
    {
        ASSERT_TRUE( take_step( *tRestart, tDeltaTime ) );
        EXPECT_NEAR( tRestart->get_current_on_component( tInductor ),
                     tCurrentLA( k ), 1e-12 );
        EXPECT_NEAR( tRestart->get_current_on_component( tCapacitorC ),
                     tCurrentCA( k ), 1e-12 );
        EXPECT_NEAR( tRestart->get_voltage_on_node( 1 ),
                     tVoltageA( k ), 1e-12 );
    }
}

TEST( ElectricalCircuit, WarmRestartPartialHistory )
{
    // dump after ONE accepted order-2 step: the registers are partially
    // filled ( size 1 < capacity 2 ), the fill state the big test never hits
    const real tFreq = 1e4;
    const real tDeltaTime = 1e-7;
    const uint tInductor = 1;

    FileGuard tPath( "test_circuit_partial_history.hdf5" );

    std::unique_ptr< ElectricalCircuit > tTwin = create_rlc_circuit( tFreq );
    ASSERT_TRUE( take_step( *tTwin, tDeltaTime ) );
    ASSERT_TRUE( take_step( *tTwin, tDeltaTime ) );

    {
        std::unique_ptr< ElectricalCircuit > tSource = create_rlc_circuit( tFreq );
        ASSERT_TRUE( take_step( *tSource, tDeltaTime ) );
        save_circuit( *tSource, tPath.mPath );
    }

    std::unique_ptr< ElectricalCircuit > tRestart = create_rlc_circuit( tFreq );
    load_circuit( *tRestart, tPath.mPath );
    ASSERT_TRUE( take_step( *tRestart, tDeltaTime ) );

    EXPECT_NEAR( tRestart->get_current_on_component( tInductor ),
                 tTwin->get_current_on_component( tInductor ), 1e-12 );
    EXPECT_NEAR( tRestart->get_voltage_on_node( 1 ),
                 tTwin->get_voltage_on_node( 1 ), 1e-12 );
}

TEST( ElectricalCircuit, WarmRestartSwitchLatch )
{
    // a switch that fired before the dump must come back fired WITHOUT
    // relying on the re-fire, and a further shift must not toggle it again
    FileGuard tPath( "test_circuit_switch_latch.hdf5" );

    {
        ElectricalCircuit tSource( 2 );
        tSource.create_resistor( 1000.0, 0, 1, "R1" );
        tSource.create_switch( false, 0.5e-6, 0, 1, "S1" );
        tSource.set_timestep( 1e-6 );
        tSource.shift() ;                  // crosses t_switch: latch fires
        tSource.compute_MNA_matrix() ;     // sizes the state vectors
        ASSERT_TRUE( tSource.component( 1 )->is_closed() );
        save_circuit( tSource, tPath.mPath );
    }

    ElectricalCircuit tRestart( 2 );
    tRestart.create_resistor( 1000.0, 0, 1, "R1" );
    tRestart.create_switch( false, 0.5e-6, 0, 1, "S1" );
    load_circuit( tRestart, tPath.mPath );

    // fired state restored before any shift
    EXPECT_TRUE( tRestart.component( 1 )->is_closed() );

    // a further shift past t_switch must not toggle again
    tRestart.set_timestep( 1e-6 );
    tRestart.shift() ;
    EXPECT_TRUE( tRestart.component( 1 )->is_closed() );
}

TEST( ElectricalCircuit, WarmRestartBeforeFirstShift )
{
    // a dump taken before any shift ( empty histories ) must load into an
    // object that is coherent for a stamp without a prior shift: the loader
    // seeds the BDF1 companions from the restored delta_time
    const real tFreq = 1e4;
    const real tDeltaTime = 1e-7;

    FileGuard tPath( "test_circuit_pre_shift.hdf5" );

    {
        std::unique_ptr< ElectricalCircuit > tSource = create_rlc_circuit( tFreq );
        tSource->set_timestep( tDeltaTime );
        tSource->compute_MNA_matrix() ;    // sizes the vectors; no shift yet
        save_circuit( *tSource, tPath.mPath );
    }

    // twin: the never-dumped cold circuit, first step
    std::unique_ptr< ElectricalCircuit > tTwin = create_rlc_circuit( tFreq );
    tTwin->set_timestep( tDeltaTime );
    tTwin->compute_MNA_matrix() ;
    ASSERT_TRUE( take_step( *tTwin, tDeltaTime ) );

    // restart: load the empty-history dump, then take the same first step
    std::unique_ptr< ElectricalCircuit > tRestart = create_rlc_circuit( tFreq );
    load_circuit( *tRestart, tPath.mPath );
    ASSERT_TRUE( take_step( *tRestart, tDeltaTime ) );

    EXPECT_NEAR( tRestart->get_current_on_component( 1 ),
                 tTwin->get_current_on_component( 1 ), 1e-12 );
    EXPECT_NEAR( tRestart->get_voltage_on_node( 1 ),
                 tTwin->get_voltage_on_node( 1 ), 1e-12 );
}

TEST( ElectricalCircuit, WarmRestartRefusesOldFormat )
{
    // a pre-v2 dump ( vectors only, no per-component data ) must be refused
    // rather than silently cold-starting the histories
    FileGuard tPath( "test_circuit_old_format.hdf5" );

    {
        belfem::HDF5 tFile( tPath.mPath, belfem::FileMode::NEW );
        tFile.create_group( "circuit" );
        hid_t tGroup = tFile.active_group();
        herr_t tStatus = 0;
        real tTime = 1e-6;
        real tDeltaTime = 1e-6;
        Vector< real > tX( 2, 1.0 );
        Vector< real > tPrevX( 2, 0.5 );
        belfem::hdf5::save_scalar_to_file( tGroup, "time", tTime, tStatus );
        belfem::hdf5::save_scalar_to_file( tGroup, "delta_time", tDeltaTime, tStatus );
        belfem::hdf5::save_vector_to_file( tGroup, "x", tX, tStatus );
        belfem::hdf5::save_vector_to_file( tGroup, "prev_x", tPrevX, tStatus );
    }

    std::unique_ptr< ElectricalCircuit > tCircuit = create_vr_circuit();
    EXPECT_THROW( load_circuit( *tCircuit, tPath.mPath ), std::runtime_error );
}

#endif // BELFEM_HDF5

//-----------------------------------------------------------------------------
// library semantics: an order-2 series RLC trace must not be
// poisoned by a rejected-and-retried timestep. Primary gate: agreement
// with a never-rejected twin on the common time grid after the retry.
// Analytic resonance amplitude is a loose sanity bound only.
//
// Pre-registered reading of a RED result here: suspected inductor
// current loss through ElectricalCircuit::shift_back() ( companions not
// reverted before compute_current() ) — a NEW finding to report, not a
// license to widen the fix under test.
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, RLCRingAcrossRejectedStep )
{
    const real tFreq = 1e4;
    const real tDeltaTime = 1e-7;          // 1000 steps per period
    const uint tNumSteps = 3000;           // three periods
    const uint tRejectAt = 1500;           // forced rejection mid-run
    const uint tInductor = 1;              // component index of L1

    // twin A: never rejected
    std::unique_ptr< ElectricalCircuit > tTwin = create_rlc_circuit( tFreq );
    Vector< real > tCurrentA( tNumSteps, 0.0 );
    for ( uint k = 0; k < tNumSteps; ++k )
    {
        ASSERT_TRUE( take_step( *tTwin, tDeltaTime ) );
        tCurrentA( k ) = tTwin->get_current_on_component( tInductor );
    }

    // run B: identical, but step tRejectAt is solved, rejected, and
    // re-integrated as two half steps ( the demo/controller reject path );
    // the grid realigns at the end of the second half step
    std::unique_ptr< ElectricalCircuit > tRun = create_rlc_circuit( tFreq );
    Vector< real > tCurrentB( tNumSteps, 0.0 );
    for ( uint k = 0; k < tNumSteps; ++k )
    {
        if ( k == tRejectAt )
        {
            ASSERT_TRUE( take_step( *tRun, tDeltaTime ) );
            tRun->shift_back();            // reject the converged attempt
            ASSERT_TRUE( take_step( *tRun, 0.5 * tDeltaTime ) );
            ASSERT_TRUE( take_step( *tRun, 0.5 * tDeltaTime ) );
        }
        else
        {
            ASSERT_TRUE( take_step( *tRun, tDeltaTime ) );
        }
        tCurrentB( k ) = tRun->get_current_on_component( tInductor );
    }

    // analytic sanity: at resonance the steady-state amplitude is V/R = 1 A
    real tPeak = 0.0;
    for ( uint k = tNumSteps - 1000; k < tNumSteps; ++k )
    {
        tPeak = std::max( tPeak, std::abs( tCurrentA( k ) ) );
    }
    EXPECT_NEAR( tPeak, 1.0, 0.05 );

    // primary gate: after the retry the two traces share the time grid and
    // must agree to far better than the truncation left by one halved step.
    // check immediately after the retry and over the remaining trace
    real tMaxDiff = 0.0;
    for ( uint k = tRejectAt; k < tNumSteps; ++k )
    {
        tMaxDiff = std::max( tMaxDiff,
                             std::abs( tCurrentB( k ) - tCurrentA( k ) ) );
    }
    EXPECT_LT( tMaxDiff, 1e-3 ); // vs a 1 A carrier: truncation-level only
}

//-----------------------------------------------------------------------------
// resistive: 1 V / 50 Hz sine across 1 ohm, node 1 = ground. No
// discretization error exists, so the gates are exact ( absolute floor
// because the window contains the sine's zeros ). The source sample is
// read from the component, not re-derived, so a source-function defect
// and a stamping defect cannot cancel. This is the test that pins the
// MNA branch-current convention: node rows are "current leaving the
// node", so the source branch current is the NEGATIVE of the delivered
// current while the resistor's is positive. Assert the sign; do not
// absorb it into std::abs.
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, ResistiveSinePinsCurrentConvention )
{
    const real tFreq = 50.0;
    const real tR = 1.0;
    const real tDeltaTime = 50e-6;         // 400 steps per period
    const uint tNumSteps = 400;            // one period
    const uint tSource = 0;
    const uint tResistor = 1;

    ElectricalCircuit tCircuit( 2 );       // ground is the LAST index
    SourceFunction * tFun = new SourceFunction();
    tFun->set_periodic( SourceFunctionType::Sine, 1.0, 1.0 / tFreq, 0.0 );
    tCircuit.create_voltage_source( tFun, 0, 1, "V1" );
    tCircuit.create_resistor( tR, 0, 1, "R1" );

    for ( uint k = 0; k < tNumSteps; ++k )
    {
        ASSERT_TRUE( take_step( tCircuit, tDeltaTime ) ) << "step " << k + 1;

        // get_value() holds the sample written by the last shift(), i.e.
        // the end-of-step source value the stamp used
        const real tWant = tCircuit.component( tSource )->get_value();
        const real tV = tCircuit.get_voltage_on_node( 0 )
                      - tCircuit.get_voltage_on_node( 1 );
        EXPECT_NEAR( tV, tWant, mixed_tol( tWant, 1e-12 ) ) << "step " << k + 1;

        const real tI = tV / tR;
        EXPECT_NEAR( tCircuit.get_current_on_component( tSource ),
                     -tI, mixed_tol( tI, 1e-12 ) ) << "step " << k + 1;
        EXPECT_NEAR( tCircuit.get_current_on_component( tResistor ),
                     tI, mixed_tol( tI, 1e-12 ) ) << "step " << k + 1;
    }
}

//-----------------------------------------------------------------------------
// Frederic's circuit, two regimes. Before the switch fires the source
// sees R + jwL; afterwards jwL + R || ( 1 / jwC ). Both amplitudes are
// computed from the live component values, not from literals, so a
// transcription slip cannot be the thing that passes. The post-switch
// amplitude is 51 % higher, so the gate fails loudly if the switch never
// fires or fires in the wrong regime. Windows are full 400-step periods:
// a half period can miss the peak depending on phase.
//
// Step numbering is 1-based in the comments ( step n ends at n * dt ).
// Step 500 is already post-switch: shift() latches the switch at the
// end-of-step time 25.0 ms and compute_MNA_matrix() then stamps it
// closed, so the pre-switch window stops at step 499.
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, FredericCircuitTwoRegimes )
{
    const real tDeltaTime = 50e-6;         // 400 steps per 50 Hz period
    const uint tNumSteps = 1600;           // 80 ms
    const uint tSource = 0;
    const uint tInductor = 1;
    const uint tResistor = 2;
    const uint tSwitch = 3;
    const uint tCapacitor = 4;

    std::unique_ptr< ElectricalCircuit > tCircuit = create_frederic_circuit();

    // expected amplitudes from the component values
    using cplx = std::complex< real >;
    const cplx tJ( 0.0, 1.0 );
    const real tOmega = 2.0 * belfem::constant::pi * gFredericFreq;
    const real tL = tCircuit->component( tInductor )->get_value();
    const real tR = tCircuit->component( tResistor )->get_value();
    const real tC = tCircuit->component( tCapacitor )->get_value();
    const real tPreAmp = 1.0 / std::abs( tR + tJ * tOmega * tL );
    const cplx tZc = 1.0 / ( tJ * tOmega * tC );
    const real tPostAmp = 1.0 / std::abs( tJ * tOmega * tL
                                          + tR * tZc / ( tR + tZc ) );

    Vector< real > tCurrent( tNumSteps, 0.0 );   // source branch current
    for ( uint k = 0; k < tNumSteps; ++k )
    {
        ASSERT_TRUE( take_step( *tCircuit, tDeltaTime ) ) << "step " << k + 1;
        tCurrent( k ) = tCircuit->get_current_on_component( tSource );

        // step 100 ends at 5.00 ms, the source's positive peak
        if ( k + 1 == 100 )
        {
            EXPECT_GT( tCircuit->component( tSource )->get_value(), 0.99 );
        }

        // the switch fires at the end of step 500: its branch current is
        // pinned to zero while open and nonzero from the first closed step
        if ( k + 1 == 499 )
        {
            EXPECT_LT( std::abs( tCircuit->get_current_on_component( tSwitch ) ), 1e-12 );
        }
        else if ( k + 1 == 500 )
        {
            EXPECT_GT( std::abs( tCircuit->get_current_on_component( tSwitch ) ), 1e-6 );
        }
    }

    // sign, once, at a known-nonzero sample: at step 100 the RL current
    // lags the source peak by only 3.6 degrees, so the delivered current
    // is ~ +0.996 A and the source branch current is its negative
    EXPECT_LT( tCurrent( 99 ), -0.9 );

    // last pre-switch period: steps 100 ... 499
    real tPeakPre = 0.0;
    for ( uint k = 99; k < 499; ++k )
    {
        tPeakPre = std::max( tPeakPre, std::abs( tCurrent( k ) ) );
    }
    EXPECT_NEAR( tPeakPre, tPreAmp, 0.01 * tPreAmp );

    // last post-switch period: steps 1201 ... 1600 ( ~11 to 17 RC after
    // the switch, so the transient has settled )
    real tPeakPost = 0.0;
    for ( uint k = 1200; k < tNumSteps; ++k )
    {
        tPeakPost = std::max( tPeakPost, std::abs( tCurrent( k ) ) );
    }
    EXPECT_NEAR( tPeakPost, tPostAmp, 0.01 * tPostAmp );
}

//-----------------------------------------------------------------------------
// diode bridge with the load the demo lacked: 10 V / 60 Hz across 0-1,
// four Shockley diodes, 100 ohm from node 2 to ground ( node 3 ). With
// v(0) > v(1) the path is 0 -> D1 -> 2 -> R -> 3 -> D2 -> 1; with
// v(1) > v(0) it is 1 -> D4 -> 2 -> R -> 3 -> D3 -> 0. Both deliver + to
// node 2, which is what makes it full-wave.
//
// This is the first transient exercise of the diode Newton. Convergence
// is driven with omega = 1 ( no relaxation ); the worst-case iteration
// count is reported as the evidence for or against adding one. IEEE
// overflow of exp( dV / Vt ) needs dV > 18.4 V, which a 10 V source
// cannot reach.
//
// The source carries a half-step phase offset so that its zero crossings
// fall BETWEEN samples. This circuit has no reactive element, so at a
// sample where the source is ~1e-13 V the exact solution is ~1e-14 in
// every unknown, and ElectricalCircuit::residual() -- norm( rhs ) over
// norm( J * x ) -- divides roundoff by roundoff: the Newton plateaus at
// ~1e-7 forever although the solution is converged to machine precision.
// Measured 2026-09-04 with phase 0: steps 100, 200, 300 and 400 fail,
// every other step converges; a halved timestep fails the same way. A
// commutation ( conducting pair swap ) costs about 9 iterations here.
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, DiodeBridgeFullWaveRectifies )
{
    const real tFreq = 60.0;
    const real tAmplitude = 10.0;
    const real tIs = 0.1e-3;
    const real tVt = 0.026;
    const real tR = 100.0;
    const uint tStepsPerPeriod = 200;
    const real tDeltaTime = 1.0 / ( tFreq * tStepsPerPeriod );
    const uint tNumSteps = 2 * tStepsPerPeriod;

    // half a step of phase: zero crossings sit midway between samples
    const real tPhase = belfem::constant::pi / tStepsPerPeriod;

    ElectricalCircuit tCircuit( 4 );       // ground = node 3
    SourceFunction * tFun = new SourceFunction();
    tFun->set_periodic( SourceFunctionType::Sine, tAmplitude, 1.0 / tFreq, tPhase );
    tCircuit.create_voltage_source( tFun, 0, 1, "V1" );  // house order: V first
    tCircuit.create_diode( tIs, tVt, 0, 2, "D1" );
    tCircuit.create_diode( tIs, tVt, 3, 1, "D2" );
    tCircuit.create_diode( tIs, tVt, 3, 0, "D3" );
    tCircuit.create_diode( tIs, tVt, 1, 2, "D4" );
    tCircuit.create_resistor( tR, 2, 3, "R1" );

    Vector< real > tOut( tNumSteps, 0.0 );  // v(2) - v(3)
    uint tMaxIterations = 0;
    for ( uint k = 0; k < tNumSteps; ++k )
    {
        uint tIterations = 0;
        ASSERT_TRUE( take_step( tCircuit, tDeltaTime, &tIterations ) )
                << "step " << k + 1;
        tMaxIterations = std::max( tMaxIterations, tIterations );
        tOut( k ) = tCircuit.get_voltage_on_node( 2 )
                  - tCircuit.get_voltage_on_node( 3 );

        // rectification: the floor is 1e-6, not tighter, because the
        // Newton tolerance is a relative residual that bounds no node
        // voltage to 1e-9 V at the commutation instants
        EXPECT_GE( tOut( k ), -1e-6 ) << "step " << k + 1;
    }
    RecordProperty( "max_newton_iterations", static_cast< int >( tMaxIterations ) );
    std::printf( "diode bridge: worst-case Newton iterations per step = %u\n",
                 tMaxIterations );

    // expected peak: two diodes in series each drop Vt * ln( I / Is + 1 ),
    // so the peak output is the fixed point of
    //   v = V0 - 2 Vt ln( v / ( R Is ) + 1 )        ( 9.6426 V )
    real tPeakWant = tAmplitude;
    for ( uint i = 0; i < 100; ++i )
    {
        tPeakWant = tAmplitude
                  - 2.0 * tVt * std::log( tPeakWant / ( tR * tIs ) + 1.0 );
    }

    // second period only: peak and number of maxima
    real tPeak = 0.0;
    for ( uint k = tStepsPerPeriod; k < tNumSteps; ++k )
    {
        tPeak = std::max( tPeak, tOut( k ) );
    }
    EXPECT_NEAR( tPeak, tPeakWant, 0.01 * tPeakWant );

    // full-wave: two maxima per source period; the threshold at half the
    // observed peak keeps numerical wobble near zero from manufacturing
    // extra maxima ( a half-wave result has one ). With the half-step
    // phase each true peak sits midway between two samples that are equal
    // up to roundoff; the strict/non-strict pair below counts such a tie
    // exactly once whichever side roundoff favors
    uint tNumMaxima = 0;
    for ( uint k = tStepsPerPeriod; k + 1 < tNumSteps; ++k )
    {
        if ( tOut( k ) > tOut( k - 1 ) && tOut( k ) >= tOut( k + 1 )
             && tOut( k ) > 0.5 * tPeak )
        {
            ++tNumMaxima;
        }
    }
    EXPECT_EQ( tNumMaxima, 2u );
}

//-----------------------------------------------------------------------------
// netlist equivalence: Frederic's circuit as a .cir deck built through
// NgspiceCircuitFactory must reproduce the hand-built run step by step,
// not just at the end, or a transient divergence that re-converges would
// pass. The deck's SPICE numbers evaluate to the same arithmetic the
// builder writes, so the gate is 1e-12 with an absolute floor. If it
// misses, do not loosen it: check ( 1 ) that the order directives took
// effect ( the factory default is BDF1 ), ( 2 ) the arithmetic forms,
// ( 3 ) whether the switch fired on the same step.
//-----------------------------------------------------------------------------

TEST( ElectricalCircuit, NetlistTwinMatchesFredericCircuit )
{
    const real tDeltaTime = 50e-6;
    const uint tNumSteps = 1600;
    const uint tNumComponents = 5;
    const uint tSource = 0;
    const uint tSwitch = 3;

    // node names chosen so first-appearance packing reproduces the
    // hand-built indices; resolved through node_index() below anyway
    const std::string tDeck =
        "* Frederic's circuit, netlist twin of create_frederic_circuit()\n"
        "V1 n0 0 SIN(0 1 50)\n"
        "L1 n0 n1 200u\n"
        "R1 n1 0 1\n"
        "* belfem: switch S1 n+=n1 n-=n2 state=open t_switch=24.975m\n"
        "C1 n2 0 3.183m\n"
        "* belfem: order L1 3\n"
        "* belfem: order C1 3\n"
        ".end\n";
    const std::string tPath = write_scratch( "frederic_twin.cir", tDeck );

    std::unique_ptr< ElectricalCircuit > tHand = create_frederic_circuit();

    NgspiceCircuitFactory tFactory( tPath );
    ASSERT_EQ( tFactory.node_index( "n0" ), index_t( 0 ) );
    ASSERT_EQ( tFactory.node_index( "n1" ), index_t( 1 ) );
    ASSERT_EQ( tFactory.node_index( "n2" ), index_t( 2 ) );
    ASSERT_EQ( tFactory.node_index( "0" ),  index_t( 3 ) );
    std::unique_ptr< ElectricalCircuit > tTwin( tFactory.circuit() );

    // components are created in merged line order, so index 0 is the
    // voltage source on both sides; the parser folds labels to lowercase
    ASSERT_EQ( tHand->component( tSource )->get_label(), "v1" );
    ASSERT_EQ( tTwin->component( tSource )->get_label(), "v1" );

    real tMaxError = 0.0;                  // scaled by the mixed tolerance
    uint tFirstBadStep = 0;
    for ( uint k = 0; k < tNumSteps; ++k )
    {
        ASSERT_TRUE( take_step( *tHand, tDeltaTime ) ) << "hand-built, step " << k + 1;
        ASSERT_TRUE( take_step( *tTwin, tDeltaTime ) ) << "netlist, step " << k + 1;

        real tStepError = 0.0;
        for ( index_t n = 0; n < 3; ++n )
        {
            const real tWant = tHand->get_voltage_on_node( n );
            tStepError = std::max( tStepError,
                    std::abs( tTwin->get_voltage_on_node( n ) - tWant )
                    / mixed_tol( tWant, 1.0 ) );
        }
        for ( index_t c = 0; c < tNumComponents; ++c )
        {
            const real tWant = tHand->get_current_on_component( c );
            tStepError = std::max( tStepError,
                    std::abs( tTwin->get_current_on_component( c ) - tWant )
                    / mixed_tol( tWant, 1.0 ) );
        }
        if ( tStepError > 1e-12 && tFirstBadStep == 0 )
        {
            tFirstBadStep = k + 1;
        }
        tMaxError = std::max( tMaxError, tStepError );

        // both runs must fire the switch at the end of step 500
        if ( k + 1 == 499 )
        {
            EXPECT_LT( std::abs( tHand->get_current_on_component( tSwitch ) ), 1e-12 );
            EXPECT_LT( std::abs( tTwin->get_current_on_component( tSwitch ) ), 1e-12 );
        }
        else if ( k + 1 == 500 )
        {
            EXPECT_GT( std::abs( tHand->get_current_on_component( tSwitch ) ), 1e-6 );
            EXPECT_GT( std::abs( tTwin->get_current_on_component( tSwitch ) ), 1e-6 );
        }
    }
    EXPECT_LE( tMaxError, 1e-12 )
            << "netlist twin separates from the hand-built run, first at step "
            << tFirstBadStep << " ( max scaled error " << tMaxError << " )";
}
