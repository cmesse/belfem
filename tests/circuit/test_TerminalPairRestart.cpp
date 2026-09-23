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
 * A restart restores the FEM terminal pair. The pair is not an
 * unknown-current dof, so update_components() cannot rebuild its current
 * from x: FEMTwoTerminals::load_state must set the active current from the
 * restored mIn itself, or the first Jacobian after the restart reads the
 * current of a cold circuit. Asserted on get_current(), not on the
 * round-tripped dataset, which would pass without the restore.
 *
 * The three history registers ( h, i, v ) and the pair ( in, vn ) must also
 * survive a load followed by a save unchanged.
 *
 * Scratch files use absolute paths under the test temp directory: Ascii and
 * HDF5 would otherwise resolve against a working directory ctest changes.
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <memory>
#include <string>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_ElectricalCircuit.hpp"
#include "cl_FEMTwoTerminals.hpp"

#ifdef BELFEM_HDF5
#include "cl_HDF5.hpp"
#include "hdf5_tools.hpp"

using belfem::real ;
using belfem::uint ;
using belfem::index_t ;
using belfem::Vector ;
using belfem::electronics::ElectricalCircuit ;
using belfem::electronics::FEMTwoTerminals ;

namespace
{
    const real gDeltaTime = 1.0e-3 ;

    // the component index of the terminal pair in make_circuit(); the
    // restart datasets carry it as their prefix
    const std::string gPairPrefix = "c001_" ;

    // an absolute scratch path, removed before use and on scope exit
    class ScratchFile
    {
        const std::string mPath ;

    public:

        explicit ScratchFile( const std::string & aName ) :
            mPath( ::testing::TempDir() + aName )
        {
            std::remove( mPath.c_str() );
        }

        ScratchFile( const ScratchFile & ) = delete ;
        ScratchFile( ScratchFile && ) = delete ;
        ScratchFile & operator=( const ScratchFile & ) = delete ;
        ScratchFile & operator=( ScratchFile && ) = delete ;

        ~ScratchFile()
        {
            std::remove( mPath.c_str() );
        }

        const std::string &
        path() const
        {
            return mPath ;
        }
    };

    // a 1 kOhm resistor fed by one FEM terminal pair between node 0 and
    // ground, which is the last node ( 1 ) of a two-node circuit
    std::unique_ptr< ElectricalCircuit >
    make_circuit()
    {
        std::unique_ptr< ElectricalCircuit > tCircuit( new ElectricalCircuit( 2 ) );
        tCircuit->create_resistor( 1000.0, 0, 1, "R1" );
        tCircuit->create_terminal_pair( 0, 1, "TP" );
        return tCircuit ;
    }

    FEMTwoTerminals *
    pair_of( ElectricalCircuit & aCircuit )
    {
        return aCircuit.terminal_pairs()( 0 );
    }

    // one accepted step in controller order
    bool
    take_step( ElectricalCircuit & aCircuit )
    {
        aCircuit.set_timestep( gDeltaTime );
        aCircuit.shift();
        aCircuit.compute_MNA_matrix();
        aCircuit.set_omega( 1.0 );

        real tEpsilon = BELFEM_REAL_MAX ;
        for ( uint tIt = 0; tIt < 50 && tEpsilon > 1.0e-9; ++tIt )
        {
            aCircuit.compute_jacobian_and_rhs();
            aCircuit.solve();
            tEpsilon = aCircuit.residual();
        }
        return tEpsilon <= 1.0e-9 ;
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

    // a dump whose pair holds I = 12.5 A, V = 0.4 V, and whose registers
    // hold the pushed step before it ( I = 3 A, V = 0.1 V, dt )
    void
    write_dump( const std::string & aPath )
    {
        std::unique_ptr< ElectricalCircuit > tSource = make_circuit();
        pair_of( *tSource )->set_IV( 3.0, 0.1 );
        ASSERT_TRUE( take_step( *tSource ) );
        pair_of( *tSource )->set_IV( 12.5, 0.4 );
        save_circuit( *tSource, aPath );
    }

    void
    read_register( const std::string & aPath,
                   const std::string & aLabel,
                   Vector< real >    & aSamples )
    {
        belfem::HDF5 tFile( aPath, belfem::FileMode::OPEN_RDONLY );
        hid_t tGroup = tFile.select_group( "circuit" );
        herr_t tStatus = 0 ;
        belfem::hdf5::load_vector_from_file( tGroup, gPairPrefix + aLabel, aSamples, tStatus );
    }

    real
    read_scalar( const std::string & aPath, const std::string & aLabel )
    {
        belfem::HDF5 tFile( aPath, belfem::FileMode::OPEN_RDONLY );
        hid_t tGroup = tFile.select_group( "circuit" );
        herr_t tStatus = 0 ;
        real tValue = BELFEM_QUIET_NAN ;
        belfem::hdf5::load_scalar_from_file( tGroup, gPairPrefix + aLabel, tValue, tStatus );
        return tValue ;
    }
}

//------------------------------------------------------------------------------

TEST( ElectricalCircuit, TerminalPairRestoresItsActiveCurrent )
{
    ScratchFile tDump( "test_terminal_pair_restore.hdf5" );
    write_dump( tDump.path() );
    ASSERT_FALSE( HasFatalFailure() );

    std::unique_ptr< ElectricalCircuit > tRestart = make_circuit();
    ASSERT_EQ( pair_of( *tRestart )->get_current(), 0.0 );

    load_circuit( *tRestart, tDump.path() );

    // the ACTIVE current, the one the Jacobian reads
    EXPECT_EQ( pair_of( *tRestart )->get_current(), 12.5 );
}

//------------------------------------------------------------------------------

TEST( ElectricalCircuit, TerminalPairStateSurvivesLoadAndSave )
{
    ScratchFile tDumpA( "test_terminal_pair_state_a.hdf5" );
    ScratchFile tDumpB( "test_terminal_pair_state_b.hdf5" );
    write_dump( tDumpA.path() );
    ASSERT_FALSE( HasFatalFailure() );

    {
        std::unique_ptr< ElectricalCircuit > tRestart = make_circuit();
        load_circuit( *tRestart, tDumpA.path() );
        save_circuit( *tRestart, tDumpB.path() );
    }

    // the first dump holds what write_dump put there; otherwise the
    // comparison below could agree on two empty registers
    Vector< real > tA ;
    read_register( tDumpA.path(), "i", tA );
    ASSERT_EQ( tA.length(), 1u );
    EXPECT_EQ( tA( 0 ), 3.0 );
    read_register( tDumpA.path(), "v", tA );
    ASSERT_EQ( tA.length(), 1u );
    EXPECT_EQ( tA( 0 ), 0.1 );
    read_register( tDumpA.path(), "h", tA );
    ASSERT_EQ( tA.length(), 1u );
    EXPECT_EQ( tA( 0 ), gDeltaTime );

    // pure serialization round trip: bitwise
    for ( const std::string tLabel : { "h", "i", "v" } )
    {
        Vector< real > tB ;
        read_register( tDumpA.path(), tLabel, tA );
        read_register( tDumpB.path(), tLabel, tB );
        ASSERT_EQ( tB.length(), tA.length() ) << "register " << tLabel ;
        for ( index_t k = 0; k < tA.length(); ++k )
        {
            EXPECT_EQ( tB( k ), tA( k ) ) << "register " << tLabel << " sample " << k ;
        }
    }

    EXPECT_EQ( read_scalar( tDumpB.path(), "in" ), 12.5 );
    EXPECT_EQ( read_scalar( tDumpB.path(), "vn" ), 0.4 );
}

#endif
