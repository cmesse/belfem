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
 * Unit tests for the hybrid circuit path ( plan decision 6b ):
 * circuit { file : <netlist> } in the deck, lumped topology from the
 * netlist, terminal pairs and output from the deck with node references
 * resolved as netlist node NAMES ( O5 ). Also one classic-path regression
 * that pins the read_terminal_pair/read_output refactor as
 * behavior-preserving.
 *
 * All scratch files use absolute paths: Ascii resolves relative paths
 * against getenv("PWD"), which is stale under ctest.
 */

#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_ElectricalCircuitFactory.hpp"
#include "cl_ElectricalCircuit.hpp"
#include "Component_Enums.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"

using namespace belfem;
using electronics::ElectricalCircuitFactory;
using electronics::ElectricalCircuit;
using electronics::ComponentType;

namespace
{
    // write aContent to an absolute scratch path and return it
    string
    write_scratch( const string & aName, const string & aContent )
    {
        const string tPath = ::testing::TempDir() + aName;
        std::ofstream tFile( tPath );
        tFile << aContent;
        return tPath;
    }

    // drive the circuit through one converged solve, controller-style
    real
    solve_and_probe( ElectricalCircuit * aCircuit, const index_t aNode )
    {
        aCircuit->set_timestep( 1.0e-3 );
        aCircuit->shift();
        aCircuit->compute_MNA_matrix();

        real tEpsilon = BELFEM_REAL_MAX;
        uint tIteration = 0;
        while ( tIteration++ < 20 && tEpsilon > 1.0e-10 )
        {
            aCircuit->set_omega( 1.0 );
            aCircuit->compute_jacobian_and_rhs();
            aCircuit->solve();
            tEpsilon = aCircuit->residual();
        }
        return aCircuit->get_voltage_on_node( aNode );
    }
}

//------------------------------------------------------------------------------
// the hybrid path
//------------------------------------------------------------------------------

TEST( HybridFactory, NetlistDeckBuildsAndSolves )
{
    const string tCir = write_scratch( "hybrid_divider.cir",
        "divider deck\n"
        "V1 in 0 DC 10\n"
        "R1 in mid 1k\n"
        "R2 mid 0 1k\n"
        ".tran 0.1m 10m\n"
        ".end\n" );

    const string tConf = write_scratch( "hybrid_divider.conf",
        "circuit\n"
        "{\n"
        "    file : " + tCir + " ;\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    ElectricalCircuitFactory tFactory( tConf, tBCs );

    ElectricalCircuit * tCircuit =
            dynamic_cast< ElectricalCircuit * >( tFactory.circuit() );
    ASSERT_NE( tCircuit, nullptr );

    // in, mid pack first; ground is last
    EXPECT_EQ( tCircuit->number_of_nodes(), 3u );
    EXPECT_EQ( tCircuit->number_of_components(), 3u );

    // v(mid) = 5 V through the real solve
    EXPECT_NEAR( solve_and_probe( tCircuit, 1 ), 5.0, 1.0e-9 );

    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
}

TEST( HybridFactory, TerminalPairByNodeName )
{
    const string tCir = write_scratch( "hybrid_tp.cir",
        "coupled deck\n"
        "Is tape 0 SIN(0 500 10)\n"
        "L1 tape 0 2.5u\n"
        "* belfem: order L1 2\n"
        ".end\n" );

    const string tConf = write_scratch( "hybrid_tp.conf",
        "circuit\n"
        "{\n"
        "    file : " + tCir + " ;\n"
        "    topology\n"
        "    {\n"
        "        terminal pair\n"
        "        {\n"
        "            label : Z1 ;\n"
        "            node + : tape ;\n"
        "            node - : gnd ;\n"
        "            input curves : [1,2] ;\n"
        "            output curves : [3,4] ;\n"
        "        }\n"
        "    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    ElectricalCircuitFactory tFactory( tConf, tBCs );

    ElectricalCircuit * tCircuit =
            dynamic_cast< ElectricalCircuit * >( tFactory.circuit() );
    ASSERT_NE( tCircuit, nullptr );

    // Is, L1 from the netlist, then the deck's terminal pair
    ASSERT_EQ( tCircuit->number_of_components(), 3u );
    EXPECT_EQ( tCircuit->component( 2 )->component_type(),
               ComponentType::TERMINALPAIR );

    // hybrid mode is case-insensitive: the deck's Z1 folds like the
    // netlist instance names, so the currents lookup cannot miss
    EXPECT_EQ( tCircuit->component( 2 )->get_label(), "z1" );

    // the name-resolved pair nodes: tape packs first, ground is last
    EXPECT_EQ( tCircuit->component( 2 )->get_node_plus()->index(),  0u );
    EXPECT_EQ( tCircuit->component( 2 )->get_node_minus()->index(), 1u );

    // the terminal pair spawned exactly one boundary condition
    ASSERT_EQ( tBCs.size(), 1u );
    EXPECT_EQ( tBCs( 0 )->label(), "z1" );

    for ( auto * tBC : tBCs ) { delete tBC; }
    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
}

TEST( HybridFactory, OutputLabelsSurviveSaveTimestep )
{
    // the regression the audit round asked for: a mixed-case currents list
    // naming a netlist instance AND a deck terminal pair must survive
    // save_timestep's label lookup ( a one-sided fold failed here )
    const string tCir = write_scratch( "hybrid_out.cir",
        "output deck\n"
        "V1 in 0 DC 10\n"
        "R1 in mid 1k\n"
        "R2 mid 0 1k\n"
        ".end\n" );

    const string tOut = ::testing::TempDir() + "hybrid_out_results.txt";
    const string tConf = write_scratch( "hybrid_out.conf",
        "circuit\n{\n"
        "    file : " + tCir + " ;\n"
        "    output\n    {\n"
        "        file : " + tOut + " ;\n"
        "        currents : V1, R1 ;\n"
        "        voltages : mid, 0 ;\n"
        "    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    ElectricalCircuitFactory tFactory( tConf, tBCs );
    ElectricalCircuit * tCircuit =
            dynamic_cast< ElectricalCircuit * >( tFactory.circuit() );
    ASSERT_NE( tCircuit, nullptr );

    EXPECT_NEAR( solve_and_probe( tCircuit, 1 ), 5.0, 1.0e-9 );
    tCircuit->save_timestep();

    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
    std::remove( tOut.c_str() );
}

TEST( HybridFactory, RejectsCrossSourceDuplicateLabel )
{
    // a deck label colliding with a netlist instance ( case-only included )
    const string tCir = write_scratch( "hybrid_dup.cir",
        "t\nR1 a 0 1k\n.end\n" );
    const string tConf = write_scratch( "hybrid_dup.conf",
        "circuit\n{\n"
        "    file : " + tCir + " ;\n"
        "    topology\n    {\n"
        "        terminal pair\n        {\n"
        "            label : R1 ;\n"
        "            node + : a ;\n"
        "            node - : 0 ;\n"
        "            input curves : [1] ;\n"
        "        }\n    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_THROW( ElectricalCircuitFactory tF( tConf, tBCs ),
                  std::runtime_error );
    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
}

//------------------------------------------------------------------------------
// hybrid refusals
//------------------------------------------------------------------------------

TEST( HybridFactory, RejectsNumberOfNodesWithFile )
{
    const string tCir = write_scratch( "hybrid_nn.cir",
        "t\nR1 a 0 1k\n.end\n" );
    const string tConf = write_scratch( "hybrid_nn.conf",
        "circuit\n{\n"
        "    file : " + tCir + " ;\n"
        "    number of nodes : 2 ;\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_THROW( ElectricalCircuitFactory tF( tConf, tBCs ),
                  std::runtime_error );
    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
}

TEST( HybridFactory, RejectsLumpedElementsInDeckTopology )
{
    const string tCir = write_scratch( "hybrid_lump.cir",
        "t\nR1 a 0 1k\n.end\n" );
    const string tConf = write_scratch( "hybrid_lump.conf",
        "circuit\n{\n"
        "    file : " + tCir + " ;\n"
        "    topology\n    {\n"
        "        resistor\n        {\n"
        "            node + : a ;\n"
        "            node - : 0 ;\n"
        "            value : 5 Ohm ;\n"
        "        }\n    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_THROW( ElectricalCircuitFactory tF( tConf, tBCs ),
                  std::runtime_error );
    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
}

TEST( HybridFactory, RejectsUnknownNodeName )
{
    const string tCir = write_scratch( "hybrid_unknown.cir",
        "t\nR1 a 0 1k\n.end\n" );
    const string tConf = write_scratch( "hybrid_unknown.conf",
        "circuit\n{\n"
        "    file : " + tCir + " ;\n"
        "    topology\n    {\n"
        "        terminal pair\n        {\n"
        "            node + : nowhere ;\n"
        "            node - : 0 ;\n"
        "            input curves : [1] ;\n"
        "        }\n    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_THROW( ElectricalCircuitFactory tF( tConf, tBCs ),
                  std::runtime_error );
    std::remove( tCir.c_str() );
    std::remove( tConf.c_str() );
}

//------------------------------------------------------------------------------
// classic-path regression: the refactor is behavior-preserving
//------------------------------------------------------------------------------

TEST( HybridFactory, ClassicDeckStillBuildsAndSolves )
{
    const string tConf = write_scratch( "classic_divider.conf",
        "circuit\n"
        "{\n"
        "    number of nodes : 3 ;\n"
        "    topology\n"
        "    {\n"
        "        voltage source\n"
        "        {\n"
        "            label : V1 ;\n"
        "            node + : 0 ;\n"
        "            node - : 2 ;\n"
        "            type : constant ;\n"
        "            amplitude : 10 V ;\n"
        "        }\n"
        "        resistor\n"
        "        {\n"
        "            label : R1 ;\n"
        "            node + : 0 ;\n"
        "            node - : 1 ;\n"
        "            value : 1000 Ohm ;\n"
        "        }\n"
        "        resistor\n"
        "        {\n"
        "            label : R2 ;\n"
        "            node + : 1 ;\n"
        "            node - : 2 ;\n"
        "            value : 1000 Ohm ;\n"
        "        }\n"
        "    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    ElectricalCircuitFactory tFactory( tConf, tBCs );

    ElectricalCircuit * tCircuit =
            dynamic_cast< ElectricalCircuit * >( tFactory.circuit() );
    ASSERT_NE( tCircuit, nullptr );
    EXPECT_EQ( tCircuit->number_of_components(), 3u );

    // classic indices: node 1 is the divider midpoint, node 2 is ground
    EXPECT_NEAR( solve_and_probe( tCircuit, 1 ), 5.0, 1.0e-9 );

    std::remove( tConf.c_str() );
}

//------------------------------------------------------------------------------
//  'expk' is refused, not ignored
//------------------------------------------------------------------------------

// The key was parsed and stored for years and never reached the waveform.
// Dropping the read in silence would have left exactly the class this
// check exists to close, so it is a hard error now.
TEST( HybridFactory, ExpkOnSigmoidIsRefused )
{
    const string tConf = write_scratch( "expk_sigmoid.conf",
        "circuit\n"
        "{\n"
        "    number of nodes : 2 ;\n"
        "    topology\n"
        "    {\n"
        "        voltage source\n"
        "        {\n"
        "            label : Vs ;\n"
        "            node + : 0 ;\n"
        "            node - : 1 ;\n"
        "            type : sigmoid ;\n"
        "            amplitude : 1 V ;\n"
        "            period : 1 s ;\n"
        "            offset : 0 s ;\n"
        "            expk : 99 ;\n"
        "        }\n"
        "        resistor { label : R1 ; node + : 0 ; node - : 1 ; value : 1 Ohm ; }\n"
        "    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_THROW( ElectricalCircuitFactory tF( tConf, tBCs ), std::runtime_error );
}

//------------------------------------------------------------------------------

// The refusal sits BEFORE the type dispatch on purpose. 'expk' was only ever
// read on the sigmoid branch, so a sigmoid-only check would have left it
// silently ignored on every other source shape -- the same defect, moved.
TEST( HybridFactory, ExpkOnNonSigmoidIsAlsoRefused )
{
    const string tConf = write_scratch( "expk_sine.conf",
        "circuit\n"
        "{\n"
        "    number of nodes : 2 ;\n"
        "    topology\n"
        "    {\n"
        "        voltage source\n"
        "        {\n"
        "            label : Vs ;\n"
        "            node + : 0 ;\n"
        "            node - : 1 ;\n"
        "            type : sine ;\n"
        "            amplitude : 1 V ;\n"
        "            frequency : 50 Hz ;\n"
        "            expk : 99 ;\n"
        "        }\n"
        "        resistor { label : R1 ; node + : 0 ; node - : 1 ; value : 1 Ohm ; }\n"
        "    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_THROW( ElectricalCircuitFactory tF( tConf, tBCs ), std::runtime_error );
}

//------------------------------------------------------------------------------

// The guard must not fire on a deck that never mentions the key.
TEST( HybridFactory, SigmoidWithoutExpkStillBuilds )
{
    const string tConf = write_scratch( "sigmoid_ok.conf",
        "circuit\n"
        "{\n"
        "    number of nodes : 2 ;\n"
        "    topology\n"
        "    {\n"
        "        voltage source\n"
        "        {\n"
        "            label : Vs ;\n"
        "            node + : 0 ;\n"
        "            node - : 1 ;\n"
        "            type : sigmoid ;\n"
        "            amplitude : 1 V ;\n"
        "            period : 1 s ;\n"
        "            offset : 0 s ;\n"
        "            fuzzyness : 1e-4 ;\n"
        "        }\n"
        "        resistor { label : R1 ; node + : 0 ; node - : 1 ; value : 1 Ohm ; }\n"
        "    }\n"
        "}\n" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs;
    EXPECT_NO_THROW( ElectricalCircuitFactory tF( tConf, tBCs ) );
}
