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
 * Adjacency-graph regression tests for ElectricalCircuit::compute_adjacency().
 *
 * WHAT THIS IS NOT: it is not a validation of the heap overflow these cases
 * were written for. That defect under-sized the vertex containers while
 * inserting the SAME set, and the degrees below are INSERTION counts, so a
 * surviving process reports them unchanged either way. That is why the old
 * suite passed green with the overflow live, and why no count test can be the
 * gate. The discriminating gate is valgrind reporting zero invalid
 * reads/writes in the circuit suite, and that is not a gtest.
 *
 * Do not read that as "all of these would have passed pre-fix". The shortfall
 * was not a uniform one slot -- two sources on one node were short by two, a
 * fully grounded source was not short at all, and a fully grounded switch
 * inserted its diagonal into a null buffer. The old suite simply did not
 * contain these cases.
 *
 * WHAT THIS IS: a pin on the adjacency graph's degree sequence. The sizing and
 * the population inside compute_adjacency() were once two copies of a filter
 * that drifted apart; the in-function assertion compares those two against
 * each other and therefore cannot see a drift that moves both. These expected
 * vectors are independent of both.
 *
 * What that does and does not catch, measured by mutating the implementation:
 *   - a drift that changes what gets INSERTED ( dropping a type from the
 *     allowlist ) fails here in any build
 *   - a drift that changes only the SIZING is INVISIBLE here, because the
 *     degrees are insertion counts. Such a drift under-allocates and is caught
 *     only by the post-fill assertion in compute_adjacency(), which exists just
 *     in an assertion-enabled build. In a release tree it can pass green
 * So this suite must be run in a DEBUG build to guard the second class -- the
 * class the original overflow belonged to.
 *
 * The degree sequence is more discriminating than the matrix nnz, which is
 * merely the SUM of these degrees: the divider [3,2,1] and the live-to-live
 * source [2,2,2] are different graphs that both sum to 6. It does not prove
 * graph identity, though -- the divider and the parallel-passive case share a
 * degree sequence AND a graph, and differ only in the construction path that
 * produced it.
 *
 * Every expected vector below was MEASURED against the implementation rather
 * than derived by hand. Both auditors independently re-derived the full set
 * from the insert rules; the degenerate and type-coverage cases were their
 * recommendations, measured afterwards.
 *
 * The graph is mVertices: the non-ground nodes in index order, then one branch
 * vertex per unknown-current component in creation order. Ground is the last
 * NODE and is not in the graph. Node degree is
 *   1 (self) + |unique allowlisted neighbors| + (unknown-current attachments)
 * and branch degree is (non-ground terminal count) + (1 if SWITCH).
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_ElectricalCircuit.hpp"
#include "cl_SourceFunction.hpp"

using belfem::real;
using belfem::uint;
using belfem::index_t;
using belfem::Cell;
using belfem::SourceFunction;
using belfem::electronics::ElectricalCircuit;

namespace
{
//-----------------------------------------------------------------------------

    // a DC source. The circuit takes ownership of the function
    SourceFunction *
    dc_source( const real aValue )
    {
        SourceFunction * tFunction = new SourceFunction();
        tFunction->set_constant( aValue );
        return tFunction;
    }

//-----------------------------------------------------------------------------

    // build the adjacency and compare the whole degree sequence. Only
    // compute_adjacency() is called -- never compute_MNA_matrix() or solve().
    // Some topologies below are structurally accepted by the creators but are
    // NOT electrically valid circuits ( a shorted source, a source with both
    // terminals on ground ); building their adjacency is defined, solving them
    // is not, so nothing here solves
    void
    expect_degrees( ElectricalCircuit & aCircuit,
                    const Cell< uint > & aExpected )
    {
        aCircuit.compute_adjacency();

        // length first: degree_of_graph_vertex() is only bounds-checked in an
        // assertion-enabled build, and the release tree has none
        ASSERT_EQ( aCircuit.number_of_graph_vertices(), aExpected.size() );

        for ( uint k = 0; k < aExpected.size(); ++k )
        {
            EXPECT_EQ( aCircuit.degree_of_graph_vertex( k ), aExpected( k ) )
                    << "degree mismatch at graph vertex " << k;
        }
    }

//-----------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------
// The topology that overflowed: the source's other terminal is GROUND, so it
// contributed nothing to the old size yet still inserted its branch vertex
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, VoltageSourceToGround )
{
    ElectricalCircuit tCircuit( 3 );                 // in(0) mid(1) gnd(2)
    tCircuit.create_voltage_source( dc_source( 10.0 ), 0, 2, "V1" );
    tCircuit.create_resistor( 1000.0, 0, 1, "R1" );
    tCircuit.create_resistor( 1000.0, 1, 2, "R2" );

    expect_degrees( tCircuit, { 3, 2, 1 } );
}

//-----------------------------------------------------------------------------
// The case that did NOT overflow: with both terminals live the old count and
// the inserts happened to balance one for one
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, VoltageSourceBetweenLiveNodes )
{
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 1, "V1" );
    tCircuit.create_resistor( 1.0, 0, 2, "R1" );
    tCircuit.create_resistor( 1.0, 1, 2, "R2" );

    // a source is NOT an allowlisted neighbor, so b is not adjacent to a
    expect_degrees( tCircuit, { 2, 2, 2 } );
}

//-----------------------------------------------------------------------------
// The switch diagonal: a switch inserts itself into its own branch vertex, and
// that slot has to be reserved whether the switch is open or closed
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, SwitchDiagonalIsReserved )
{
    // creation order fixes the graph layout: [ a, b, S1, V1 ]
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_switch( true, 0.0, 0, 1, "S1" );
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 2, "V1" );
    tCircuit.create_resistor( 1.0, 1, 2, "R1" );

    expect_degrees( tCircuit, { 3, 2, 3, 1 } );
}

//-----------------------------------------------------------------------------
// The diagonal does not depend on the switch state. compute_adjacency() keys on
// the component TYPE, never on is_closed(), and the graph is built once -- so a
// future is_closed() guard around the diagonal insert would leave the closed
// case above green and re-break the open one
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, OpenSwitchAlsoReservesItsDiagonal )
{
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_switch( false, 0.0, 0, 1, "S1" );
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 2, "V1" );
    tCircuit.create_resistor( 1.0, 1, 2, "R1" );

    expect_degrees( tCircuit, { 3, 2, 3, 1 } );
}

//-----------------------------------------------------------------------------
// Two DISTINCT sources on one node get two distinct branch vertices
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, TwoSourcesOnOneNode )
{
    ElectricalCircuit tCircuit( 2 );                 // a(0) gnd(1)
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 1, "V1" );
    tCircuit.create_voltage_source( dc_source( 2.0 ), 0, 1, "V2" );
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );

    expect_degrees( tCircuit, { 3, 1, 1 } );
}

//-----------------------------------------------------------------------------
// The SAME branch vertex inserted twice, because both terminals of one source
// sit on one node. This is the case that would collapse if the unknown-current
// attachments were ever deduplicated -- the distinct-sources case above would
// survive such a change, so it does not pin this
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, SourceShortedAcrossOneNode )
{
    ElectricalCircuit tCircuit( 2 );                 // a(0) gnd(1)
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 0, "V1" );
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );

    expect_degrees( tCircuit, { 3, 2 } );
}

//-----------------------------------------------------------------------------
// Components in parallel name the same neighbor twice and must be deduped
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, ParallelPassivesAreDeduplicated )
{
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );
    tCircuit.create_resistor( 2.0, 0, 1, "R2" );
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 2, "V1" );

    expect_degrees( tCircuit, { 3, 2, 1 } );
}

//-----------------------------------------------------------------------------
// A node with two DISTINCT neighbors, which the parallel case cannot cover
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, ChainNodeWithTwoNeighbors )
{
    ElectricalCircuit tCircuit( 4 );                 // 0 1 2 gnd(3)
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 3, "V1" );
    tCircuit.create_inductor( 1.0e-6, 2, 0, 1, "L1" );
    tCircuit.create_capacitor( 1.0e-6, 2, 1, 2, "C1" );
    tCircuit.create_resistor( 1.0, 2, 3, "R1" );

    expect_degrees( tCircuit, { 3, 3, 2, 1 } );
}

//-----------------------------------------------------------------------------
// A switch with only one live terminal: branch degree 2, not 3
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, SwitchWithOneLiveTerminal )
{
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_switch( true, 0.0, 1, 2, "S1" );
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );
    tCircuit.create_voltage_source( dc_source( 1.0 ), 0, 2, "V1" );

    expect_degrees( tCircuit, { 3, 3, 2, 1 } );
}

//-----------------------------------------------------------------------------
// Every allowlisted passive type produces one neighbor edge. Each type sits on
// its OWN node pair on purpose: put them all in parallel and unique() collapses
// them to a single neighbor, so dropping DIODE or SUPERCONDUCTOR from the
// allowlist would leave the degrees unchanged and this test would not notice.
// Chained, every type owns an edge and removing any one of them breaks it
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, EachPassiveTypeIsIndependentlyPinned )
{
    ElectricalCircuit tCircuit( 6 );                 // 0 1 2 3 4 gnd(5)
    tCircuit.create_capacitor( 1.0e-6, 2, 0, 1, "C1" );
    tCircuit.create_inductor( 1.0e-6, 2, 1, 2, "L1" );
    tCircuit.create_diode( 1.0e-12, 0.026, 2, 3, "D1" );
    tCircuit.create_superconductor( 100.0, 20.0, 1.0e-4, 1.0, 3, 4, "SC1" );
    tCircuit.create_resistor( 1.0, 4, 5, "R1" );

    expect_degrees( tCircuit, { 2, 3, 3, 3, 2 } );
}

//-----------------------------------------------------------------------------
// A terminal pair IS an allowlisted neighbor although it receives no MNA
// stamp, so it contributes a structural zero. Removing it from the allowlist
// would change the sparsity pattern
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, TerminalPairIsANeighbor )
{
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_terminal_pair( 0, 1, "Z1" );
    tCircuit.create_resistor( 1.0, 1, 2, "R1" );

    expect_degrees( tCircuit, { 2, 2 } );
}

//-----------------------------------------------------------------------------
// A current source is neither a neighbor nor an unknown-current dof, so it
// contributes no graph entry at all
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, CurrentSourceContributesNoEdge )
{
    ElectricalCircuit tCircuit( 3 );                 // a(0) b(1) gnd(2)
    tCircuit.create_current_source( dc_source( 1.0 ), 0, 1, "I1" );
    tCircuit.create_resistor( 1.0, 0, 2, "R1" );
    tCircuit.create_resistor( 1.0, 1, 2, "R2" );

    expect_degrees( tCircuit, { 1, 1 } );
}

//-----------------------------------------------------------------------------
// Degenerate: a source with BOTH terminals on ground owns a zero-degree branch
// vertex. Its container is allocated with size zero, which is a null buffer
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, FullyGroundedSourceHasEmptyBranch )
{
    ElectricalCircuit tCircuit( 2 );                 // a(0) gnd(1)
    tCircuit.create_voltage_source( dc_source( 1.0 ), 1, 1, "V1" );
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );

    expect_degrees( tCircuit, { 1, 0 } );
}

//-----------------------------------------------------------------------------
// Degenerate: the same for a switch, whose diagonal survives grounding
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, FullyGroundedSwitchKeepsItsDiagonal )
{
    ElectricalCircuit tCircuit( 2 );                 // a(0) gnd(1)
    tCircuit.create_switch( true, 0.0, 1, 1, "S1" );
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );

    expect_degrees( tCircuit, { 1, 1 } );
}

//-----------------------------------------------------------------------------
// Degenerate: no unknown-current components at all, so the graph is the nodes
//-----------------------------------------------------------------------------
TEST( CircuitAdjacency, NoUnknownCurrents )
{
    ElectricalCircuit tCircuit( 2 );                 // a(0) gnd(1)
    tCircuit.create_resistor( 1.0, 0, 1, "R1" );

    expect_degrees( tCircuit, { 1 } );
}
