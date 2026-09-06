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
 * Unit tests for src/circuit/cl_NgspiceCircuitFactory.{hpp,cpp}: netlist IR
 * -> ElectricalCircuit under the O5/O8/O9 contracts, including one solved
 * DC divider as the end-to-end analytic gate.
 */

#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_NetlistParser.hpp"
#include "cl_NgspiceCircuitFactory.hpp"
#include "cl_ElectricalCircuit.hpp"
#include "Component_Enums.hpp"

using namespace belfem;
using electronics::NetlistParser;
using electronics::NgspiceCircuitFactory;
using electronics::ElectricalCircuit;
using electronics::ComponentType;

//------------------------------------------------------------------------------
// node map ( O9 packing, O5 lookup )
//------------------------------------------------------------------------------

TEST( NgspiceFactory, NodeMapPacksFirstAppearanceGroundLast )
{
    Cell< string > tLines = {
        "title",
        "R1 out mid 1k",
        "R2 mid 0 1k",
        "V1 in gnd 5",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );

    // out, mid, in pack 0..2 in first appearance; ground is last
    EXPECT_EQ( tFactory.number_of_nodes(), 4u );
    EXPECT_EQ( tFactory.node_index( "out" ), 0u );
    EXPECT_EQ( tFactory.node_index( "mid" ), 1u );
    EXPECT_EQ( tFactory.node_index( "in" ),  2u );
    // "0" and "gnd" are one node, on the last index
    EXPECT_EQ( tFactory.node_index( "0" ),   3u );
    EXPECT_EQ( tFactory.node_index( "gnd" ), 3u );

    EXPECT_THROW( tFactory.node_index( "nowhere" ), std::runtime_error );

    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );
    EXPECT_EQ( tCircuit->number_of_nodes(), 4u );
    EXPECT_EQ( tCircuit->number_of_components(), 3u );

    // a second handover is refused
    EXPECT_THROW( tFactory.circuit(), std::runtime_error );
}

TEST( NgspiceFactory, GroundlessDeckIsRefused )
{
    Cell< string > tLines = { "t", "R1 a b 1k", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tLines ) ) ),
                  std::runtime_error );
}

TEST( NgspiceFactory, ComponentsInNetlistOrderWithLabels )
{
    Cell< string > tLines = {
        "title",
        "Is 1 0 SIN(0 500 10)",
        "L1 1 0 2.5u",
        "* belfem: order L1 2",
        "C1 1 0 100",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );

    ASSERT_EQ( tCircuit->number_of_components(), 3u );
    EXPECT_EQ( tCircuit->component( 0 )->component_type(),
               ComponentType::CURRENTSOURCE );
    EXPECT_EQ( tCircuit->component( 0 )->get_label(), "is" );
    EXPECT_EQ( tCircuit->component( 1 )->component_type(),
               ComponentType::INDUCTOR );
    EXPECT_EQ( tCircuit->component( 2 )->component_type(),
               ComponentType::CAPACITOR );
}

TEST( NgspiceFactory, DirectiveComponentsAreCreated )
{
    Cell< string > tLines = {
        "title",
        "R1 1 0 1k",
        "* belfem: superconductor SC1 n+=1 n-=0 Ic=200 n=20 Ec=1e-4 length=1",
        "* belfem: switch S1 n+=1 n-=0 state=open t_switch=20m",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );

    ASSERT_EQ( tCircuit->number_of_components(), 3u );
    EXPECT_EQ( tCircuit->component( 1 )->component_type(),
               ComponentType::SUPERCONDUCTOR );
    EXPECT_EQ( tCircuit->component( 1 )->get_label(), "sc1" );
    EXPECT_EQ( tCircuit->component( 2 )->component_type(),
               ComponentType::SWITCH );
}

TEST( NgspiceFactory, DiodeWithModel )
{
    Cell< string > tLines = {
        "title",
        "V1 in 0 5",
        "D1 in 0 DMOD",
        ".model DMOD D( is=1e-14 n=1.5 )",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );
    ASSERT_EQ( tCircuit->number_of_components(), 2u );
    EXPECT_EQ( tCircuit->component( 1 )->component_type(),
               ComponentType::DIODE );
}

//------------------------------------------------------------------------------
// the analytic end-to-end gate: a solved DC divider
//------------------------------------------------------------------------------

TEST( NgspiceFactory, SolvedDcDivider )
{
    // V1 drives 10 V between "in" and ground; two equal resistors divide
    // it: v(mid) = 5 V. This also pins the V-source polarity convention:
    // the FIRST node of the card is the + terminal, as in SPICE.
    Cell< string > tLines = {
        "divider",
        "V1 in 0 DC 10",
        "R1 in mid 1k",
        "R2 mid 0 1k",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    const index_t tIn  = tFactory.node_index( "in" );
    const index_t tMid = tFactory.node_index( "mid" );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );

    // the controller's per-step sequence ( shift, stamp, iterate )
    tCircuit->set_timestep( 1.0e-3 );
    tCircuit->shift();
    tCircuit->compute_MNA_matrix();

    real tEpsilon = BELFEM_REAL_MAX;
    uint tIteration = 0;
    while ( tIteration++ < 20 && tEpsilon > 1.0e-10 )
    {
        tCircuit->set_omega( 1.0 );
        tCircuit->compute_jacobian_and_rhs();
        tCircuit->solve();
        tEpsilon = tCircuit->residual();
    }

    EXPECT_LT( tEpsilon, 1.0e-10 );
    EXPECT_NEAR( tCircuit->get_voltage_on_node( tIn ),  10.0, 1.0e-9 );
    EXPECT_NEAR( tCircuit->get_voltage_on_node( tMid ),  5.0, 1.0e-9 );
}

TEST( NgspiceFactory, GndOnlyDeck )
{
    // ground spelled only "gnd", never "0"
    Cell< string > tLines = {
        "title",
        "V1 a gnd 5",
        "R1 a gnd 1k",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    EXPECT_EQ( tFactory.number_of_nodes(), 2u );
    EXPECT_EQ( tFactory.node_index( "a" ),   0u );
    EXPECT_EQ( tFactory.node_index( "gnd" ), 1u );
    EXPECT_EQ( tFactory.node_index( "0" ),   1u );
}

TEST( NgspiceFactory, NodeIndexFoldsCase )
{
    // the O5 lookup folds like the parser map ( deck-side values are not
    // guaranteed lower case )
    Cell< string > tLines = { "t", "R1 in 0 1k", ".end" };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    EXPECT_EQ( tFactory.node_index( "IN" ),  0u );
    EXPECT_EQ( tFactory.node_index( "GND" ), 1u );
}

TEST( NgspiceFactory, KeywordValueForm )
{
    Cell< string > tLines = { "t", "R1 in 0 r=1k", ".end" };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );
    ASSERT_EQ( tCircuit->number_of_components(), 1u );
    EXPECT_DOUBLE_EQ( tCircuit->component( 0 )->get_value(), 1.0e3 );
}

TEST( NgspiceFactory, SwitchBeforeVoltageSourceAccepted )
{
    // the dof-counter fix made the layout legal; the v1 refusal is gone
    Cell< string > tLines = {
        "t",
        "* belfem: switch S1 n+=1 n-=0 state=open t_switch=1",
        "V1 1 0 5",
        "R1 1 0 1k",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );
    EXPECT_EQ( tCircuit->number_of_components(), 3u );

    // V before switch stays legal too
    Cell< string > tLegal = {
        "t",
        "V1 1 0 5",
        "R1 1 0 1k",
        "* belfem: switch S1 n+=1 n-=0 state=open t_switch=1",
        ".end"
    };
    NetlistParser tParsedLegal( tLegal );
    NgspiceCircuitFactory tFactoryLegal( tParsedLegal );
    std::unique_ptr< ElectricalCircuit > tCircuitLegal( tFactoryLegal.circuit() );
    EXPECT_EQ( tCircuitLegal->number_of_components(), 3u );
}

TEST( NgspiceFactory, RejectsNonphysicalDeviceValues )
{
    // n=0 would divide the Shockley law by zero
    Cell< string > tN0 = {
        "t", "V1 1 0 5", "D1 1 0 M1", ".model M1 D( n=0 )", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tN0 ) ) ),
                  std::runtime_error );

    // negative saturation current
    Cell< string > tIsNeg = {
        "t", "V1 1 0 5", "D1 1 0 M1", ".model M1 D( is=-1e-14 )", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tIsNeg ) ) ),
                  std::runtime_error );

    // superconductor with a zero length divides the power law
    Cell< string > tSc = {
        "t", "R1 1 0 1k",
        "* belfem: superconductor SC1 n+=1 n-=0 Ic=200 n=20 Ec=1e-4 length=0",
        ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tSc ) ) ),
                  std::runtime_error );

    // an order token that would wrap a 32-bit accumulator
    Cell< string > tWrap = {
        "t", "L1 1 0 1u", "* belfem: order L1 4294967297", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tWrap ) ) ),
                  std::runtime_error );
}

TEST( NgspiceFactory, CurrentSourcePolarity )
{
    // SPICE convention: positive current flows from n+ THROUGH the source
    // to n-, i.e. it is drawn from n+ and delivered into n-. With
    // "I1 0 mid DC 1" through a 1 kOhm to ground, v(mid) = +1000 V.
    Cell< string > tLines = {
        "polarity",
        "I1 0 mid DC 1",
        "R1 mid 0 1k",
        ".end"
    };
    NetlistParser tParsed( tLines );
    NgspiceCircuitFactory tFactory( tParsed );
    const index_t tMid = tFactory.node_index( "mid" );
    std::unique_ptr< ElectricalCircuit > tCircuit( tFactory.circuit() );

    tCircuit->set_timestep( 1.0e-3 );
    tCircuit->shift();
    tCircuit->compute_MNA_matrix();

    real tEpsilon = BELFEM_REAL_MAX;
    uint tIteration = 0;
    while ( tIteration++ < 20 && tEpsilon > 1.0e-10 )
    {
        tCircuit->set_omega( 1.0 );
        tCircuit->compute_jacobian_and_rhs();
        tCircuit->solve();
        tEpsilon = tCircuit->residual();
    }

    EXPECT_LT( tEpsilon, 1.0e-10 );
    EXPECT_NEAR( tCircuit->get_voltage_on_node( tMid ), 1.0e3, 1.0e-6 );
}

//------------------------------------------------------------------------------
// refusals
//------------------------------------------------------------------------------

TEST( NgspiceFactory, RejectsBadOrderDirectives )
{
    // order on a resistor
    Cell< string > tOnR = {
        "t", "R1 1 0 1k", "* belfem: order R1 2", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tOnR ) ) ),
                  std::runtime_error );

    // unknown instance
    Cell< string > tUnknown = {
        "t", "L1 1 0 1u", "* belfem: order L9 2", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tUnknown ) ) ),
                  std::runtime_error );

    // order out of the BDF range
    Cell< string > tSeven = {
        "t", "L1 1 0 1u", "* belfem: order L1 7", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tSeven ) ) ),
                  std::runtime_error );

    // "2k" is a SPICE number but not a plain integer
    Cell< string > tSuffix = {
        "t", "L1 1 0 1u", "* belfem: order L1 2k", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tSuffix ) ) ),
                  std::runtime_error );

    // duplicate
    Cell< string > tDup = {
        "t", "L1 1 0 1u",
        "* belfem: order L1 2", "* belfem: order L1 3", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tDup ) ) ),
                  std::runtime_error );
}

TEST( NgspiceFactory, RejectsValueProblems )
{
    // ic= is value-changing and refused
    Cell< string > tIc = { "t", "C1 1 0 100n ic=5", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tIc ) ) ),
                  std::runtime_error );

    // both positional and keyword value
    Cell< string > tBoth = { "t", "R1 1 0 5k r=5k", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tBoth ) ) ),
                  std::runtime_error );

    // no value at all
    Cell< string > tNone = { "t", "R1 1 0", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tNone ) ) ),
                  std::runtime_error );

    // zero-valued resistor would put 1/R = inf into the MNA
    Cell< string > tZero = { "t", "R1 1 0 0", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tZero ) ) ),
                  std::runtime_error );
}

TEST( NgspiceFactory, RejectsUnsupportedSources )
{
    // SIN with a nonzero offset cannot be represented
    Cell< string > tOffset = { "t", "V1 1 0 SIN(1 5 50)", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tOffset ) ) ),
                  std::runtime_error );

    // PULSE is Phase 6
    Cell< string > tPulse = { "t", "V1 1 0 PULSE(0 5 0 1n 1n 1u 2u)", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tPulse ) ) ),
                  std::runtime_error );

    // a source without any value
    Cell< string > tBare = { "t", "V1 1 0", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tBare ) ) ),
                  std::runtime_error );
}

TEST( NgspiceFactory, RejectsDirectiveProblems )
{
    // unknown kind
    Cell< string > tKind = {
        "t", "R1 1 0 1k", "* belfem: varistor VX n+=1 n-=0", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tKind ) ) ),
                  std::runtime_error );

    // missing key
    Cell< string > tMissing = {
        "t", "R1 1 0 1k",
        "* belfem: superconductor SC1 n+=1 n-=0 Ic=200 n=20 Ec=1e-4", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tMissing ) ) ),
                  std::runtime_error );

    // unknown key
    Cell< string > tExtra = {
        "t", "R1 1 0 1k",
        "* belfem: switch S1 n+=1 n-=0 state=open t_switch=1m bogus=1", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tExtra ) ) ),
                  std::runtime_error );

    // bad switch state
    Cell< string > tState = {
        "t", "R1 1 0 1k",
        "* belfem: switch S1 n+=1 n-=0 state=ajar t_switch=1m", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tState ) ) ),
                  std::runtime_error );
}

TEST( NgspiceFactory, RejectsDuplicateLabelsAndBadModels )
{
    // two cards with one instance name
    Cell< string > tDup = { "t", "R1 1 0 1k", "R1 1 0 2k", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tDup ) ) ),
                  std::runtime_error );

    // diode without its model
    Cell< string > tNoModel = { "t", "V1 1 0 5", "D1 1 0 DMOD", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tNoModel ) ) ),
                  std::runtime_error );

    // model of the wrong type
    Cell< string > tWrongType = {
        "t", "V1 1 0 5", "D1 1 0 M1", ".model M1 NMOS( l=1u )", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tWrongType ) ) ),
                  std::runtime_error );

    // model parameter outside the v1 subset
    Cell< string > tParam = {
        "t", "V1 1 0 5", "D1 1 0 M1", ".model M1 D( is=1e-14 rs=2 )", ".end" };
    EXPECT_THROW( NgspiceCircuitFactory tF( ( NetlistParser( tParam ) ) ),
                  std::runtime_error );
}
