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
 * Unit tests for src/circuit/cl_NetlistParser.{hpp,cpp}: the lexer that
 * turns an ngspice netlist into the IR the netlist factory consumes.
 * Comment/continuation/title rules follow the ngspice manual; the v1
 * hard-error policy is todo/ngspice_parser_plan.md §12.2.
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_NetlistParser.hpp"

using namespace belfem;
using electronics::NetlistParser;

//------------------------------------------------------------------------------
// basics
//------------------------------------------------------------------------------

TEST( NetlistParser, ShowcaseShapedDeck )
{
    // the examples/tapestack_circuit lumped part as it will look after conversion
    Cell< string > tLines = {
        "* tape stack demo",
        "Is 1 0 SIN(0 500 10)",
        "L1 1 0 2.5u",
        "C1 1 0 100",
        "* belfem: order L1 2",
        "* belfem: order C1 2",
        ".tran 0.1m 150m",
        ".end"
    };
    NetlistParser tParser( tLines );

    EXPECT_EQ( tParser.title(), "* tape stack demo" );
    ASSERT_EQ( tParser.elements().size(), 3u );
    EXPECT_EQ( tParser.elements()( 0 ).mType, 'i' );
    EXPECT_EQ( tParser.elements()( 0 ).mName, "is" );
    EXPECT_EQ( tParser.elements()( 1 ).mType, 'l' );
    EXPECT_EQ( tParser.elements()( 2 ).mType, 'c' );

    // SIN(0 500 10) normalizes to positional values
    const auto & tSource = tParser.elements()( 0 );
    ASSERT_EQ( tSource.mValues.size(), 4u );
    EXPECT_EQ( tSource.mValues( 0 ), "sin" );
    EXPECT_EQ( tSource.mValues( 1 ), "0" );
    EXPECT_EQ( tSource.mValues( 2 ), "500" );
    EXPECT_EQ( tSource.mValues( 3 ), "10" );

    // two order directives
    ASSERT_EQ( tParser.directives().size(), 2u );
    EXPECT_EQ( tParser.directives()( 0 ).mKind, "order" );
    ASSERT_EQ( tParser.directives()( 0 ).mArgs.size(), 2u );
    EXPECT_EQ( tParser.directives()( 0 ).mArgs( 0 ), "l1" );
    EXPECT_EQ( tParser.directives()( 0 ).mArgs( 1 ), "2" );

    // .tran stored ( hphirun ignores it later; the parser keeps it )
    ASSERT_EQ( tParser.controls().size(), 1u );
    EXPECT_EQ( tParser.controls()( 0 ).mKeyword, "tran" );
    ASSERT_EQ( tParser.controls()( 0 ).mArgs.size(), 2u );
    EXPECT_EQ( tParser.controls()( 0 ).mArgs( 0 ), "0.1m" );

    EXPECT_TRUE( tParser.has_end() );

    // node names in first-appearance order ( O9 )
    ASSERT_EQ( tParser.node_names().size(), 2u );
    EXPECT_EQ( tParser.node_names()( 0 ), "1" );
    EXPECT_EQ( tParser.node_names()( 1 ), "0" );
}

TEST( NetlistParser, TitleLineIsAlwaysConsumed )
{
    // a deck that "starts" with a card loses that card to the title --
    // the classic failure the manual warns about
    Cell< string > tLines = {
        "R1 1 0 5k",
        "R2 1 0 5k",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mName, "r2" );
}

TEST( NetlistParser, EndIsOptionalAtEof )
{
    Cell< string > tLines = {
        "title",
        "R1 in out 5k"
    };
    NetlistParser tParser( tLines );
    EXPECT_EQ( tParser.elements().size(), 1u );
    EXPECT_FALSE( tParser.has_end() );
}

TEST( NetlistParser, ContentAfterEndIsIgnored )
{
    Cell< string > tLines = {
        "title",
        "R1 1 0 5k",
        ".end",
        "Qxx 1 2 3 broken_card_that_would_error"
    };
    NetlistParser tParser( tLines );
    EXPECT_EQ( tParser.elements().size(), 1u );
    EXPECT_TRUE( tParser.has_end() );
}

//------------------------------------------------------------------------------
// comments and continuations
//------------------------------------------------------------------------------

TEST( NetlistParser, CommentForms )
{
    Cell< string > tLines = {
        "title",
        "* a full-line comment",
        "R1 1 0 5k $ dollar comment",
        "R2 1 0 5k ; semicolon comment",
        "R3 1 0 5k // slash comment",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 3u );
    for ( index_t k = 0; k < 3; ++k )
    {
        EXPECT_EQ( tParser.elements()( k ).mValues.size(), 1u );
        EXPECT_EQ( tParser.elements()( k ).mValues( 0 ), "5k" );
    }
}

TEST( NetlistParser, SlashCommentOnlyAtTokenStart )
{
    // '//' inside a token is not a comment
    Cell< string > tLines = {
        "title",
        "V1 a//b 0 5",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mNodes( 0 ), "a//b" );
}

TEST( NetlistParser, ContinuationJoins )
{
    Cell< string > tLines = {
        "title",
        "Vin 3 0",
        "+ SIN(0 1 100MEG",
        "+ 1NS 1E10)",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 1u );
    const auto & tSource = tParser.elements()( 0 );
    ASSERT_EQ( tSource.mValues.size(), 6u );
    EXPECT_EQ( tSource.mValues( 0 ), "sin" );
    EXPECT_EQ( tSource.mValues( 3 ), "100meg" );
    EXPECT_EQ( tSource.mValues( 5 ), "1e10" );
}

TEST( NetlistParser, CommentBetweenCardAndContinuation )
{
    Cell< string > tLines = {
        "title",
        "Vin 3 0",
        "* interleaved comment",
        "",
        "+ 5",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 1u );
    ASSERT_EQ( tParser.elements()( 0 ).mValues.size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mValues( 0 ), "5" );
}

//------------------------------------------------------------------------------
// tokenizing details
//------------------------------------------------------------------------------

TEST( NetlistParser, KwargsWithAndWithoutSpaces )
{
    Cell< string > tLines = {
        "title",
        "R1 1 0 r=5k",
        "R2 1 0 r = 5k",
        "C1 1 0 100n ic = 0",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 3u );
    EXPECT_EQ( tParser.elements()( 0 ).mKwargs( "r" ), "5k" );
    EXPECT_EQ( tParser.elements()( 1 ).mKwargs( "r" ), "5k" );
    EXPECT_EQ( tParser.elements()( 2 ).mValues( 0 ), "100n" );
    EXPECT_EQ( tParser.elements()( 2 ).mKwargs( "ic" ), "0" );
}

TEST( NetlistParser, CaseFoldingAndCrlf )
{
    Cell< string > tLines = {
        "title",
        "RLOAD IN GND 5K\r",
        ".END\r"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mName, "rload" );
    EXPECT_EQ( tParser.elements()( 0 ).mNodes( 0 ), "in" );
    EXPECT_EQ( tParser.elements()( 0 ).mNodes( 1 ), "gnd" );
    EXPECT_EQ( tParser.elements()( 0 ).mValues( 0 ), "5k" );
    EXPECT_TRUE( tParser.has_end() );
}

TEST( NetlistParser, DiodeWithModelCard )
{
    Cell< string > tLines = {
        "title",
        "D1 anode cathode DMOD",
        ".model DMOD D( is=1e-14 n=1.5 )",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mType, 'd' );
    EXPECT_EQ( tParser.elements()( 0 ).mValues( 0 ), "dmod" );
    ASSERT_EQ( tParser.models().size(), 1u );
    EXPECT_EQ( tParser.models()( 0 ).mName, "dmod" );
    EXPECT_EQ( tParser.models()( 0 ).mType, "d" );
    EXPECT_EQ( tParser.models()( 0 ).mKwargs( "is" ), "1e-14" );
    EXPECT_EQ( tParser.models()( 0 ).mKwargs( "n" ), "1.5" );
}

TEST( NetlistParser, DirectiveWithKwargs )
{
    Cell< string > tLines = {
        "title",
        "R1 1 0 5k",
        "* belfem: superconductor SC1 n+=1 n-=0 Ic=200 n=20 Ec=1e-4 length=1",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.directives().size(), 1u );
    const auto & tDirective = tParser.directives()( 0 );
    EXPECT_EQ( tDirective.mKind, "superconductor" );
    ASSERT_EQ( tDirective.mArgs.size(), 1u );
    EXPECT_EQ( tDirective.mArgs( 0 ), "sc1" );
    EXPECT_EQ( tDirective.mKwargs( "ic" ), "200" );
    EXPECT_EQ( tDirective.mKwargs( "n+" ), "1" );
    EXPECT_EQ( tDirective.mKwargs( "n-" ), "0" );
}

TEST( NetlistParser, NodeNamesFirstAppearanceOrder )
{
    Cell< string > tLines = {
        "title",
        "R1 out mid 1k",
        "R2 mid 0 1k",
        "V1 in 0 5",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.node_names().size(), 4u );
    EXPECT_EQ( tParser.node_names()( 0 ), "out" );
    EXPECT_EQ( tParser.node_names()( 1 ), "mid" );
    EXPECT_EQ( tParser.node_names()( 2 ), "0" );
    EXPECT_EQ( tParser.node_names()( 3 ), "in" );
}

TEST( NetlistParser, DirectiveNodesJoinTheNodeList )
{
    // a node that exists ONLY on a directive must still enter the O9
    // packing list, in written order
    Cell< string > tLines = {
        "title",
        "* belfem: superconductor SC1 n+=tape n-=0 Ic=200 n=20 Ec=1e-4 length=1",
        "R1 c 0 5k",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.node_names().size(), 3u );
    EXPECT_EQ( tParser.node_names()( 0 ), "tape" );
    EXPECT_EQ( tParser.node_names()( 1 ), "0" );
    EXPECT_EQ( tParser.node_names()( 2 ), "c" );
}

TEST( NetlistParser, ElementNodesPrecedeLaterDirectiveNodes )
{
    // O9 regression: a card's nodes enter the list before a LATER
    // directive's nodes -- the directive must not jump the queue while
    // the card is still pending
    Cell< string > tLines = {
        "title",
        "R1 a 0 1k",
        "* belfem: switch S1 n+=b n-=0 state=open t_switch=1",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.node_names().size(), 3u );
    EXPECT_EQ( tParser.node_names()( 0 ), "a" );
    EXPECT_EQ( tParser.node_names()( 1 ), "0" );
    EXPECT_EQ( tParser.node_names()( 2 ), "b" );
}

TEST( NetlistParser, DirectiveEndsContinuationChain )
{
    // a directive is a statement: it completes the pending card, so a
    // '+' after it is an orphan
    Cell< string > tLines = {
        "title",
        "Vin 3 0",
        "* belfem: order L1 2",
        "+ 5",
        ".end"
    };
    EXPECT_THROW( NetlistParser tP( tLines ), std::runtime_error );
}

TEST( NetlistParser, DirectiveMarkerVariants )
{
    // whitespace before the colon is tolerated; "belfem" without a colon
    // is an ordinary comment; a directive-shaped first line is the title
    Cell< string > tLines = {
        "* belfem: order L1 2",
        "* belfem : order C1 2",
        "* belfem is great",
        "L1 1 0 2.5u",
        "C1 1 0 100",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.directives().size(), 1u );
    EXPECT_EQ( tParser.directives()( 0 ).mArgs( 0 ), "c1" );
    EXPECT_EQ( tParser.title(), "* belfem: order L1 2" );
}

TEST( NetlistParser, SwitchDirectiveShape )
{
    // plan §5 shape for the timed switch
    Cell< string > tLines = {
        "title",
        "R1 1 2 5k",
        "* belfem: switch S1 n+=1 n-=2 state=open t_switch=20m",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.directives().size(), 1u );
    const auto & tSwitch = tParser.directives()( 0 );
    EXPECT_EQ( tSwitch.mKind, "switch" );
    EXPECT_EQ( tSwitch.mArgs( 0 ), "s1" );
    EXPECT_EQ( tSwitch.mKwargs( "state" ), "open" );
    EXPECT_EQ( tSwitch.mKwargs( "t_switch" ), "20m" );
}

TEST( NetlistParser, ParensAndCommasAreSeparators )
{
    Cell< string > tLines = {
        "title",
        "R1(1)(0)(5k)",
        "R2,1,0,5k",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements().size(), 2u );
    for ( index_t k = 0; k < 2; ++k )
    {
        EXPECT_EQ( tParser.elements()( k ).mNodes( 0 ), "1" );
        EXPECT_EQ( tParser.elements()( k ).mNodes( 1 ), "0" );
        EXPECT_EQ( tParser.elements()( k ).mValues( 0 ), "5k" );
    }
}

TEST( NetlistParser, SlashInsideValueTokenSurvives )
{
    // the threat the comment-strip guards against: '//' inside a token
    Cell< string > tLines = {
        "title",
        "V1 1 0 1e-//2",
        ".end"
    };
    NetlistParser tParser( tLines );
    ASSERT_EQ( tParser.elements()( 0 ).mValues.size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mValues( 0 ), "1e-//2" );
}

//------------------------------------------------------------------------------
// hard errors
//------------------------------------------------------------------------------

TEST( NetlistParser, RejectsUnsupportedElements )
{
    // controlled switch
    Cell< string > tSwitch = { "t", "S1 1 0 2 0 SW", ".end" };
    EXPECT_THROW( NetlistParser tP( tSwitch ), std::runtime_error );

    // subcircuit instance
    Cell< string > tSub = { "t", "X1 1 0 opamp", ".end" };
    EXPECT_THROW( NetlistParser tP( tSub ), std::runtime_error );

    // transistor
    Cell< string > tQ = { "t", "Q1 1 2 3 BJT", ".end" };
    EXPECT_THROW( NetlistParser tP( tQ ), std::runtime_error );

    // behavioral source
    Cell< string > tB = { "t", "B1 1 0 V=V(2)*3", ".end" };
    EXPECT_THROW( NetlistParser tP( tB ), std::runtime_error );
}

TEST( NetlistParser, RejectsUnsupportedControlCards )
{
    Cell< string > tSubckt = { "t", "R1 1 0 5k", ".subckt foo 1 2", ".end" };
    EXPECT_THROW( NetlistParser tP( tSubckt ), std::runtime_error );

    Cell< string > tControl = { "t", "R1 1 0 5k", ".control", ".end" };
    EXPECT_THROW( NetlistParser tP( tControl ), std::runtime_error );

    Cell< string > tIc = { "t", "R1 1 0 5k", ".ic v(1)=5", ".end" };
    EXPECT_THROW( NetlistParser tP( tIc ), std::runtime_error );

    Cell< string > tParam = { "t", "R1 1 0 5k", ".param x=5", ".end" };
    EXPECT_THROW( NetlistParser tP( tParam ), std::runtime_error );
}

TEST( NetlistParser, EndVariantsDoNotTerminate )
{
    // .endc/.ends/.endif must NOT match .end -- they reach the
    // unsupported-control error instead of silently ending the deck
    Cell< string > tEndc = { "t", "R1 1 0 5k", ".endc", ".end" };
    EXPECT_THROW( NetlistParser tP( tEndc ), std::runtime_error );

    Cell< string > tEnds = { "t", "R1 1 0 5k", ".ends", ".end" };
    EXPECT_THROW( NetlistParser tP( tEnds ), std::runtime_error );

    Cell< string > tEndif = { "t", "R1 1 0 5k", ".endif", ".end" };
    EXPECT_THROW( NetlistParser tP( tEndif ), std::runtime_error );

    // and .end with trailing tokens is a hard error, not a termination
    Cell< string > tJunk = { "t", "R1 1 0 5k", ".end nonsense" };
    EXPECT_THROW( NetlistParser tP( tJunk ), std::runtime_error );
}

TEST( NetlistParser, RejectsMalformedParameters )
{
    Cell< string > tTrailing = { "t", "R1 1 0 5k=", ".end" };
    EXPECT_THROW( NetlistParser tP( tTrailing ), std::runtime_error );

    Cell< string > tLeading = { "t", "R1 1 0 =5", ".end" };
    EXPECT_THROW( NetlistParser tP( tLeading ), std::runtime_error );

    Cell< string > tDouble = { "t", "R1 1 0 r==5", ".end" };
    EXPECT_THROW( NetlistParser tP( tDouble ), std::runtime_error );

    // last-wins would silently change the circuit
    Cell< string > tDup = { "t", "R1 1 0 r=5k r=10k", ".end" };
    EXPECT_THROW( NetlistParser tP( tDup ), std::runtime_error );

    // a card of only separators normalizes to zero tokens
    Cell< string > tParen = { "t", "R1 1 0 5k", "(", ".end" };
    EXPECT_THROW( NetlistParser tP( tParen ), std::runtime_error );
}

TEST( NetlistParser, DollarNodeNameIsCutAsComment )
{
    // default-mode ngspice: '$' starts a comment anywhere, so a
    // PSPICE-style node "$G_VDD" truncates the card and errors loudly
    // ( PSPICE compatibility decks are documented as unsupported )
    Cell< string > tLines = { "t", "R1 $G_VDD 0 5k", ".end" };
    EXPECT_THROW( NetlistParser tP( tLines ), std::runtime_error );
}

TEST( NetlistParser, ErrorMessagesCarrySourceAndLine )
{
    // the §12.2 contract: source, line and card in every parser error
    Cell< string > tLines = { "t", "Q1 1 2 3 BJT", ".end" };
    try
    {
        NetlistParser tParser( tLines );
        FAIL() << "expected a parse error";
    }
    catch ( const std::runtime_error & aError )
    {
        const string tWhat( aError.what() );
        EXPECT_NE( tWhat.find( "<memory>:2" ), string::npos ) << tWhat;
        // the card text is quoted as written ( raw, unfolded )
        EXPECT_NE( tWhat.find( "Q1 1 2 3 BJT" ), string::npos ) << tWhat;
    }
}

TEST( NetlistParser, FileConstructorSmoke )
{
    // one on-disk round trip through the Ascii path. The path must be
    // ABSOLUTE: Ascii resolves relative paths against getenv("PWD")
    // ( cl_Ascii.cpp ), which is a shell variable and goes stale when
    // ctest/make change the working directory -- an ofstream would then
    // write where Ascii does not look
    const string tPath = ::testing::TempDir() + "netlist_parser_smoke.cir";
    {
        std::ofstream tFile( tPath );
        tFile << "smoke test deck\n"
              << "R1 in 0 5k ; load\n"
              << ".end\n";
    }
    NetlistParser tParser( tPath );
    std::remove( tPath.c_str() );

    ASSERT_EQ( tParser.elements().size(), 1u );
    EXPECT_EQ( tParser.elements()( 0 ).mName, "r1" );
    EXPECT_EQ( tParser.elements()( 0 ).mValues( 0 ), "5k" );
    EXPECT_TRUE( tParser.has_end() );
}

TEST( NetlistParser, RejectsStructuralErrors )
{
    // continuation with nothing to continue
    Cell< string > tOrphan = { "t", "+ 5", ".end" };
    EXPECT_THROW( NetlistParser tP( tOrphan ), std::runtime_error );

    // element card with too few tokens
    Cell< string > tShort = { "t", "R1 1", ".end" };
    EXPECT_THROW( NetlistParser tP( tShort ), std::runtime_error );

    // kwarg where a node is required
    Cell< string > tNode = { "t", "R1 r=5 0 5k", ".end" };
    EXPECT_THROW( NetlistParser tP( tNode ), std::runtime_error );

    // empty deck and title-only deck
    Cell< string > tEmpty;
    EXPECT_THROW( NetlistParser tP( tEmpty ), std::runtime_error );
    Cell< string > tTitleOnly = { "t" };
    EXPECT_THROW( NetlistParser tP( tTitleOnly ), std::runtime_error );

    // empty directive
    Cell< string > tDirective = { "t", "R1 1 0 5k", "* belfem:", ".end" };
    EXPECT_THROW( NetlistParser tP( tDirective ), std::runtime_error );
}
