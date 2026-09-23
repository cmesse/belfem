/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * The thin-shell ghost switch as the deck states it: the one reader shared
 * by the thin-shell factory, the controller and the mesh cache tag
 * ( fn_FEM_ghost_switch.hpp ), and the tag's view of it.
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <string>
#include <unistd.h>

#include "typedefs.hpp"
#include "cl_InputFile.hpp"
#include "fn_FEM_ghost_switch.hpp"
#include "fn_mesh_config_tag.hpp"

namespace
{
    // a minimal deck in the test temp directory; the parser needs a file.
    // The path is absolute, so the case does not depend on the working
    // directory. Removed again by the Deck guard below so the test leaves
    // no litter
    std::string
    write_deck( const std::string & aName, const std::string & aSolverBody )
    {
        const std::string tPath = ::testing::TempDir() + aName ;
        std::ofstream tFile( tPath );
        tFile << "solver\n{\n" << aSolverBody << "\n}\n";
        tFile.close();
        return tPath ;
    }

    struct Deck
    {
        std::string mName ;
        Deck( const std::string & aName, const std::string & aSolverBody ) :
            mName( write_deck( aName, aSolverBody ) ) {}
        ~Deck() { std::remove( mName.c_str() ); }
    };
}

TEST( GhostSwitch, AbsentBlockIsOff )
{
    Deck tGuard( "ghost_absent.conf",
        "    nonlinear magnetic\n    {\n        tolerance : 1e-7 ;\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_DOUBLE_EQ( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ), -1.0 );
    EXPECT_FALSE( belfem::fem::ghost_facets_requested( tDeck ) );
}

TEST( GhostSwitch, EtaZeroIsOff )
{
    Deck tGuard( "ghost_zero.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n"
        "            eta   : 0 ;\n            k_reg : 1e-3 Ohm ;\n        }\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_DOUBLE_EQ( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ), 0.0 );
    EXPECT_FALSE( belfem::fem::ghost_facets_requested( tDeck ) );
}

TEST( GhostSwitch, EtaPositiveIsOn )
{
    Deck tGuard( "ghost_on.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n"
        "            eta   : 4 ;\n            k_reg : 1e-3 Ohm ;\n        }\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_DOUBLE_EQ( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ), 4.0 );
    EXPECT_TRUE( belfem::fem::ghost_facets_requested( tDeck ) );
}

TEST( GhostSwitch, BlockWithoutEtaIsAnError )
{
    // k_reg alone cannot say whether the ghost is wanted
    Deck tGuard( "ghost_kreg_only.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n"
        "            k_reg : 1e-3 Ohm ;\n        }\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_ANY_THROW( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ) );
}

TEST( GhostSwitch, NegativeEtaIsAnError )
{
    Deck tGuard( "ghost_negative.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n"
        "            eta : -1 ;\n        }\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_ANY_THROW( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ) );
}

TEST( GhostSwitch, MagneticSectionWinsOverPlainNonlinear )
{
    // the alias rule of every key in the section: a block under the losing
    // `nonlinear` is ignored when `nonlinear magnetic` exists
    Deck tGuard( "ghost_alias.conf",
        "    nonlinear\n    {\n        nitsche ghost penalty\n        {\n            eta : 4 ;\n        }\n    }\n"
        "    nonlinear magnetic\n    {\n        tolerance : 1e-7 ;\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_DOUBLE_EQ( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ), -1.0 );
}

TEST( GhostSwitch, NoSolverSectionIsOff )
{
    const std::string tPath = ::testing::TempDir() + "ghost_nosolver.conf" ;
    std::ofstream tFile( tPath );
    tFile << "mesh\n{\n    file : x.msh ;\n    unit : mm ;\n}\n";
    tFile.close();

    belfem::InputFile tDeck( tPath );
    EXPECT_FALSE( belfem::fem::ghost_facets_requested( tDeck ) );
    std::remove( tPath.c_str() );
}

TEST( GhostSwitch, PlainNonlinearSectionIsReadWhenAlone )
{
    // the alias rule's other half: with no `nonlinear magnetic`, `nonlinear`
    // is the winning section and its block counts
    Deck tGuard( "ghost_plain.conf",
        "    nonlinear\n    {\n        nitsche ghost penalty\n        {\n            eta : 4 ;\n        }\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_DOUBLE_EQ( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ), 4.0 );
}

TEST( GhostSwitch, DimensionedEtaIsAnError )
{
    // eta is dimensionless; "4 mOhm" must not be silently SI-scaled to 0.004
    Deck tGuard( "ghost_dimensioned.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n            eta : 4 mOhm ;\n        }\n    }" );
    belfem::InputFile tDeck( tGuard.mName );

    EXPECT_ANY_THROW( belfem::fem::read_ghost_eta( tDeck.section( "solver" ) ) );
}

TEST( GhostSwitch, MeshCacheTagSeesTheSwitch )
{
    // two decks that differ ONLY in eta must produce different cache tags:
    // the switch changes the discretization the .bfm stores, so a flipped
    // deck must miss the cache instead of loading the other layout
    Deck tGuardOn( "ghost_tag_on.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n            eta : 4 ;\n        }\n    }" );
    belfem::InputFile tOn( tGuardOn.mName );
    Deck tGuardOff( "ghost_tag_off.conf",
        "    nonlinear magnetic\n    {\n        nitsche ghost penalty\n        {\n            eta : 0 ;\n        }\n    }" );
    belfem::InputFile tOff( tGuardOff.mName );

    const std::string tTextOn  = belfem::fem::maxwell::mesh_config_text( tOn );
    const std::string tTextOff = belfem::fem::maxwell::mesh_config_text( tOff );

    EXPECT_NE( tTextOn.find( "thinshell.ghost = on" ),  std::string::npos );
    EXPECT_NE( tTextOff.find( "thinshell.ghost = off" ), std::string::npos );
    EXPECT_NE( belfem::fem::maxwell::mesh_config_tag( tOn ), belfem::fem::maxwell::mesh_config_tag( tOff ) );
}

//------------------------------------------------------------------------------
// the tag's sensitivity to the ids that feed cut construction. Each pair of
// decks below differs in one place and must produce different texts and
// different tags
//------------------------------------------------------------------------------

namespace
{
    // a whole deck, verbatim; removed again on destruction
    struct RawDeck
    {
        std::string mName ;
        RawDeck( const std::string & aName, const std::string & aBody ) :
            mName( ::testing::TempDir() + std::to_string( ( long ) getpid() ) + "_" + aName )
        {
            std::ofstream tFile( mName );
            tFile << aBody ;
            tFile.close();
        }
        ~RawDeck() { std::remove( mName.c_str() ); }
    };

    std::string
    current_deck( const std::string & aInputTerminals )
    {
        return "boundary conditions\n{\n    current\n    {\n"
               "        input terminals : " + aInputTerminals + " ;\n"
               "        output terminals : 2,4,6 ;\n"
               "        amplitude : 90 A ;\n    }\n}\n" ;
    }

    // two unlabelled terminal pairs under circuit/topology, as in
    // examples/tapestack_circuit minus the labels
    std::string
    circuit_deck( const std::string & aFirstInput, const std::string & aSecondInput )
    {
        return "circuit\n{\n    topology\n    {\n"
               "        terminal pair\n        {\n"
               "            input curves  : " + aFirstInput + " ;\n"
               "            output curves : [7,8] ;\n        }\n"
               "        terminal pair\n        {\n"
               "            input curves  : " + aSecondInput + " ;\n"
               "            output curves : [9,10] ;\n        }\n"
               "    }\n}\n" ;
    }

    void
    expect_distinct( const std::string & aNameA, const std::string & aBodyA,
                     const std::string & aNameB, const std::string & aBodyB )
    {
        RawDeck tGuardA( aNameA, aBodyA );
        RawDeck tGuardB( aNameB, aBodyB );
        belfem::InputFile tA( tGuardA.mName );
        belfem::InputFile tB( tGuardB.mName );

        const std::string tTextA = belfem::fem::maxwell::mesh_config_text( tA );
        const std::string tTextB = belfem::fem::maxwell::mesh_config_text( tB );

        EXPECT_NE( tTextA, tTextB ) << "text:\n" << tTextA ;
        EXPECT_NE( belfem::fem::maxwell::mesh_config_tag( tA ),
                   belfem::fem::maxwell::mesh_config_tag( tB ) )
            << "texts differ but the tags collide:\n" << tTextA << "---\n" << tTextB ;
    }
}

TEST( MeshConfigTag, LastIdOfAListChangesTheTag )
{
    // an id list starts with a digit but is not a number; it must be
    // compared verbatim, not truncated at the first separator
    expect_distinct( "tag_ids_11.conf", current_deck( "1,3,5,7,9,11" ),
                     "tag_ids_13.conf", current_deck( "1,3,5,7,9,13" ) );
}

TEST( MeshConfigTag, OrderWithinAListChangesTheTag )
{
    // same characters, same length: a hash that ignored order would collide
    expect_distinct( "tag_order_a.conf", current_deck( "1,3,5,7,9,11,13" ),
                     "tag_order_b.conf", current_deck( "1,3,5,7,9,13,11" ) );
}

TEST( MeshConfigTag, SwappedUnlabelledPairsChangeTheTag )
{
    // the ids of two unlabelled terminal pairs swapped between them: the
    // lines differ only in which sibling they belong to
    expect_distinct( "tag_pairs_a.conf", circuit_deck( "[1,2]", "[5,6]" ),
                     "tag_pairs_b.conf", circuit_deck( "[5,6]", "[1,2]" ) );
}

TEST( MeshConfigTag, CircuitInputCurvesChangeTheTag )
{
    // the only difference is one terminal pair's input curves
    expect_distinct( "tag_curves_a.conf", circuit_deck( "[1,2]", "[5,6]" ),
                     "tag_curves_b.conf", circuit_deck( "[1,2]", "[5,11]" ) );
}
