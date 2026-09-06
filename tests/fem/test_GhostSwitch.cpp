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

#include "typedefs.hpp"
#include "cl_InputFile.hpp"
#include "fn_FEM_ghost_switch.hpp"
#include "fn_mesh_config_tag.hpp"

namespace
{
    // a minimal deck in the test temp directory; the parser needs a file.
    // The path must be ABSOLUTE: Ascii resolves relative paths against
    // getenv("PWD") ( cl_Ascii.cpp ), which goes stale when ctest changes
    // the working directory -- the ofstream would then write where Ascii
    // does not look. Removed again by the Deck guard below so the test
    // leaves no litter
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
