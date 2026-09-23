/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * A terminal pair creates ONE lumped component, and the controller pairs
 * boundary conditions with components by position, so its input curves must
 * form exactly one bracket group. Two groups, written as two brackets or as
 * an unbracketed list ( where every id is its own group ), are refused by
 * name when the factory reads the deck, before any mesh is touched. The
 * one-group spelling is the control; the same deck with "[1,2]" builds in
 * HybridFactory.TerminalPairByNodeName.
 *
 * The hybrid path folds labels to lower case, so the deck's Z1 is named z1
 * in the refusal. All scratch files use absolute paths.
 */

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_ElectricalCircuitFactory.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"

using namespace belfem ;
using electronics::ElectricalCircuitFactory ;

namespace
{
    // the netlist and the deck of one case, removed again on scope exit
    struct HybridDeck
    {
        std::string mCir ;
        std::string mConf ;

        explicit HybridDeck( const std::string & aInputCurves )
        {
            const std::string tStem = ::testing::TempDir() + "terminal_pair_groups_"
                + ::testing::UnitTest::GetInstance()->current_test_info()->name();

            mCir  = tStem + ".cir" ;
            mConf = tStem + ".conf" ;

            std::ofstream tCir( mCir );
            tCir << "coupled deck\n"
                    "Is tape 0 SIN(0 500 10)\n"
                    "L1 tape 0 2.5u\n"
                    "* belfem: order L1 2\n"
                    ".end\n" ;
            tCir.close();

            std::ofstream tConf( mConf );
            tConf << "circuit\n"
                     "{\n"
                     "    file : " << mCir << " ;\n"
                     "    topology\n"
                     "    {\n"
                     "        terminal pair\n"
                     "        {\n"
                     "            label : Z1 ;\n"
                     "            node + : tape ;\n"
                     "            node - : gnd ;\n"
                     "            input curves : " << aInputCurves << " ;\n"
                     "            output curves : [3,4] ;\n"
                     "        }\n"
                     "    }\n"
                     "}\n" ;
        }

        ~HybridDeck()
        {
            std::remove( mCir.c_str() );
            std::remove( mConf.c_str() );
        }
    };

    // requires the factory to refuse the deck with the terminal-pair group
    // guard, naming the lowercased label and the group count
    void
    expect_group_refusal( const std::string & aInputCurves )
    {
        HybridDeck tDeck( aInputCurves );

        Cell< fem::PhysicalBoundaryCondition * > tBCs ;

        bool tThrown = false ;
        try
        {
            ElectricalCircuitFactory tFactory( tDeck.mConf, tBCs );
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;

            const std::string tMessage( aError.what() );
            EXPECT_NE( tMessage.find( "terminal pair 'z1'" ), std::string::npos )
                << "message was: " << tMessage ;
            EXPECT_NE( tMessage.find( "defines 2 bracket groups" ), std::string::npos )
                << "message was: " << tMessage ;
        }
        EXPECT_TRUE( tThrown ) << "input curves : " << aInputCurves ;

        // the guard fires before this pair's condition is created
        EXPECT_EQ( tBCs.size(), 0u );
        for ( auto * tBC : tBCs ) { delete tBC ; }
    }
}

TEST( TerminalPairGroups, TwoBracketGroupsAreRefused )
{
    expect_group_refusal( "[1],[2]" );
}

TEST( TerminalPairGroups, UnbracketedListIsOneGroupPerIdAndRefused )
{
    expect_group_refusal( "1,2" );
}

TEST( TerminalPairGroups, OneBracketGroupBuilds )
{
    HybridDeck tDeck( "[1,2]" );

    Cell< fem::PhysicalBoundaryCondition * > tBCs ;
    EXPECT_NO_THROW( { ElectricalCircuitFactory tFactory( tDeck.mConf, tBCs ); } );

    EXPECT_EQ( tBCs.size(), 1u );
    for ( auto * tBC : tBCs ) { delete tBC ; }
}
