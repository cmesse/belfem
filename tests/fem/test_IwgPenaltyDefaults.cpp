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
 * The penalty slots of IWG_Maxwell: 0 = ghost eta, 1 = ghost k_reg,
 * 2 = Coulomb-gauge chi. The constructor holds fallbacks for paths that never
 * reach Controller::set_params, and the controller writes slot 0 and slot 2
 * even when the deck has no block, because the fallbacks are not the
 * production defaults ( an absent block means 0 ). Both halves of that are
 * pinned: the fallbacks exist, and a write replaces them.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "en_Maxwell_Formulations.hpp"
#include "support/cl_TS_TestStack.hpp"

using namespace belfem;

//------------------------------------------------------------------------------

TEST( IwgPenaltyDefaults, ConstructorFallbacksAreTheDocumentedTriple )
{
    const fem::IWG_Maxwell tIWG( fem::maxwell::Formulation::HPhi,
                                 ModelDimensionality::ThreeD,
                                 false, false );

    EXPECT_EQ( tIWG.penalty( 0 ), 4.0 );
    EXPECT_EQ( tIWG.penalty( 1 ), 1.0e-3 );
    EXPECT_EQ( tIWG.penalty( 2 ), 1.0e-4 );
}

//------------------------------------------------------------------------------

// set_penalty reads the rank of the IWG's field, so this needs the full
// kernel chain; the ghost stack never writes a penalty itself
TEST( IwgPenaltyDefaults, SetPenaltyOverwritesTheFallback )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::test::TS_TestGhostStack tStack ;
    fem::IWG_Maxwell * tIWG = tStack.iwg();

    ASSERT_EQ( tIWG->penalty( 0 ), 4.0 );
    ASSERT_EQ( tIWG->penalty( 2 ), 1.0e-4 );

    // what the controller writes for an absent gauge block and an absent
    // ghost block
    tIWG->set_penalty( 0.0, 2 );
    tIWG->set_penalty( 0.0, 0 );

    EXPECT_EQ( tIWG->penalty( 2 ), 0.0 );
    EXPECT_EQ( tIWG->penalty( 0 ), 0.0 );

    // the neighbor slot is untouched
    EXPECT_EQ( tIWG->penalty( 1 ), 1.0e-3 );
}
