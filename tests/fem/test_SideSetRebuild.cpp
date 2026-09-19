/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * A sideset calculator that changes its integration order after
 * construction must evaluate its edge functions at the points of the
 * rebuilt lookup tables, not at the points of the tables it was built with.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "support/cl_TS_TestStack.hpp"

using namespace belfem;

TEST( SideSetRebuild, EdgeFunctionsFollowTheNewOrder )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }
#if ! BELFEM_ASSERTIONS_ACTIVE
    GTEST_SKIP() << "the stale precompute is caught by the bounds assertion, "
                    "which this build compiles out" ;
#endif

    fem::test::TS_TestGhostStack tStack( 0.1, false );
    fem::Calculator * tCalc = tStack.calculator();
    fem::Element * tGhost = tStack.ghost_element();
    const uint tMasterIndex = tGhost->facet()->index_on_master();

    // the fixture built the group's tables at the automatic order
    const uint tPointsBefore =
        tCalc->group()->master_integration( tMasterIndex )->weights().length();

    // an order whose facet rule has more points than the automatic one
    const uint tOrder = 13 ;
    tCalc->set_integration_order( tOrder );
    const uint tPointsAfter =
        tCalc->group()->master_integration( tMasterIndex )->weights().length();
    ASSERT_GT( tPointsAfter, tPointsBefore )
        << "the test needs an order that enlarges the facet rule" ;

    // the element link pairs the edge function with the rebuilt table
    tCalc->link( tGhost );

    // the last point of the new rule exists in the edge function's
    // precomputed factors only if they were computed after the rebuild
    EXPECT_NO_THROW( tCalc->Em( tPointsAfter - 1 ) )
        << "master edge function precomputed from the tables of the old order" ;
    EXPECT_NO_THROW( tCalc->Es( tPointsAfter - 1 ) )
        << "slave edge function precomputed from the tables of the old order" ;
}

//------------------------------------------------------------------------------

// the calculator's workspace vectors survive an order change: the map
// hands out the same object, so a reference bound to it stays valid. This
// checks the map on a sideset calculator; it does not exercise the
// MaxwellData references, which exist on block calculators only
TEST( SideSetRebuild, VectorsSurviveAnOrderChange )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::test::TS_TestGhostStack tStack( 0.1, false );
    fem::Calculator * tCalc = tStack.calculator();
    const string & tLabel = tStack.iwg()->all_fields()( 0 );
    ASSERT_TRUE( tCalc->vector_exists( tLabel ) );

    const Vector< real > * tBefore = & tCalc->vector( tLabel );
    tCalc->set_integration_order( 13 );
    const Vector< real > * tAfter = & tCalc->vector( tLabel );

    EXPECT_EQ( tBefore, tAfter )
        << "workspace vector " << tLabel << " was replaced by the order change" ;
}
