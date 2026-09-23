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
 * The resistivity window [ gRhoMin, gRhoMax ] as the communicator leaves it.
 * The clamp in calculator::MaxwellData::compute_rho() passes both globals to
 * std::clamp, which is undefined on an inverted window, and checks the
 * ordering nowhere; the ordering is a precondition on every writer
 * ( globals.hpp ). This pins the one writer in the open tree.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

#include "typedefs.hpp"
#include "globals.hpp"

using namespace belfem;

TEST( RhoWindow, CommunicatorLeavesAnOrderedWideWindow )
{
    // the test main has run Communicator::init(), which sets the globals
    ASSERT_TRUE( std::isfinite( gRhoMin ) );
    ASSERT_TRUE( std::isfinite( gRhoMax ) );
    // fatal: an inverted window must not reach std::clamp below
    ASSERT_LT( gRhoMin, gRhoMax );

    // the upper bound only catches the power law's runaway; it must sit
    // orders of magnitude above the physical range checked below
    EXPECT_GE( gRhoMax, 1e6 );

    // the clamp and its flag, evaluated with the calculator's two
    // expressions, are the identity on [ 1e-12, 1e3 ) Ohm*m
    for ( int e = -12 ; e <= 2 ; ++e )
    {
        for ( const real tMantissa : { 1.0, 2.5, 7.9 } )
        {
            const real tRho = tMantissa * std::pow( 10.0, e ) ;
            const bool tClamped = ( tRho < gRhoMin ) || ( tRho > gRhoMax ) ;

            EXPECT_FALSE( tClamped ) << "rho = " << tRho ;
            EXPECT_EQ( std::clamp( tRho, gRhoMin, gRhoMax ), tRho ) << "rho = " << tRho ;
        }
    }
}
