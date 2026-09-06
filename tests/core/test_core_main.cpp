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
 * Test main for core unit tests.
 */

#include <gtest/gtest.h>
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "assert.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 5 );

int
main( int    argc,
      char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // error paths are unit-tested with EXPECT_THROW, so this binary asks
    // for the throwing reaction even when built with NDEBUG
    belfem::assert::set_throw_on_error( true );
    // error paths are exercised on purpose here - keep them out of syslog
    belfem::assert::set_syslog_on_error( false );

    // start test session
    testing::InitGoogleTest( &argc, argv );

    // run the tests
    int aResult = RUN_ALL_TESTS();

    // close communicator
    gComm.finalize();

    // return the test result
    return aResult;
}
