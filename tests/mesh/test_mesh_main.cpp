/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Test main for mesh module tests.
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
    gComm.init( argc, argv );

    // error paths are unit-tested with EXPECT_THROW, so this binary asks
    // for the throwing reaction even when built with NDEBUG
    belfem::assert::set_throw_on_error( true );
    // error paths are exercised on purpose here - keep them out of syslog
    belfem::assert::set_syslog_on_error( false );
    testing::InitGoogleTest( &argc, argv );
    int aResult = RUN_ALL_TESTS();
    gComm.finalize();
    return aResult;
}
