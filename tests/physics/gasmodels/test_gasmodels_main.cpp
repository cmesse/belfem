//
// Created by Christian Messe on 12.09.19.
//


#include <gtest/gtest.h>
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "assert.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 3 );

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