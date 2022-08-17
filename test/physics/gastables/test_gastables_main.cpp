//
// Created by Christian Messe on 2019-08-16.
//

#include <gtest/gtest.h>
#include "cl_Communicator.hpp"

#include "cl_Logger.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 3 );

int
main( int    argc,
      char * argv[] )
{
    // create communicator
    gComm = belfem::Communicator( argc, argv );

    // start test session
    testing::InitGoogleTest( &argc, argv );

    // run the tests
    int aResult = RUN_ALL_TESTS();

    // close communicator
    gComm.finalize();

    // return the test result
    return aResult;
}