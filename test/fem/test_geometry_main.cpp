//
// Created by christian on 3/14/23.
//
#include <gtest/gtest.h>

#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "cl_HDF5.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 5 );
belfem::HDF5         * gDatabase ;

int
main( int    argc,
      char * argv[] )
{
    // create communicator
    gComm = belfem::Communicator( argc, argv );

    int aResult = 0 ;

    if( gComm.rank() == 0 )
    {
        // start test session
        testing::InitGoogleTest( &argc, argv );

        // load the file
        gDatabase = new belfem::HDF5( "/tmp/test_geometry.hdf5", belfem::FileMode::OPEN_RDONLY );


        // run the tests
        aResult = RUN_ALL_TESTS();

        delete gDatabase ;
    }

    // close communicator
    gComm.finalize();

    // return the test result
    return aResult;
}