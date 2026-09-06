# set the name of the test
set( TESTNAME io )

# list the test sources
set( SOURCES
        test_HDF5.cpp
        test_filetools.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
