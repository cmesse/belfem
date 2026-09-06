# set the name of the test
set( TESTNAME core )
set( TESTLABELS fast )

# list the test sources
set( SOURCES
        test_StringTools.cpp
        test_AvailableMemory.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
