# set the name of the test
set( TESTNAME core )

# list the test sources
set( SOURCES
        stringtools.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )