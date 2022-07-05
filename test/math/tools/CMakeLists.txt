# List source files
set( TESTNAME mathtools )

set( SOURCES
        fn_cardano.cpp
        fn_cubic_bezier.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/math/tools )

#set ( LIBLIST
#        <add lib here> )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )