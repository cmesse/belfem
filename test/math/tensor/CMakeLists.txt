# List source files
set( TESTNAME tensor )

set( SOURCES
        fn_ddot_contract42.cpp
        fn_ddot_contract44.cpp
        fn_rotate.cpp
        fn_kelvin_christoffel.cpp
        fn_invert_symmetric.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/math/tensor )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
#set ( LIBLIST
#        <add lib here> )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )