# List source files
set( TESTNAME spline )

set( SOURCES
        cl_Spline.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )