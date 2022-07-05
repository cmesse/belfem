# List source files
set( TESTNAME linalg )

set( SOURCES
        cl_Vector.cpp
        cl_Matrix.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )