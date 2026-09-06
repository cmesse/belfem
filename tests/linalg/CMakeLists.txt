# set the name of the test
set( TESTNAME linalg )
set( TESTLABELS fast )

# list the test sources
set( SOURCES
        test_Vector.cpp
        test_Matrix.cpp
        test_LinalgSolvers.cpp
        test_lapack.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
