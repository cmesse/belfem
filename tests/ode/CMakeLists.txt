# set the name of the test
set( TESTNAME ode )

# ODE module needs its own include directory
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )

# list the test sources
set( SOURCES
        test_BDF.cpp
        test_ODE_Integrator.cpp
        )

set( LIBLIST ode )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
