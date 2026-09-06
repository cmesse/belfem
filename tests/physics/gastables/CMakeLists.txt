# List source files
set( TESTNAME gastables )

set( SOURCES
        cl_GT_InputThermo.cpp
        cl_GT_InputTransport.cpp
        cl_GT_InputData.cpp
        cl_GT_HeatPoly.cpp
        cl_GT_RefGas.cpp
        cl_GT_RefGas_Consistency.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )

set ( LIBLIST
        gastables )

# whole suite runs in well under a second
set( TESTLABELS fast )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )