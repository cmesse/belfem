# set the name of the test
set( TESTNAME mesh )
set( TESTLABELS fast )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )

# list the test sources
set( SOURCES
        test_MeshEntities.cpp
        test_Mesh.cpp
        )

set( LIBLIST graph bezier interpolation integration mesh )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
