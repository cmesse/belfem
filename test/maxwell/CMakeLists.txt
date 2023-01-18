set( TESTNAME maxwell )

set( SOURCES
        test_tri3.cpp
        test_tri6.cpp
        test_tet4.cpp
        test_tet10.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
if ( USE_GASMODELS )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
endif()
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/edgeinterpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/maxwell )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )

if ( USE_GASMODELS )

    set ( LIBLIST
            sparse
            bezier
            mesh
            materials
            integration
            ode
            interpolation
            edgeinterpolation
            gastables
            gasmodels
            iwg
            kernel
            maxwell
            postproc
            integration )
else()
    set ( LIBLIST
            sparse
            bezier
            mesh
            materials
            integration
            ode
            interpolation
            edgeinterpolation
            iwg
            kernel
            maxwell
            kernel
            postproc
            integration )
endif()

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_maxwell.hdf5
        DESTINATION /tmp )