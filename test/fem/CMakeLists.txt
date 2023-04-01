
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_database.hdf5
        DESTINATION /tmp )

set( TESTNAME geometry )

set( SOURCES
        fn_normal_tri3.cpp
        fn_normal_tri6.cpp
        fn_normal_quad4.cpp
        fn_normal_quad9.cpp
        fn_normal_tet4.cpp
        fn_normal_tet10.cpp
        fn_normal_hex8.cpp
        fn_normal_hex27.cpp
        cl_IntegrationData_Interface.cpp
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
            postproc
            kernel
            maxwell
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
            postproc
            kernel
            maxwell
            kernel
            integration )
endif()


# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
