set( LIBNAME interpolation )

set( SOURCES
        cl_EdgeFunctionFactory.cpp
        nedelec/cl_EF_TRI3.cpp
        nedelec/cl_EF_TRI6.cpp
        nedelec/cl_EF_TET4.cpp
        nedelec/cl_EF_TET10.cpp
        nedelec/cl_EF_LINE3.cpp
        nedelec/cl_EF_HEX8.cpp
        nedelec/cl_EF_HEX8TS.cpp
        nedelec/cl_EF_HEX8TB.cpp
        nedelec/cl_EF_QUAD4TS.cpp
        nedelec/cl_EF_PENTA6TS.cpp
        cl_IF_InterpolationFunctionFactory.cpp
        cl_IF_IntegrationData.cpp
        cl_IF_ElementMapper.cpp
        fn_IF_initialize_integration_points.cpp
        fn_IF_initialize_integration_points_on_facet.cpp
        fn_IF_initialize_shape_function.cpp
        )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )

set ( LIBLIST
        mesh
        integration
        interpolation
        kernel
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES )
    set( EXECNAME interpolation )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()