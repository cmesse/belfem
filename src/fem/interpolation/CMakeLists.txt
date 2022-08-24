set( LIBNAME interpolation )

set( SOURCES
        cl_IF_InterpolationFunctionFactory.cpp
        cl_IF_IntegrationData.cpp
        fn_IF_initialize_integration_points.cpp
        fn_IF_initialize_integration_points_on_facet.cpp
        fn_IF_initialize_shape_function.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )

set ( LIBLIST
        mesh
        interpolation
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES )
    set( EXECNAME interpolation )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()