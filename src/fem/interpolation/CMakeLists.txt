set( LIBNAME interpolation )

set( SOURCES
        cl_IF_InterpolationFunctionFactory.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )

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