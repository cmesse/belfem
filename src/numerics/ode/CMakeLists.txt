set( LIBNAME ode )

set( SOURCES
        cl_ODE_Integrator.cpp
        fn_ODE_RK45.cpp )


set ( LIBLIST
        ode )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )

if( USE_EXAMPLES )
    set( EXECNAME ode )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()