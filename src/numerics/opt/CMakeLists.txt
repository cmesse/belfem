set( LIBNAME opt )

set( SOURCES
        cl_Optimizer.cpp )


set ( LIBLIST
        opt )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

include_directories( ${BELFEM_SOURCE_DIR}/numerics/opt )

# the example needs a working solver backend, so it is only built with NLOPT
if( USE_EXAMPLES AND USE_NLOPT )
    set( EXECNAME opt )
    set( MAIN     opt.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()
