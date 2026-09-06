# -------------------------------------------------------------------------
# GPERFTOOLS : Make sure that paths are set correctly.
# -------------------------------------------------------------------------

if (USE_PROFILER)
    list( APPEND BELFEM_DEFS "BELFEM_PROFILER" )
    belfem_link_libraries( 2 "-lprofiler" )
    if( NOT USE_DEBUG )
        set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -g" )
        set( BELFEM_CFLAGS "${BELFEM_CFLAGS} -g" )
    endif()
endif()