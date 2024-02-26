# -------------------------------------------------------------------------
# GPERFTOOLS : Make sure that paths are set correctly.
# -------------------------------------------------------------------------

if (USE_PROFILER)
    list( APPEND BELFEM_DEFS "BELFEM_PROFILER" )
    set( BELFEM_PERFORMANCELIBS "-lprofiler")
    #set( BELFEM_PERFORMANCELIBS "-lprofiler -ltcmalloc")
    if( NOT USE_DEBUG )
        set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -g" )
        set( BELFEM_CFLAGS "${BELFEM_CFLAGS} -g" )
    endif()
endif()