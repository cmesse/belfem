if( USE_EXODUS )
    if( USE_HDF5 )
        list( APPEND BELFEM_DEFS "BELFEM_EXODUS" )
        set( BELFEM_IO_LIBS "-lexodus -lnetcdf ${BELFEM_IO_LIBS}" )
    else()
        message(FATAL_ERROR "Turn on HDF5 if you want to link against ExodusII." )
    endif()
endif()