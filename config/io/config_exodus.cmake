if( USE_EXODUS )
    if( NOT USE_HDF5 )
        message(FATAL_ERROR "Turn on HDF5 if you want to link against ExodusII." )
    endif()

    belfem_find_package(
            EXODUS
            REQUIRED
            HEADERS exodusII.h
            LIBRARIES exodus netcdf )

    list( APPEND BELFEM_DEFS "BELFEM_EXODUS" )
    list( APPEND BELFEM_INCLUDES "${EXODUS_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${EXODUS_LIB_DIR}" )

    # exodus (rank 7) sits above netcdf (rank 6) in the dependency graph
    list( GET EXODUS_LIBRARIES 0 EXODUS_EXODUS_LIBRARY )
    list( GET EXODUS_LIBRARIES 1 EXODUS_NETCDF_LIBRARY )
    belfem_link_libraries( 7 ${EXODUS_EXODUS_LIBRARY} )
    belfem_link_libraries( 6 ${EXODUS_NETCDF_LIBRARY} )
endif()
