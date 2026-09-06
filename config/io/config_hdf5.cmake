if( USE_HDF5 )
    if( NOT USE_MPI )
        message(FATAL_ERROR "Turn on MPI if you want to link against HDF5" )
    endif()

    belfem_find_package(
            HDF5
            REQUIRED
            HEADERS hdf5.h hdf5_hl.h
            LIBRARIES hdf5_hl hdf5
            INCLUDE_SUFFIXES . hdf5 )

    list( APPEND BELFEM_DEFS "BELFEM_HDF5" )
    list( APPEND BELFEM_INCLUDES "${HDF5_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${HDF5_LIB_DIR}" )
    belfem_link_libraries( 5 ${HDF5_LIBRARIES} )
    belfem_link_libraries( 0 "-lz" "-ldl" )
endif()
