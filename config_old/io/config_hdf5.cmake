if( USE_HDF5 )
    list( APPEND BELFEM_DEFS "BELFEM_HDF5" )
    set( BELFEM_IO_LIBS "-lhdf5_hl -lhdf5 -lz -ldl ${BELFEM_IO_LIBS}" )
endif()