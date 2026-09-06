if( USE_PARPACK )
    if( NOT USE_ARPACK )
        message(FATAL_ERROR "Turn on ARPACK if you want to link against PARPACK" )
    endif()
    if( NOT USE_MPI )
        message(FATAL_ERROR "Turn on MPI if you want to link against PARPACK" )
    endif()

    list( APPEND BELFEM_DEFS "BELFEM_PARPACK" )

    # rank 7, one above arpack: parpack calls into arpack, and the link line is
    # assembled from rank 10 down to 0
    belfem_link_libraries( 7 "-lparpack" )
endif ()
