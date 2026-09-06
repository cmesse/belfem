if( USE_ARPACK )
    list( APPEND BELFEM_DEFS "BELFEM_ARPACK" )
    belfem_link_libraries( 6 "-larpack" )
endif ()