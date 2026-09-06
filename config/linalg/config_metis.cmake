if( USE_METIS )
    # guard against duplicates: config_superlu.cmake also sets BELFEM_METIS
    list(FIND BELFEM_DEFS "BELFEM_METIS" ITEM_INDEX )
    if(ITEM_INDEX EQUAL -1)
        list(APPEND BELFEM_DEFS "BELFEM_METIS")
    endif()
    list(FIND BELFEM_DEFS "BELFEM_PARMETIS" ITEM_INDEX )
    if(ITEM_INDEX EQUAL -1)
        list(APPEND BELFEM_DEFS "BELFEM_PARMETIS")
    endif()

    belfem_find_package(
            METIS
            REQUIRED
            HEADERS metis.h parmetis.h
            LIBRARIES parmetis metis
            INCLUDE_SUFFIXES . metis )

    list( APPEND BELFEM_INCLUDES "${METIS_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${METIS_LIB_DIR}" )

    # parmetis (rank 5) sits above metis (rank 3) in the dependency graph
    list( GET METIS_LIBRARIES 0 METIS_PARMETIS_LIBRARY )
    list( GET METIS_LIBRARIES 1 METIS_METIS_LIBRARY )
    belfem_link_libraries( 5 ${METIS_PARMETIS_LIBRARY} )
    belfem_link_libraries( 3 ${METIS_METIS_LIBRARY} )
endif ()
