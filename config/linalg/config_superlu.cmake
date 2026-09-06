if( USE_SUPERLU )
    belfem_find_package(
            SUPERLU
            REQUIRED
            HEADERS slu_ddefs.h
            LIBRARIES superlu
            INCLUDE_SUFFIXES . superlu )

    list( APPEND BELFEM_DEFS "BELFEM_SUPERLU" )

    list(FIND BELFEM_DEFS "BELFEM_METIS" ITEM_INDEX )
    if(ITEM_INDEX EQUAL -1)
        list(APPEND BELFEM_DEFS "BELFEM_METIS")
    endif()

    list( APPEND BELFEM_INCLUDES "${SUPERLU_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${SUPERLU_LIB_DIR}" )
    belfem_link_libraries( 4 ${SUPERLU_LIBRARIES} )

    if( USE_MATRIX_ARMADILLO )
        # BELFEM ships its own sequential SuperLU wrapper (cl_SolverSUPERLU).
        # Disable Armadillo's bundled SuperLU to avoid a fatal include-guard
        # collision on superlu_enum_consts.h (Armadillo includes it inside
        # namespace arma::superlu, defining the global guard and hiding fact_t
        # etc. from the global <slu_ddefs.h>).
        list( APPEND BELFEM_DEFS "ARMA_DONT_USE_SUPERLU" )
    endif ()
endif()
