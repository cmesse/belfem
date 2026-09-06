if( USE_MUMPS )
    if( NOT USE_MPI )
        message(FATAL_ERROR "Turn on MPI if you want to link against MUMPS" )
    endif()
    if( NOT USE_SCOTCH )
        message(FATAL_ERROR "Turn on SCOTCH if you want to link against MUMPS" )
    endif()
    if( NOT USE_METIS )
        message(FATAL_ERROR "Turn on METIS if you want to link against MUMPS" )
    endif()

    belfem_find_package(
            MUMPS
            REQUIRED
            HEADERS dmumps_c.h
            LIBRARIES dmumps mumps_common esmumps pord
            INCLUDE_SUFFIXES . mumps )

    list( APPEND BELFEM_DEFS "BELFEM_MUMPS" )
    list( APPEND BELFEM_INCLUDES "${MUMPS_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${MUMPS_LIB_DIR}" )

    belfem_link_libraries( 6 ${MUMPS_LIBRARIES} )

    # the fortran interface of MPI; assume openmpi (mpich needs -lmpifort)
    belfem_link_libraries( 4 "-lmpi_mpifh" "-lmpi_usempif08" "-lmpi_usempi_ignore_tkr" )
endif()
