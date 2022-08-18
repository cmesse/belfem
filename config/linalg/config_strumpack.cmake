if( USE_STRUMPACK )
    if ( NOT USE_MPI )
        message(FATAL_ERROR "Turn on MPI if you want to link against STRUMPACK" )
    endif()
    if ( NOT USE_OPENMP )
        message(FATAL_ERROR "Turn on OpenMP if you want to link against STRUMPACK" )
    endif()
    list( APPEND BELFEM_DEFS "BELFEM_STRUMPACK" )
    if( USE_MUMPS )
        set( BELFEM_STRUMPACK_LIBS "-L$ENV{STRUMPACK_DIR}/lib -lstrumpack")
    else()
        set( BELFEM_STRUMPACK_LIBS "-L$ENV{STRUMPACK_DIR}/lib -lstrumpack -lparmetis -lptscotch -lptscotcherr -lscotch -lscotcherr")
    endif()

    list( APPEND BELFEM_INCLUDES "$ENV{STRUMPACK_DIR}/include" )
    set( BELFEM_FCFLAGS   "${BELFEM_FCFLAGS} -I$ENV{STRUMPACK_DIR}/include" )

    list( APPEND BELFEM_RPATH $ENV{STRUMPACK_DIR}/lib )
    link_directories($ENV{STRUMPACK_DIR}/lib)


    set( BELFEM_MATRIX_LIBS "${BELFEM_STRUMPACK_LIBS} ${BELFEM_MATRIX_LIBS}")

    #if( BELFEM_USE_INTEL ) <-- use this if mpiicc (intel mpi)
    #else()
    # for mpich, use -lmpifort instead of -lmpi_mpifh
    #set( BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -lmpi_mpifh -lmpi_usempif08 -lmpi_usempi_ignore_tkr" )
    #endif()

endif()