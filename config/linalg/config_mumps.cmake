if( USE_MUMPS )
    if ( USE_MPI )
    else()
        message(FATAL_ERROR "Turn on MPI if you want to link against MUMPS" )
    endif()

    list( APPEND BELFEM_DEFS "BELFEM_MUMPS" )
    set( BELFEM_MUMPS_LIBS "-ldmumps -lmumps_common -lesmumps -lpord -lstrumpack -lparmetis -lptscotch -lptscotcherr -lscotch -lscotcherr")


    list( APPEND BELFEM_INCLUDES "$ENV{TPLS}/include/mumps" )
    set( BELFEM_FCFLAGS   "${BELFEM_FCFLAGS} -I$ENV{TPLS}/include/mumps" )


    set( BELFEM_MATRIX_LIBS "${BELFEM_MUMPS_LIBS} ${BELFEM_MATRIX_LIBS}")


    #if( BELFEM_USE_INTEL ) <-- use this if mpiicc (intel mpi)
    #else()
        # for mpich, use -lmpifort instead of -lmpi_mpifh
    set( BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -lmpi_mpifh -lmpi_usempif08 -lmpi_usempi_ignore_tkr" )
    #endif()

endif()