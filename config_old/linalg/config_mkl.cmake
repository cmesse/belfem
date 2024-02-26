if( USE_MKL )
    include( ${BELFEM_CONFIG_DIR}/system/find_mkl.cmake )
    # common setup
    list( APPEND BELFEM_DEFS "BELFEM_MKL" )
    # we link against 32 bit indices
    #list( APPEND BELFEM_DEFS "MKL_ILP64" )
    list( APPEND BELFEM_INCLUDES "${BELFEM_MKLROOT}/include" )

    #see  https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

    #if( ${COMPILER_ID} EQUAL 2 ) # Intel
    #    if( APPLE )
    #        set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib" )
    #    elseif( UNIX )
    #        set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${MKLROOT}/lib/intel64 -Wl,-rpath,${MKLROOT}/lib/intel64" )
    #    endif()
    #endif()

    if( APPLE )
        SET( BELFEM_MKL_LIBDIR "${BELFEM_MKLROOT}/lib")
        # core
    elseif( UNIX )
        SET( BELFEM_MKL_LIBDIR  "${BELFEM_MKLROOT}/lib/intel64" )
    else()
        message(FATAL_ERROR "Unknown Operating System" )
    endif()

    #list ( APPEND BELFEM_RPATH "${BELFEM_MKL_LIBDIR}" )
    #set( BELFEM_MATRIX_LIBS "-lmkl_rt -lpthread -lm -ldl" )

    if( USE_PARDISO )
        list( APPEND BELFEM_DEFS "BELFEM_PARDISO" )
    endif()


    # define interface
    if( USE_MKL_64BIT_API )
        set( BELFEM_MATRIX_LIBS "${BELFEM_MKL_LIBDIR}/libmkl_intel_ilp64.a ")
        list( APPEND BELFEM_DEFS "MKL_ILP64" )
        if( NOT BELFEM_USE_INTEL )
            set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -m64" )
        endif()
    else()
        set( BELFEM_MATRIX_LIBS "${BELFEM_MKL_LIBDIR}/libmkl_intel_lp64.a ")
    endif()

    # add thread
    if( USE_OPENMP )
        if( ${COMPILER_ID} EQUAL 1 AND NOT BELFEM_USE_CLANG )
            set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_gnu_thread.a ")
        else()
            set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_intel_thread.a ")
        endif()
    else()
        set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_sequential.a ")
    endif()

    # add core
    set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_core.a ")

    # BLACS
    if( USE_MPI )
        if( APPLE ) # assume mpich as default
            if( USE_MKL_64BIT_API )
                set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_blacs_mpich_ilp64.a ")
            else()
                set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_blacs_mpich_lp64.a ")
            endif()
        else() # assume openmp
            if( USE_MKL_64BIT_API )
                set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_blacs_openmpi_ilp64.a ")
            else()
                set( BELFEM_MATRIX_LIBS " ${BELFEM_MATRIX_LIBS} ${BELFEM_MKL_LIBDIR}/libmkl_blacs_openmpi_lp64.a ")
            endif()
        endif()
    endif()

    # group for gnu compiler
    if( NOT APPLE )
        set(BELFEM_MATRIX_LIBS "-Wl,--start-group ${BELFEM_MATRIX_LIBS} -Wl,--end-group" )
    endif()

    # add scalapack
    if( USE_MKL_64BIT_API )
        set( BELFEM_MATRIX_LIBS "${BELFEM_MKL_LIBDIR}/libmkl_scalapack_ilp64.a ${BELFEM_MATRIX_LIBS}" )
    else()
        set( BELFEM_MATRIX_LIBS "${BELFEM_MKL_LIBDIR}/libmkl_scalapack_lp64.a ${BELFEM_MATRIX_LIBS}" )
    endif()

    if( USE_OPENMP )
        if( ${COMPILER_ID} EQUAL 1 )
            set( BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -lgomp" )
        else()
            set( BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -liomp5" )
        endif()
    endif()
    #message( "MKL FLAGS: ${BELFEM_MATRIX_LIBS}")
    set(BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -lpthread -lm -ldl" )
else()
    list( APPEND BELFEM_DEFS "BELFEM_NETLIB" )

    set( BELFEM_MATRIX_LIBS " -llapack -lfspblas -lcblas -lblas")

    if( NOT APPLE )
        set(BELFEM_MATRIX_LIBS "-Wl,--start-group ${BELFEM_MATRIX_LIBS} -Wl,--end-group" )
    endif()
    if( USE_MPI )
        set( BELFEM_MATRIX_LIBS "-lscalapack ${BELFEM_MATRIX_LIBS}")
    endif()

    if( USE_OPENMP )
        if( ${COMPILER_ID} EQUAL 1 )
            set( APPEND BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -lgomp" )
        else()
            set( APPEND BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -liomp5" )
        endif()
    endif()

    set( BELFEM_MATRIX_LIBS "${BELFEM_MATRIX_LIBS} -lpthread -ldl" )

    if( USE_PARDISO )
        message(FATAL_ERROR "Turn on MKL if you want to use Pardiso" )
    endif()
endif()