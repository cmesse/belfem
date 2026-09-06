if( USE_MKL )
    if( APPLE )
        message(FATAL_ERROR "MKL is not supported on macOS anymore. Use the netlib backend instead." )
    endif()

    include( ${BELFEM_CONFIG_DIR}/system/find_mkl.cmake )
    # common setup
    list( APPEND BELFEM_DEFS "BELFEM_MKL" )
    list( APPEND BELFEM_INCLUDES "${BELFEM_MKLROOT}/include" )

    # see https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
    if( IS_DIRECTORY "${BELFEM_MKLROOT}/lib/intel64" )
        set( BELFEM_MKL_LIBDIR "${BELFEM_MKLROOT}/lib/intel64" )
    else()
        # oneAPI 2024+ dropped the intel64 subdirectory
        set( BELFEM_MKL_LIBDIR "${BELFEM_MKLROOT}/lib" )
    endif()

    # link dynamically; finalize_compiler.cmake turns this into -L and rpath
    list( APPEND BELFEM_RPATH "${BELFEM_MKL_LIBDIR}" )

    if( USE_PARDISO )
        list( APPEND BELFEM_DEFS "BELFEM_PARDISO" )
    endif()

    # interface layer
    if( USE_MKL_64BIT_API )
        set( BELFEM_MKL_SUFFIX "ilp64" )
        list( APPEND BELFEM_DEFS "MKL_ILP64" )
        # list( APPEND BELFEM_DEFS "BELFEM_INT64" ) already set in main CMakeLists.txt
        if( NOT BELFEM_USE_INTEL )
            set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -m64" )
        endif()
    else()
        set( BELFEM_MKL_SUFFIX "lp64" )
    endif()

    # the mkl block is atomic: internal order follows the intel link
    # advisor; --no-as-needed is needed on linkers that default
    # to --as-needed (Debian)
    set( BELFEM_MKL_LIBS
            "-Wl,--no-as-needed"
            "-lmkl_scalapack_${BELFEM_MKL_SUFFIX}"
            "-lmkl_intel_${BELFEM_MKL_SUFFIX}" )

    # threading layer
    if( USE_OPENMP )
        if( ${COMPILER_ID} EQUAL 1 AND NOT BELFEM_USE_CLANG )
            list( APPEND BELFEM_MKL_LIBS "-lmkl_gnu_thread" )
        else()
            list( APPEND BELFEM_MKL_LIBS "-lmkl_intel_thread" )
        endif()
    else()
        list( APPEND BELFEM_MKL_LIBS "-lmkl_sequential" )
    endif()

    # core
    list( APPEND BELFEM_MKL_LIBS "-lmkl_core" )

    # BLACS, assume openmpi
    if( USE_MPI )
        list( APPEND BELFEM_MKL_LIBS "-lmkl_blacs_openmpi_${BELFEM_MKL_SUFFIX}" )
    endif()

    belfem_link_libraries( 1 ${BELFEM_MKL_LIBS} )
else()
    list( APPEND BELFEM_DEFS "BELFEM_NETLIB" )

    # scalapack sits above the mpi layer in the dependency graph
    if( USE_MPI )
        belfem_link_libraries( 5 "-lscalapack" )
    endif()

    # No -lcblas here. CBLAS is a separate component of the netlib distribution and
    # a reference BLAS/LAPACK prefix built without it ships no libcblas, so the flag
    # broke the link outright on such a prefix. It never provided a symbol either:
    # Armadillo and Blaze go to the Fortran interface on this path, and where SLATE
    # (via STRUMPACK) does want cblas_*gemm_batch the provider is MKL or OpenBLAS,
    # never netlib CBLAS, which has no batched GEMM at all.
    set( BELFEM_BLAS_LIBS "-llapack" "-lblas")

    if( NOT APPLE )
        list(INSERT BELFEM_BLAS_LIBS 0 "-Wl,--start-group" )
        list(APPEND BELFEM_BLAS_LIBS "-Wl,--end-group" )
    endif()

    belfem_link_libraries( 1 ${BELFEM_BLAS_LIBS} )

    if( USE_PARDISO )
        message(FATAL_ERROR "Turn on MKL if you want to use Pardiso" )
    endif()
endif()

# runtime and system libraries come last
#
# Note this block runs even when USE_MKL is OFF: CMakeLists.txt includes this file
# unconditionally. For the GNU compilers the OpenMP runtime already arrives through
# -fopenmp on the link line, so naming it here as well produced
#     ld: warning: ignoring duplicate libraries: '-lgomp'
# The Intel branch is left as it was; -qopenmp handling there is untested.
if( USE_OPENMP )
    if( ${COMPILER_ID} EQUAL 1 )
        if( NOT CMAKE_CXX_FLAGS MATCHES "-fopenmp" AND NOT BELFEM_CXXFLAGS MATCHES "-fopenmp" )
            belfem_link_libraries( 0 "-lgomp" )
        endif()
    else()
        belfem_link_libraries( 0 "-liomp5" )
    endif()
endif()
belfem_link_libraries( 0 "-lpthread" "-lm" "-ldl" )
