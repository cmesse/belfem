if( USE_STRUMPACK OR USE_STRUMPACK_CUDA )
    if ( NOT USE_MPI )
        message(FATAL_ERROR "Turn on MPI if you want to link against STRUMPACK" )
    endif()
    if ( NOT USE_OPENMP )
        message(FATAL_ERROR "Turn on OpenMP if you want to link against STRUMPACK" )
    endif()
    if ( NOT USE_SCOTCH )
        message(FATAL_ERROR "Turn on SCOTCH if you want to link against STRUMPACK" )
    endif()

    set( STRUMPACK_LIBRARY_NAMES
            strumpack
            sbutterflypack
            dbutterflypack
            cbutterflypack
            zbutterflypack
            slate_lapack_api
            slate )

    if( USE_STRUMPACK_ZFP )
        list( APPEND STRUMPACK_LIBRARY_NAMES zfp )
    endif()

    belfem_find_package(
            STRUMPACK
            REQUIRED
            HEADERS StrumpackSparseSolver.hpp
            LIBRARIES ${STRUMPACK_LIBRARY_NAMES} )

    list( APPEND BELFEM_DEFS "BELFEM_STRUMPACK" )
    list( APPEND BELFEM_INCLUDES "${STRUMPACK_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${STRUMPACK_LIB_DIR}" )

    # the stack (strumpack, butterflypack, slate, zfp) is registered as one
    # block: internal order already descends the dependency graph, and its
    # only consumer is petsc (rank 9)
    belfem_link_libraries( 8 ${STRUMPACK_LIBRARIES} )

    if( USE_STRUMPACK_CUDA )
        belfem_link_libraries( 1 "-L$ENV{CUDA_DIR}/cuda/lib64" "-L$ENV{CUDA_DIR}/comm_libs/nvshmem/lib" "-L$ENV{CUDA_DIR}/math_libs/lib64" "-lcusolver_static" "-lcusparse_static" "-lcublas_static" "-lcublasLt_static" "-lnvshmem" "-lcudadevrt" "-lcudart_static" "-lculibos" "-lrt" )
    endif()
endif()
