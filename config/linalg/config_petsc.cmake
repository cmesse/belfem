if( USE_PETSC )
    if( NOT USE_MPI )
        message(FATAL_ERROR "Turn on MPI if you want to link against PETSc" )
    endif()

    list( APPEND BELFEM_DEFS "BELFEM_PETSC" )

    if( DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "" )
        # classic in-tree build: headers live in both the source and the arch
        # include directories, the library under the arch directory
        if( NOT DEFINED ENV{PETSC_DIR} )
            message(FATAL_ERROR "PETSC_ARCH is set but PETSC_DIR is not" )
        endif()
        list( APPEND BELFEM_INCLUDES "$ENV{PETSC_DIR}/include" )
        list( APPEND BELFEM_INCLUDES "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include" )
        list( APPEND BELFEM_RPATH "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib" )
        belfem_link_libraries( 9 "-lpetsc" )
    else()
        # prefix install, e.g. SCLS
        belfem_find_package(
                PETSC
                REQUIRED
                HEADERS petsc.h
                LIBRARIES petsc
                INCLUDE_SUFFIXES . petsc )

        list( APPEND BELFEM_INCLUDES "${PETSC_INCLUDE_DIR}" )
        list( APPEND BELFEM_RPATH "${PETSC_LIB_DIR}" )
        belfem_link_libraries( 9 ${PETSC_LIBRARIES} )
    endif()

    # only needed if PETSc was built with X support
    #belfem_link_libraries( 0 "-lX11" "-lmpfr" )
endif()
