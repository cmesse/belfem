if( USE_PETSC )

    if ( USE_MPI )
    else()
        message(FATAL_ERROR "Turn of MPI if you want to link against PETSC" )
    endif()

    # test if PETSC environent variable was set
    if( DEFINED ENV{PETSC_DIR} )

        if( DEFINED ENV{PETSC_ARCH} )
            list( APPEND BELFEM_INCLUDES "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include" )
            list( APPEND BELFEM_INCLUDES "$ENV{PETSC_DIR}/include" )
            list( APPEND BELFEM_RPATH $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib )
            if( APPLE )
                list( APPEND BELFEM_RPATH /opt/X11/lib )
            endif()
        else()
            message(FATAL_ERROR "Please set the PETSC_ARCH environment variable." )
        endif()

    else()
        message(FATAL_ERROR "Please set the PETSC_DIR environment variable." )
    endif()

    list( APPEND BELFEM_DEFS "BELFEM_PETSC" )

    if( BELFEM_USE_INTEL )
        set( BELFEM_MATRIX_LIBS "-lpetsc ${BELFEM_MATRIX_LIBS} -lX11")
    else()
        set( BELFEM_MATRIX_LIBS "-lpetsc ${BELFEM_MATRIX_LIBS} -lX11 -lmpfr")
    endif()
endif()