if( USE_MPI )
    # test if MPI_HOME is set
    if( DEFINED ENV{MPI_HOME} )
        # for summary output
        set( BELFEM_MPIHOME $ENV{MPI_HOME} )
        if( NOT $ENV{MPI_HOME} STREQUAL $ENV{TPLIBS} )
            # step 3: find library directory
            if( IS_DIRECTORY $ENV{MPI_HOME}/lib64 )
                list( APPEND BELFEM_RPATH $ENV{MPI_HOME}/lib64 )
            elseif( IS_DIRECTORY $ENV{MPI_HOME}/lib )
                list( APPEND BELFEM_RPATH $ENV{MPI_HOME}/lib )
            else()
                message( FATAL_ERROR "Could not find MPI library directory" )
            endif()
        endif()
    endif()
endif()