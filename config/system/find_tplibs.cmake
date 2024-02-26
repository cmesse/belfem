# test if environment variable is set
if( DEFINED ENV{SCLS} )
    # library directory
    if( IS_DIRECTORY $ENV{SCLS}/lib64 )
        list( APPEND BELFEM_RPATH $ENV{SCLS}/lib64 )
        set( SCLSLIBDIR $ENV{SCLS}/lib64 )
        link_directories($ENV{SCLS}/lib64)
    elseif( IS_DIRECTORY $ENV{SCLS}/lib )
        list( APPEND BELFEM_RPATH $ENV{SCLS}/lib )
        set( SCLSLIBDIR $ENV{SCLS}/lib )
        link_directories($ENV{SCLS}/lib)
    else()
        message( FATAL_ERROR "Could not find SCLS library directory" )
    endif()

    # include directory
    if( IS_DIRECTORY $ENV{SCLS}/include )
        list( APPEND BELFEM_INCLUDES $ENV{SCLS}/include )
    else()
        message( FATAL_ERROR "Could not find SCLS include directory" )
    endif()

else()
    # take a guess
    if( IS_DIRECTORY /opt/tpls/sw_gcc/include )
        list( APPEND BELFEM_INCLUDES /opt/tpls/sw_gcc/include )
    endif()

    # library directory
    if( IS_DIRECTORY /opt/tpls/sw_gcc/lib64 )
        list( APPEND BELFEM_RPATH /opt/tpls/sw_gcc/lib64 )
        link_directories(/opt/tpls/sw_gcc/lib64)
    elseif( IS_DIRECTORY /opt/tpls/sw_gcc/lib )
        list( APPEND BELFEM_RPATH /opt/tpls/sw_gcc/lib )
        link_directories(/opt/tpls/sw_gcc/lib)
    else()
        message( FATAL_ERROR "SCLS environment variable is not set" )
    endif()
endif()