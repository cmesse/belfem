# test if environment variable is set
if( DEFINED ENV{TPLS} )
    # library directory
    if( IS_DIRECTORY $ENV{TPLS}/lib64 )
        list( APPEND BELFEM_RPATH $ENV{TPLS}/lib64 )
        set( TPLSLIBDIR $ENV{TPLS}/lib64 )
        link_directories($ENV{TPLS}/lib64)
    elseif( IS_DIRECTORY $ENV{TPLS}/lib )
        list( APPEND BELFEM_RPATH $ENV{TPLS}/lib )
        set( TPLSLIBDIR $ENV{TPLS}/lib )
        link_directories($ENV{TPLS}/lib)
    else()
        message( FATAL_ERROR "Could not find TPLS library directory" )
    endif()

    # include directory
    if( IS_DIRECTORY $ENV{TPLS}/include )
        list( APPEND BELFEM_INCLUDES $ENV{TPLS}/include )
    else()
        message( FATAL_ERROR "Could not find TPLS include directory" )
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
        message( FATAL_ERROR "TPLS environment variable is not set" )
    endif()
endif()