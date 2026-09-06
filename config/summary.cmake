
message( "Configuration Summary:" )

if( APPLE )
    message( "System     : macOS")
else()
    message( "System     : Linux")
endif()

    message( "Compiler   : " ${BELFEM_COMPILER} " " ${BELFEM_COMPILER_VERSION} )


if( USE_DEBUG )
    message( "Mode       : Debug")
else()
    message( "Mode       : Release")
endif()

if( USE_MPI )
    message( "Using MPI  : yes")
    message( "MPIHOME    : " ${BELFEM_MPIHOME} )
else()
    message( "Using MPI  : no")
endif()

message( "SCLS       : " $ENV{SCLS} )

# the flavor is what defaulted USE_MKL and USE_DEBUG, so print it next to them:
# it is the answer to "why is this tree Release?" without reading any cmake
if( BELFEM_SCLS_FLAVOR )
    message( "SCLS flavor: " ${BELFEM_SCLS_FLAVOR} )
endif()

if( USE_MKL )
    message( "LINALG     : MKL" )
else()
    message( "LINALG     : BLAS and LAPACK")
endif()

if( USE_MATRIX_ARMADILLO )
    message( "Matrix     : Armadillo" )
elseif( USE_MATRIX_BLAZE )
    message( "Matrix     : Blaze" )
endif()

if( USE_TEST )
    if( BELFEM_USE_CLANG )
        message( "Tests      : can not be built with clang")
    else()
        message( "Tests      : yes")
    endif()
else()
    message( "Tests      : no")
endif()

