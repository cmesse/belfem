# check if we want to use MKL in the first place
if( USE_MKL )
    # test if MKL environent variable was set
    if( DEFINED ENV{MKLROOT} )
        set( BELFEM_MKLROOT $ENV{MKLROOT} )
    else()
        # take a reasonable guess
        set( BELFEM_MKLROOT "/opt/intel/oneapi/mkl/latest" )
    endif()

    # test if directory exists
    if( NOT EXISTS ${BELFEM_MKLROOT} )
        message( FATAL_ERROR, "MKL Home Directory could not be found" )
    endif()
endif()