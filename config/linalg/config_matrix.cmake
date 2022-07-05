# -------------------------------------------------------------------------
# Flags for MATRIX Library
# -------------------------------------------------------------------------

if ( USE_MATRIX_ARMADILLO )
    list( APPEND BELFEM_DEFS "BELFEM_ARMADILLO" )
    set( BELFEM_MATRIX_LIBS "-larmadillo -lsuperlu -larpack ${BELFEM_MATRIX_LIBS}")
    # special flags fir intel
    if ( ${COMPILER_ID} EQUAL 2 )
        list( APPEND BELFEM_DEFS "ARMA_ALLOW_FAKE_GCC" )
        list( APPEND BELFEM_DEFS "ARMA_ALLOW_FAKE_CLANG" )
    endif()
endif()

if ( USE_MATRIX_BLAZE )
    list( APPEND BELFEM_DEFS "BELFEM_BLAZE" )
endif()

if ( USE_MATRIX_ARMADILLO AND USE_MATRIX_BLAZE )
    message(FATAL_ERROR "Can only set Armadillo or Blaze; not both" )
endif()