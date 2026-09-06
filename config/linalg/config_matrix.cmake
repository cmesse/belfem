# -------------------------------------------------------------------------
# Flags for MATRIX Library
# -------------------------------------------------------------------------

if ( USE_MATRIX_ARMADILLO AND USE_MATRIX_BLAZE )
    message(FATAL_ERROR "Can only set Armadillo or Blaze; not both" )
elseif( NOT ( USE_MATRIX_ARMADILLO OR USE_MATRIX_BLAZE ) )
    message(FATAL_ERROR "Select a matrix backend: Armadillo or Blaze" )
endif()

if ( USE_MATRIX_ARMADILLO )
    list( APPEND BELFEM_DEFS "BELFEM_ARMADILLO" )
    # keep -lsuperlu and -larpack here: only the SCLS Armadillo build is
    # guaranteed to carry these dependencies itself
    belfem_link_libraries( 7 "-larmadillo" "-lsuperlu" "-larpack" )
    # special flags for intel
    if ( ${COMPILER_ID} EQUAL 2 )
        list( APPEND BELFEM_DEFS "ARMA_ALLOW_FAKE_GCC" )
        list( APPEND BELFEM_DEFS "ARMA_ALLOW_FAKE_CLANG" )
    endif()
else()
    list( APPEND BELFEM_DEFS "BELFEM_BLAZE" )
endif()