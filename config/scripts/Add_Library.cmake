include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# add default includes
include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/linalg )
include_directories( ${BELFEM_SOURCE_DIR}/linalg/lapack )
if ( USE_MATRIX_ARMADILLO )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/armadillo )
elseif( USE_MATRIX_BLAZE )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/blaze )
endif()
include_directories( ${BELFEM_SOURCE_DIR}/linalg/operators )
include_directories( ${BELFEM_SOURCE_DIR}/io )

include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/sparse )

include( ${BELFEM_CONFIG_DIR}/globals.cmake )


# test date
add_library( lib${LIBPREFIX}_${LIBNAME}.${LIBSUFFIX} STATIC ${SOURCES} )

#target_link_libraries( lib${LIBPREFIX}_${LIBNAME}.${LIBSUFFIX} ${BELFEM_TPLIBS} )

# Set the output path for library
set_target_properties( lib${LIBPREFIX}_${LIBNAME}.${LIBSUFFIX} PROPERTIES OUTPUT_NAME ${LIBPREFIX}_${LIBNAME} )
