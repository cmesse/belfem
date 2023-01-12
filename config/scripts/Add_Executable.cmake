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
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )

# add default libraries
add_executable( ${EXECNAME} ${MAIN})

# add basellibs
# link executable to libs
foreach( LIBITEM ${BELFEM_LIBLIST_BASE} )
     add_dependencies( ${EXECNAME} lib${LIBPREFIX}_${LIBITEM}.${LIBSUFFIX} )
endforeach()

# link executable to libs
foreach( LIBITEM ${LIBLIST} )
    add_dependencies( ${EXECNAME} lib${LIBPREFIX}_${LIBITEM}.${LIBSUFFIX} )
endforeach()

foreach( LIBITEM ${BELFEM_LIBLIST_BASE} )
    set( BELFEM_LIBS "-l${LIBPREFIX}_${LIBITEM} ${BELFEM_LIBS}" )
endforeach()

foreach( LIBITEM ${LIBLIST} )
    set( BELFEM_LIBS "-l${LIBPREFIX}_${LIBITEM} ${BELFEM_LIBS}" )
endforeach()


# tidy up
string(STRIP "${BELFEM_LIBS}" BELFEM_LIBS)

#target_include_directories( ${EXECNAME} PRIVATE $ENV{TPLS}/include )
if ( USE_PROFILER )
    target_link_libraries( ${EXECNAME} ${BELFEM_LIBS} ${BELFEM_SOLVER_LIBS} ${BELFEM_MATRIX_LIBS} ${BELFEM_IO_LIBS} ${BELFEM_PERFORMANCELIBS} ${BELFEM_FORTRANLIBS} ${BELFEM_OPENMPLIBS} )
else()
    target_link_libraries( ${EXECNAME} ${BELFEM_LIBS} ${BELFEM_SOLVER_LIBS} ${BELFEM_MATRIX_LIBS} ${BELFEM_IO_LIBS} ${BELFEM_FORTRANLIBS} ${BELFEM_OPENMPLIBS} )
endif()


if( APPLE )
    add_custom_command(TARGET ${EXECNAME}
            POST_BUILD
            COMMAND codesign -f --deep -s "Christian Messe" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXECNAME}
            )
endif()
