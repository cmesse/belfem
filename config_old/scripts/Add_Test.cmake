include( ${BELFEM_CONFIG_DIR}/globals.cmake )

include_directories( ${SSF_SRC_DIR}/core )
include_directories( ${SSF_SRC_DIR}/comm )
include_directories( ${SSF_SRC_DIR}/containers )
include_directories( ${SSF_SRC_DIR}/linalg )
include_directories( ${BELFEM_SOURCE_DIR}/linalg/lapack )
if ( USE_MATRIX_ARMADILLO )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/armadillo )
elseif( USE_MATRIX_BLAZE )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/blaze )
endif()
include_directories( ${BELFEM_SOURCE_DIR}/linalg/operators )
include_directories( ${BELFEM_SOURCE_DIR}/io )
include_directories( ${SSF_SRC_DIR}/math/graph )
include_directories( ${SSF_SRC_DIR}/sparse )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test )

add_executable( test_${TESTNAME} test_${TESTNAME}_main.cpp ${SOURCES} )

add_test(NAME ${TESTNAME} COMMAND test_${TESTNAME} )


foreach( LIBITEM ${BELFEM_LIBLIST_BASE} )
    add_dependencies( test_${TESTNAME} lib${LIBPREFIX}_${LIBITEM}.${LIBSUFFIX} )
endforeach()

# link executable to libs
foreach( LIBITEM ${LIBLIST} )
    add_dependencies( test_${TESTNAME} lib${LIBPREFIX}_${LIBITEM}.${LIBSUFFIX} )
endforeach()

foreach( LIBITEM ${BELFEM_LIBLIST_BASE} )
    set( BELFEM_LIBS "-l${LIBPREFIX}_${LIBITEM} ${BELFEM_LIBS}" )
endforeach()

foreach( LIBITEM ${LIBLIST} )
    set( BELFEM_LIBS "-l${LIBPREFIX}_${LIBITEM} ${BELFEM_LIBS}" )
endforeach()

string(STRIP "${BELFEM_LIBS}" BELFEM_LIBS)
target_include_directories( test_${TESTNAME} BEFORE INTERFACE $ENV{TPLS}/include )
target_link_libraries(  test_${TESTNAME} ${BELFEM_LIBS} ${BELFEM_IO_LIBS} ${BELFEM_MATRIX_LIBS}  -lgtest -lgtest_main  ${BELFEM_FORTRANLIBS} )
if( NOT APPLE )
    target_link_libraries(test_${TESTNAME} -pthread )
endif()

link_directories( )

if( APPLE )
    add_custom_command(TARGET test_${TESTNAME}
            POST_BUILD
            COMMAND codesign -f --deep -s "Christian Messe" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${TESTNAME}
            )
endif()
