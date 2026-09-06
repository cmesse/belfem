set( TESTNAME physics )

include_directories( ${BELFEM_SOURCE_DIR}/math/tensor )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${CMAKE_BINARY_DIR}/src/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/database )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )

set( SOURCES
        test_YBCO.cpp
        test_JcFunctions.cpp
        test_MetalFieldDependence.cpp
        test_JcCustomGuard.cpp
        )

set( LIBLIST
        mesh
        integration
        interpolation
        bezier
        database
        materials
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )

# Gas tables and gas models. Their data lives in share/fluid and is found at
# runtime, so these tests are gated at configure time on the modules existing and
# skip themselves at runtime when the tables cannot be located.
if( USE_GASMODELS )
    add_subdirectory( gastables )
    add_subdirectory( gasmodels )
endif()

# compile-only gate: user-material API must build without a linalg backend
add_subdirectory( backendfree )
