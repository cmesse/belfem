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
        test_MetalElastic.cpp
        test_MaterialDerivatives.cpp
        test_MaterialFactoryDeck.cpp
        test_MaterialValues.cpp
        test_SplineLookupTable.cpp
        test_Database.cpp
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

# -----------------------------------------------------------------------------
# user-material plugin fixture
# -----------------------------------------------------------------------------
# The suite's only dlopen-loaded artifact. test_JcCustomGuard.cpp needs a real
# plugin to reach the two-dependency jc / n registration guards, because a
# plugin's <label>_init is the only surface that calls them. It is built into
# <build>/tests/physics/plugindata/material/ and the test passes that absolute
# path to UserDefinedMaterial, so no data-path search is involved.
#
# LIBRARY_OUTPUT_DIRECTORY keeps the module out of the global
# ${CMAKE_BINARY_DIR}/lib, which a USE_SHARED_LIBS=ON test binary carries on
# its RPATH.
#
# MODULE, not SHARED: this is only ever dlopen()ed, never linked. It also keeps
# the artifact lib<name>.so on macOS rather than .dylib, matching what a deck
# writes ( same reasoning as UserMaterialTemplate.cmake ).
#
# It links nothing on purpose: the plugin reaches the material only through
# the virtual set_user_defined_function on the object the executable hands it,
# so it needs no BELFEM symbol of its own. The executable is built -rdynamic
# anyway ( config/compiler/config_gcc.cmake, -Wl,-export_dynamic under Apple
# Clang ), which is what a real plugin that calls a non-virtual host symbol
# relies on; this fixture does not exercise that.
set( BELFEM_TEST_PLUGIN_ROOT ${CMAKE_CURRENT_BINARY_DIR}/plugindata )

add_library( test_material_plugin MODULE plugin/test_material_plugin.cpp )

set_target_properties( test_material_plugin PROPERTIES
        PREFIX "lib"
        SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY ${BELFEM_TEST_PLUGIN_ROOT}/material )

# the include closure of cl_Material.hpp, i.e. the list src/physics/materials
# builds with plus that directory itself; the circuit fixture's list would not
# compile a material translation unit
target_include_directories( test_material_plugin PRIVATE
        ${BELFEM_SOURCE_DIR}/core
        ${BELFEM_SOURCE_DIR}/comm
        ${BELFEM_SOURCE_DIR}/containers
        ${BELFEM_SOURCE_DIR}/linalg
        ${BELFEM_SOURCE_DIR}/numerics/spline
        ${BELFEM_SOURCE_DIR}/numerics/bezier
        ${BELFEM_SOURCE_DIR}/numerics/integration
        ${CMAKE_BINARY_DIR}/src/numerics/integration
        ${BELFEM_SOURCE_DIR}/math/tensor
        ${BELFEM_SOURCE_DIR}/math/tools
        ${BELFEM_SOURCE_DIR}/physics
        ${BELFEM_SOURCE_DIR}/physics/materials
        ${BELFEM_SOURCE_DIR}/physics/database
        ${BELFEM_SOURCE_DIR}/mesh
        ${CMAKE_BINARY_DIR}/generated )

if( APPLE )
    # Mach-O rejects undefined symbols at link time; the host resolves them
    target_link_options( test_material_plugin PRIVATE -Wl,-undefined,dynamic_lookup )
endif()

# the test only compiles its cases when it knows where the fixture landed, so
# a tree that somehow skips the plugin target degrades to "not built" rather
# than to a silently passing suite
target_compile_definitions( test_physics
        PRIVATE BELFEM_TEST_PLUGIN_ROOT="${BELFEM_TEST_PLUGIN_ROOT}" )

add_dependencies( test_physics test_material_plugin )
