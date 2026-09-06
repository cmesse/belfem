# set the name of the test
set( TESTNAME circuit )
set( TESTLABELS fast )

# list the test sources
set( SOURCES
        test_SpiceNumber.cpp
        test_NetlistParser.cpp
        test_NgspiceFactory.cpp
        test_HybridFactory.cpp
        test_ElectricalCircuit.cpp
        test_CircuitAdjacency.cpp
        test_SourceFunction.cpp
        test_SourcePluginPath.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/circuit )

# cl_ElectricalCircuit.hpp reaches into these ( same list as
# src/circuit/CMakeLists.txt )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )

# libbelfem_circuit and its BDF / SourceFunction dependencies are not in
# BELFEM_LIBLIST_BASE; the hybrid-factory tests additionally pull
# PhysicalBoundaryCondition, which drags the fem kernel chain
# ( same list as tests/fem )
set( LIBLIST
        circuit
        ode
        sources
        mesh
        integration
        interpolation
        kernel
        iwg
        postproc
        )

if( USE_MAXWELL )
    list( APPEND LIBLIST maxwell )
endif()

list( APPEND LIBLIST materials )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )

# -----------------------------------------------------------------------------
# source-plugin fixture
# -----------------------------------------------------------------------------
# The suite's only dlopen-loaded artifact. test_SourcePluginPath.cpp needs a
# library that is reachable ONLY through the shared data directory, so this is
# built into <build>/tests/circuit/plugindata/material/ -- not beside the test
# binary, not on the loader path. The test points gBelfemDataPath at that
# plugindata/ and asks for the bare file name; before the search-path fix that
# name reached dlopen unresolved.
#
# LIBRARY_OUTPUT_DIRECTORY below is LOAD-BEARING, not tidiness: without it the
# module would land in the global ${CMAKE_BINARY_DIR}/lib, which a
# USE_SHARED_LIBS=ON test binary carries on its RPATH -- and the bare name
# would then load WITHOUT the fix, so the gate would prove nothing.
#
# MODULE, not SHARED: this is only ever dlopen()ed, never linked. It also keeps
# the artifact lib<name>.so on macOS rather than .dylib, matching what a deck
# writes ( same reasoning as UserMaterialTemplate.cmake ).
#
# It links nothing on purpose: SourceFunction::set_user_defined is resolved out
# of the test executable at load time, which works because the executable is
# built -rdynamic ( config/compiler/config_gcc.cmake:104, and
# -Wl,-export_dynamic under Apple Clang ). That is exactly the contract a real
# plugin author gets from UserLibraryTemplate.cmake, so this fixture also keeps
# that contract under test.
set( BELFEM_TEST_PLUGIN_ROOT ${CMAKE_CURRENT_BINARY_DIR}/plugindata )

add_library( test_source_plugin MODULE plugin/test_source_plugin.cpp )

set_target_properties( test_source_plugin PROPERTIES
        PREFIX "lib"
        SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY ${BELFEM_TEST_PLUGIN_ROOT}/material )

target_include_directories( test_source_plugin PRIVATE
        ${BELFEM_SOURCE_DIR}/core
        ${BELFEM_SOURCE_DIR}/containers
        ${BELFEM_SOURCE_DIR}/linalg
        ${BELFEM_SOURCE_DIR}/numerics/sources
        ${CMAKE_BINARY_DIR}/generated )

if( APPLE )
    # Mach-O rejects undefined symbols at link time; the host resolves them
    target_link_options( test_source_plugin PRIVATE -Wl,-undefined,dynamic_lookup )
endif()

# the test only compiles its cases when it knows where the fixture landed, so
# a tree that somehow skips the plugin target degrades to "not built" rather
# than to a silently passing suite
target_compile_definitions( test_circuit
        PRIVATE BELFEM_TEST_PLUGIN_ROOT="${BELFEM_TEST_PLUGIN_ROOT}" )

add_dependencies( test_circuit test_source_plugin )
