/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

//------------------------------------------------------------------------------
// the source-function plugin resolves through the shared data path.
//
// This is the discriminating gate for the fix. Before it,
// SourceFunction::read_user_defined() handed the deck's `file :` string
// straight to dlopen, so a bare plugin name that lived only under
// $BELFEM_DATA/material could not be found, while the same name in a
// `materials` section resolved fine. These cases are RED on the pre-fix tree:
// dlopen fails, read_user_defined raises, and the test binary is built to
// throw on error ( test_circuit_main.cpp ).
//
// The plugin is built by tests/circuit/CMakeLists.txt into
// BELFEM_TEST_PLUGIN_ROOT/material/ -- deliberately NOT into the project
// library directory, which a USE_SHARED_LIBS=ON test binary would carry on its
// RPATH. That placement is what makes these cases discriminate rather than
// merely pass, and EnvironmentIsNotContaminated is its control.
//
// Not covered here: the resolver's own search order, which is unit-tested
// against the real share/material tree in tests/io/test_filetools.cpp.
//------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include "typedefs.hpp"
#include "globals.hpp"
#include "filetools.hpp"
#include "cl_SourceFunction.hpp"

using belfem::real;
using belfem::SourceFunction;
using belfem::SourceFunctionType;

#ifdef BELFEM_TEST_PLUGIN_ROOT

namespace
{
    //! the fixture library, named the way a deck would name it: no directory
    //! component at all, so dlopen cannot find it without the data-path search
    const std::string tPluginName = "libtest_source_plugin.so" ;

    const std::string tPluginLabel = "TestPluginSource" ;

    //! matches test_source_plugin.cpp
    real
    expected( const real aTime )
    {
        return 7.0 * aTime + 3.0 ;
    }

    //! RAII for the resolver root. The environment cannot be used: Add_Test
    //! forces BELFEM_DATA into every test process, and Communicator::
    //! set_globals() copies it once before any test runs, so the global is the
    //! only thing a test can move
    class ScopedDataPath
    {
        const belfem::string mSaved ;

    public:

        explicit ScopedDataPath( const belfem::string & aPath ) :
                mSaved( belfem::gBelfemDataPath )
        {
            belfem::gBelfemDataPath = aPath ;
        }

        ~ScopedDataPath()
        {
            belfem::gBelfemDataPath = mSaved ;
        }
    };
}

//------------------------------------------------------------------------------

// the fixture must not be reachable any other way, or the test below would
// pass without the fix and prove nothing
TEST( SourcePluginPath, FixtureIsNotReachableFromTheRunDirectory )
{
    EXPECT_FALSE( belfem::file_exists( tPluginName ) );

    ScopedDataPath tRoot( BELFEM_TEST_PLUGIN_ROOT );

    EXPECT_TRUE( belfem::file_exists(
            belfem::string( BELFEM_TEST_PLUGIN_ROOT ) + "/material/" + tPluginName ) );
}

//------------------------------------------------------------------------------

// RED before the fix. Precisely: dlopen() of a slash-free name searches the
// caller's DT_RPATH/DT_RUNPATH, then the loader environment path, then
// ld.so.cache and the system directories ( and the Darwin equivalents ) --
// NOT just $LD_LIBRARY_PATH, as an earlier version of this comment claimed.
// What keeps the fixture out of all of those is that it is built into
// plugindata/material/ rather than into ${CMAKE_BINARY_DIR}/lib, and carries
// a name nothing else ships. See EnvironmentIsNotContaminated below, which is
// the control for that claim rather than an assertion of it
TEST( SourcePluginPath, BareNameLoadsFromTheDataDirectory )
{
    ScopedDataPath tRoot( BELFEM_TEST_PLUGIN_ROOT );

    SourceFunction tFun ;

    ASSERT_NO_THROW( tFun.read_user_defined( tPluginName, tPluginLabel ) );

    EXPECT_EQ( tFun.type(), SourceFunctionType::UserDefined );

    // the plugin really ran, rather than the object surviving in a default
    // state: this waveform is not one any builtin type produces
    EXPECT_NEAR( tFun.compute( 0.0 ), expected( 0.0 ), 1.0e-12 );
    EXPECT_NEAR( tFun.compute( 2.5 ), expected( 2.5 ), 1.0e-12 );
}

//------------------------------------------------------------------------------

// CONTROL for the case above. With no resolver root, the bare fixture name has
// only the loader's own search left -- so if this LOADS, the fixture is
// reachable some other way ( a polluted loader path, an installed library of
// the same name ), and BareNameLoadsFromTheDataDirectory above would pass
// without the fix and prove nothing. A failure here is a loud environment
// problem, which is the point: it converts a silent false pass into a red test
TEST( SourcePluginPath, EnvironmentIsNotContaminated )
{
    ScopedDataPath tRoot( "" );

    SourceFunction tFun ;

    EXPECT_ANY_THROW( tFun.read_user_defined( tPluginName, tPluginLabel ) )
            << "the fixture loaded with no resolver root, so it is reachable "
               "through the dynamic loader's own search path. The "
               "BareNameLoadsFromTheDataDirectory case cannot discriminate in "
               "this environment." ;
}

//------------------------------------------------------------------------------

// an unresolvable name still reaches dlopen unchanged, so the error names what
// the user wrote rather than a path they never asked for
TEST( SourcePluginPath, UnresolvedNameIsReportedAsWritten )
{
    ScopedDataPath tRoot( BELFEM_TEST_PLUGIN_ROOT );

    SourceFunction tFun ;

    EXPECT_ANY_THROW( tFun.read_user_defined( "no_such_plugin.so", tPluginLabel ) );
}

//------------------------------------------------------------------------------

// a second plugin on one object would strand the first mapping
TEST( SourcePluginPath, SecondLoadIsRefused )
{
    ScopedDataPath tRoot( BELFEM_TEST_PLUGIN_ROOT );

    SourceFunction tFun ;

    ASSERT_NO_THROW( tFun.read_user_defined( tPluginName, tPluginLabel ) );
    EXPECT_ANY_THROW( tFun.read_user_defined( tPluginName, tPluginLabel ) );
}

//------------------------------------------------------------------------------

#endif // BELFEM_TEST_PLUGIN_ROOT
