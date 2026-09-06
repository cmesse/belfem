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
// belfem::search_data_file() -- the shared data-file search order.
//
// These are REGRESSION cases, not the discriminating gate. They cannot run
// red on a
// pre-fix tree, because the function under test does not exist there; a test
// calling it fails to compile rather than to assert. The discriminating gate
// for the fix lives in tests/circuit/test_SourcePluginPath.cpp, which drives
// the actual consumer.
//
// What these lock down is the search order itself, which used to live in
// material::data_file() and now backs every plugin and data-file key in the
// input contract:
//
//   1. the path as written, relative to the run directory, or absolute
//   2. the same relative path below <root>/<subdirectory>
//   3. the file name alone below <root>/<subdirectory>
//   4. otherwise the argument, unchanged
//
// The root is gBelfemDataPath. Add_Test.cmake points it at the source tree's
// share/ for every test in the suite, so share/material/bhdata.hdf5 and
// share/fluid/gasdata.inp are the fixtures below -- real shipped files, so a
// case cannot pass against a fixture it invented.
//
// NOTE on the root: it cannot be moved through the environment. Add_Test
// forces BELFEM_DATA into the test process, and Communicator::set_globals()
// copies it into gBelfemDataPath once, before any test body runs. Cases that
// need a different root therefore save, set and restore the global. For the
// same reason the empty-root case below is named for the global and not for
// "$BELFEM_DATA unset" -- on an installed tree an unset environment variable
// does not imply an empty root, because set_globals() then falls back to
// BELFEM_INSTALL_DATADIR.
//------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

#include "typedefs.hpp"
#include "globals.hpp"
#include "filetools.hpp"

namespace
{
    const std::string tSubdir = "material" ;

    //! a real shipped file, so no case passes against an invented fixture
    const std::string tShipped = "bhdata.hdf5" ;

    std::string
    shipped_path()
    {
        return belfem::gBelfemDataPath + "/" + tSubdir + "/" + tShipped ;
    }

    //! see the header note: the environment cannot move the root
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

    //! creates a file in the RUN directory and removes it again, so that the
    //! run-directory-wins case has something local to win with
    class ScopedLocalFile
    {
        const std::string mName ;

    public:

        explicit ScopedLocalFile( const std::string & aName ) :
                mName( aName )
        {
            std::ofstream tFile( mName );
            tFile << "local" ;
            tFile.close();
        }

        ~ScopedLocalFile()
        {
            std::remove( mName.c_str() );
        }
    };
}

//------------------------------------------------------------------------------

// the suite is only meaningful if the shipped fixture is actually there
TEST( SearchDataFile, FixturesArePresent )
{
    // not ASSERT_GT( size(), 0 ): size_t vs int is a sign-compare, and the
    // suite builds with -Werror
    ASSERT_FALSE( belfem::gBelfemDataPath.empty() )
            << "gBelfemDataPath is empty; Add_Test.cmake should have set BELFEM_DATA" ;

    EXPECT_TRUE( belfem::file_exists( shipped_path() ) );
}

//------------------------------------------------------------------------------

// step 1: the run directory wins over the shared database.
//
// The argument is deliberately a RUN-DIRECTORY-RELATIVE name. An absolute path
// would make this case vacuous -- an existing absolute path always satisfies
// step 1, whether or not a twin sits under the subdirectory, so the test would
// no longer be about precedence at all.
TEST( SearchDataFile, RunDirectoryWins )
{
    // same NAME as the shipped file, so both candidates exist at once
    ScopedLocalFile tLocal( tShipped );

    ASSERT_TRUE( belfem::file_exists( tShipped ) );
    ASSERT_TRUE( belfem::file_exists( shipped_path() ) );

    EXPECT_EQ( belfem::search_data_file( tShipped, tSubdir ), tShipped );
}

//------------------------------------------------------------------------------

// step 2: same relative layout below the data directory
TEST( SearchDataFile, RelativePathBelowSubdirectory )
{
    ASSERT_FALSE( belfem::file_exists( tShipped ) );

    EXPECT_EQ( belfem::search_data_file( tShipped, tSubdir ), shipped_path() );
}

//------------------------------------------------------------------------------

// step 3: the name alone, for a path that only existed on the machine the
// input file came from. This is what lets `bhfile : MatData/bhdata.hdf5 ;`
// resolve in a run directory that has no MatData of its own
TEST( SearchDataFile, NameAloneWhenTheDirectoryIsForeign )
{
    const std::string tForeign = "MatData/" + tShipped ;

    ASSERT_FALSE( belfem::file_exists( tForeign ) );
    ASSERT_FALSE( belfem::file_exists(
            belfem::gBelfemDataPath + "/" + tSubdir + "/" + tForeign ) );

    EXPECT_EQ( belfem::search_data_file( tForeign, tSubdir ), shipped_path() );
}

//------------------------------------------------------------------------------

// step 4: a miss hands the argument back untouched. This is load-bearing --
// the opener reports the name the user wrote, and dlopen keeps the platform
// loader search for a name carrying no slash
TEST( SearchDataFile, MissReturnsTheArgumentUnchanged )
{
    const std::string tAbsent = "no_such_file_anywhere.hdf5" ;

    EXPECT_EQ( belfem::search_data_file( tAbsent, tSubdir ), tAbsent );

    const std::string tAbsentPath = "no/such/directory/file.hdf5" ;

    EXPECT_EQ( belfem::search_data_file( tAbsentPath, tSubdir ), tAbsentPath );
}

//------------------------------------------------------------------------------

// an empty root skips both fallbacks, leaving step 1 as the only lookup.
// Named for the global, not for the environment variable -- see the header
TEST( SearchDataFile, EmptyRootLeavesOnlyTheRunDirectory )
{
    ScopedDataPath tRoot( "" );

    EXPECT_EQ( belfem::search_data_file( tShipped, tSubdir ), tShipped );

    // step 1 still works with no root at all
    ScopedLocalFile tLocal( tShipped );

    EXPECT_EQ( belfem::search_data_file( tShipped, tSubdir ), tShipped );
}

//------------------------------------------------------------------------------

// an empty subdirectory searches the data directory itself, with no stray
// separator. share/ holds no files of its own, so the fixture is one level
// down and reached through step 2
TEST( SearchDataFile, EmptySubdirectorySearchesTheRoot )
{
    const std::string tFluid = "fluid/gasdata.inp" ;

    ASSERT_FALSE( belfem::file_exists( tFluid ) );

    EXPECT_EQ( belfem::search_data_file( tFluid, "" ),
               belfem::gBelfemDataPath + "/" + tFluid );
}

//------------------------------------------------------------------------------

// The subdirectory carries no separator of its own. A leading or trailing
// slash on either side used to be easy to double up; these pin the join
TEST( SearchDataFile, JoinProducesNoDoubleSeparator )
{
    const std::string tFound = belfem::search_data_file( tShipped, tSubdir );

    EXPECT_EQ( tFound.find( "//" ), std::string::npos ) << tFound ;
    EXPECT_EQ( tFound, shipped_path() );
}

//------------------------------------------------------------------------------

// PINS EXISTING BEHAVIOUR, does not endorse it. file_exists() is
// std::filesystem::exists(), which is true for DIRECTORIES, so an empty name
// resolves to the subdirectory itself rather than missing. This case exists so
// that a future change to that behaviour is a deliberate decision with a
// failing test behind it, rather than an accident of a refactor
TEST( SearchDataFile, EmptyNameResolvesToTheDirectoryItself )
{
    const std::string tDir = belfem::gBelfemDataPath + "/" + tSubdir + "/" ;

    ASSERT_TRUE( belfem::file_exists( tDir ) );

    EXPECT_EQ( belfem::search_data_file( "", tSubdir ), tDir );
}

//------------------------------------------------------------------------------

// The remaining quirks from the pre-move implementation. Like the empty-name
// case above these PIN behaviour rather than endorse it
TEST( SearchDataFile, TrailingSlashMissResolvesToTheDirectory )
{
    const std::string tDir = belfem::gBelfemDataPath + "/" + tSubdir + "/" ;

    // "foo/" misses locally and below the subdirectory; step 3 then takes the
    // basename after the last slash, which is empty, and lands on the
    // directory itself
    EXPECT_EQ( belfem::search_data_file( "no_such_dir/", tSubdir ), tDir );
}

//------------------------------------------------------------------------------

TEST( SearchDataFile, MissingAbsolutePathStillFallsBackToTheName )
{
    // an absolute path that does not exist is concatenated anyway ( leaving a
    // // mid-path, which POSIX collapses ), and step 3 then finds the name
    const std::string tAbsent = "/nowhere/at/all/" + tShipped ;

    EXPECT_EQ( belfem::search_data_file( tAbsent, tSubdir ), shipped_path() );
}

//------------------------------------------------------------------------------

// NAMED FOR WHAT IT ACTUALLY TESTS. An earlier version of this case was called
// MaterialDataFileContractIsUnchanged, which both code auditors flagged: it
// never calls material::data_file(), so it cannot pin the forwarder. It pins
// the search order the forwarder delegates to, which is a different claim.
// tests/io cannot include physics/materials, so the forwarder itself is
// covered by inspection only -- one line, fn_material_data_path.cpp:44
TEST( SearchDataFile, SearchOrderBehindTheMaterialForwarder )
{
    EXPECT_EQ( belfem::search_data_file( tShipped, tSubdir ), shipped_path() );

    EXPECT_EQ( belfem::search_data_file( "MatData/" + tShipped, tSubdir ),
               shipped_path() );

    EXPECT_EQ( belfem::search_data_file( "nope.hdf5", tSubdir ), "nope.hdf5" );
}

//------------------------------------------------------------------------------
