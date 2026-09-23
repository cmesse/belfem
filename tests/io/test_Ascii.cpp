/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * Where Ascii looks for a relative path: in the process working directory,
 * regardless of the value of $PWD. Every reader that derives from Ascii ( InputFile,
 * CsvFile, GmshReader ) inherits this. Two files of the same name, one in
 * the working directory and one elsewhere, tell the two apart.
 */

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <unistd.h>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Ascii.hpp"

namespace
{
    //! sets or unsets an environment variable and restores it on scope exit
    class ScopedEnv
    {
        const std::string mName ;
        bool              mHadValue = false ;
        std::string       mSaved ;

    public:

        //! aValue == nullptr unsets the variable
        ScopedEnv( const std::string & aName, const char * aValue ) :
            mName( aName )
        {
            const char * tOld = std::getenv( mName.c_str() );
            mHadValue = tOld != nullptr ;
            if ( mHadValue ) mSaved = tOld ;

            if ( aValue == nullptr )
            {
                unsetenv( mName.c_str() );
            }
            else
            {
                setenv( mName.c_str(), aValue, 1 );
            }
        }

        ~ScopedEnv()
        {
            if ( mHadValue )
            {
                setenv( mName.c_str(), mSaved.c_str(), 1 );
            }
            else
            {
                unsetenv( mName.c_str() );
            }
        }
    };

    //! writes aContent to aPath and removes the file on scope exit
    class ScopedFile
    {
        const std::string mPath ;

    public:

        ScopedFile( const std::string & aPath, const std::string & aContent ) :
            mPath( aPath )
        {
            std::ofstream tFile( mPath );
            tFile << aContent << "\n" ;
        }

        ~ScopedFile()
        {
            std::remove( mPath.c_str() );
        }
    };

    //! the name both sentinels share; the pid keeps concurrent runs apart
    std::string
    sentinel_name()
    {
        return "belfem_test_ascii_" + std::to_string( ( long ) getpid() ) + ".txt" ;
    }

    //! the scratch directory without a trailing slash
    std::string
    other_directory()
    {
        std::string tDir = ::testing::TempDir() ;
        while ( tDir.size() > 1 && tDir.back() == '/' ) tDir.pop_back() ;
        return tDir ;
    }

    bool
    other_directory_is_the_working_directory()
    {
        std::error_code tError ;
        return std::filesystem::equivalent( std::filesystem::current_path(),
                                            std::filesystem::path( other_directory() ), tError );
    }
}

TEST( Ascii, RelativePathFollowsTheWorkingDirectory )
{
    if ( other_directory_is_the_working_directory() )
    {
        GTEST_SKIP() << "the scratch directory is the working directory" ;
    }

    // PWD names another directory that holds a file of the same name.
    // Mutation: resolve against getenv( "PWD" ) again and this reads "pwd"
    const std::string tName = sentinel_name() ;
    const ScopedFile tHere( tName, "cwd" );
    const ScopedFile tThere( other_directory() + "/" + tName, "pwd" );
    const ScopedEnv  tPwd( "PWD", other_directory().c_str() );

    const belfem::Ascii tFile( tName, belfem::FileMode::OPEN_RDONLY );
    ASSERT_EQ( tFile.length(), 1u );
    EXPECT_EQ( tFile.line( 0 ), "cwd" );
}

TEST( Ascii, UnsetPwdStillResolves )
{
    // a scrubbed environment ( env -i, cron ) has no PWD at all
    const std::string tName = sentinel_name() ;
    const ScopedFile tHere( tName, "cwd" );
    const ScopedEnv  tPwd( "PWD", nullptr );

    const belfem::Ascii tFile( tName, belfem::FileMode::OPEN_RDONLY );
    ASSERT_EQ( tFile.length(), 1u );
    EXPECT_EQ( tFile.line( 0 ), "cwd" );
}

TEST( Ascii, AbsolutePathIgnoresTheWorkingDirectory )
{
    const std::string tName = sentinel_name() ;
    const ScopedFile tThere( other_directory() + "/" + tName, "elsewhere" );
    const ScopedEnv  tPwd( "PWD", "/nonexistent" );

    const belfem::Ascii tFile( other_directory() + "/" + tName, belfem::FileMode::OPEN_RDONLY );
    ASSERT_EQ( tFile.length(), 1u );
    EXPECT_EQ( tFile.line( 0 ), "elsewhere" );
}
