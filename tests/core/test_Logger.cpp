/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the file mode of belfem::Logger: what reaches the log file
 * byte for byte, and which messages never reach it at all.
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "typedefs.hpp"
#include "cl_Logger.hpp"

namespace
{
    //! owns a scratch log file and removes it again. The path is absolute and
    //! rooted in the gtest temporary directory, so the test does not depend on
    //! the working directory the suite is started from; the pid keeps two
    //! concurrent runs of the suite apart.
    class ScopedLogFile
    {
        const std::string mPath ;

    public:

        explicit ScopedLogFile( const std::string & aName ) :
                mPath( ::testing::TempDir() + "belfem_test_logger_" + aName
                       + "_" + std::to_string( ( long ) getpid() ) + ".log" )
        {
        }

        ~ScopedLogFile()
        {
            std::remove( mPath.c_str() );
        }

        const std::string &
        path() const
        {
            return mPath ;
        }
    };

    //! reads the file back in full. The Logger closes the stream in its
    //! destructor, so it must have gone out of scope before this is called.
    //! A file that cannot be opened returns a marker no case expects, so a
    //! missing file is never mistaken for an empty one.
    std::string
    read_log( const std::string & aPath )
    {
        std::ifstream tFile( aPath.c_str(), std::ios::binary );
        if ( ! tFile.is_open() )
        {
            return "<log file could not be opened>" ;
        }
        std::ostringstream tBuffer ;
        tBuffer << tFile.rdbuf();
        return tBuffer.str();
    }
}

// =============================================================================
// what reaches the file
// =============================================================================

TEST( Logger, PercentSequencesSurviveTheStringRoute )
{
    // Logger::message is a printf-style API, so the only route that carries
    // arbitrary text through unchanged is passing it as an argument to "%s".
    // The message is then written with fprintf( stream, "%s", text ), never
    // with the formatted text as the format itself.
    // Mutation: write the formatted text as the format and this case goes red,
    // because the "%s" and the "%d" below are re-interpreted.
    const ScopedLogFile tFile( "percent" );

    {
        belfem::Logger tLog( belfem::InfoLevel::Default, tFile.path() );
        tLog.message( belfem::InfoLevel::Minimal, "%s", "50% of %s and %d" );
    }

    EXPECT_EQ( read_log( tFile.path() ), std::string( "50% of %s and %d\n" ) );
}

TEST( Logger, FormatRouteExpandsArguments )
{
    // the complementary contract: given a format and arguments, the logger
    // expands them, and "%%" is the way to a literal percent sign
    const ScopedLogFile tFile( "format" );

    {
        belfem::Logger tLog( belfem::InfoLevel::Default, tFile.path() );
        tLog.message( belfem::InfoLevel::Minimal, "%d%%", 50 );
    }

    EXPECT_EQ( read_log( tFile.path() ), std::string( "50%\n" ) );
}

// =============================================================================
// what does not reach the file
// =============================================================================

TEST( Logger, LevelFilterWritesNothing )
{
    // a message whose level is above the logger's own level is dropped before
    // it is formatted, so the file stays empty rather than merely unflushed.
    // Mutation: drop the comparison in Logger::message and this goes red.
    const ScopedLogFile tFile( "level" );

    {
        belfem::Logger tLog( belfem::InfoLevel::Default, tFile.path() );
        tLog.message( belfem::InfoLevel::Everything, "%s", "must not appear" );
    }

    EXPECT_TRUE( read_log( tFile.path() ).empty() );
}

TEST( Logger, MessageAtTheLoggersOwnLevelIsWritten )
{
    // the filter is inclusive: a message at exactly the logger's level passes.
    // Mutation: make the comparison strict and this goes red.
    const ScopedLogFile tFile( "equal" );

    {
        belfem::Logger tLog( belfem::InfoLevel::Default, tFile.path() );
        tLog.message( belfem::InfoLevel::Default, "%s", "at the boundary" );
    }

    EXPECT_EQ( read_log( tFile.path() ), std::string( "at the boundary\n" ) );
}
