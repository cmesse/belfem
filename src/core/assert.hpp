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

#ifndef BELFEM_ASSERT_HPP
#define BELFEM_ASSERT_HPP

#include <sstream>
#include <stdexcept>
#include <vector>

#include "typedefs.hpp"
#include "fn_sprint.hpp"

//------------------------------------------------------------------------------

// Whether BELFEM_ASSERT expands to a check at all. Owned by this header: tests
// that exercise an assertion-tier error path must guard on this macro rather
// than restate the condition as #ifndef NDEBUG, which is not equivalent when
// both NDEBUG and DEBUG are defined.
#if !defined( NDEBUG ) || defined( DEBUG )
#define BELFEM_ASSERTIONS_ACTIVE 1
#else
#define BELFEM_ASSERTIONS_ACTIVE 0
#endif

//------------------------------------------------------------------------------

namespace belfem
{
    namespace assert
    {

//------------------------------------------------------------------------------

        std::vector<std::string>
        wrap_lines(std::size_t aMaxWidth, const std::string & aLine);

//------------------------------------------------------------------------------

        void
        hatch_dragon( std::vector< std::string > & aDragon ) ;

//------------------------------------------------------------------------------

        void
        print_line(const std::vector<std::string> & aDragon, std::size_t & aCounter );

//------------------------------------------------------------------------------

        void
        print_line( const std::vector<std::string> & aDragon,
            const string & aLine, std::size_t & aCounter );

//------------------------------------------------------------------------------

        void
        print_line( const string & aLine );

//------------------------------------------------------------------------------

        void
        get_lines( const string & aWhat, std::vector< std::string > & aLines );

//------------------------------------------------------------------------------

        void
        print_errorbox(
            const std::string    & aLocation,
            const std::string    & aTask,
            const std::string    & aCheck,
            const std::vector< std::string  > & aMessage ) ;

//------------------------------------------------------------------------------

        void
        error_abort();

//------------------------------------------------------------------------------

        /**
         * How a failed check reacts: true throws the exception, false calls
         * error_abort(). Initialized to the build's compile-time behaviour --
         * throwing where assertions are active, aborting otherwise -- so a
         * debug run throws and a production run aborts, at ANY rank count,
         * with no caller involvement.
         *
         * DO NOT make the debug reaction depend on the rank count. It is
         * tempting: a BELFEM_ERROR inside an `if ( rank == 0 )` block throws
         * on one rank and leaves the others in their collective, and aborting
         * would end the job cleanly. It was proposed and rejected on
         * 2026-08-30. The parallel debugging workflow is one debugger per
         * rank ( mpirun launching an lldb per process, each in its own
         * terminal ), and a throw is what stops that rank's debugger with a
         * live backtrace while its peers are still inspectable. MPI_Abort
         * tears the whole job down and the developer sees nothing -- which is
         * the failure being debugged, made invisible. Production keeps the
         * abort, because there is no debugger there to serve.
         *
         * Test executables set it to true after gComm.init(), which makes
         * BELFEM_ERROR paths catchable with EXPECT_THROW in a release build as
         * well. It is a test hook, not a user-facing configuration knob:
         * a production run must keep the MPI_Abort reaction, because a throw
         * that escapes main() terminates one rank and leaves its peers blocked
         * in a collective.
         *
         * The state deliberately lives in assert.cpp. error() below is a
         * function template instantiated in the *calling* translation unit, so
         * a static in this header would give each TU its own copy and a test
         * binary could not change the reaction of an already-compiled library.
         */
        bool
        throw_on_error();

        void
        set_throw_on_error( const bool aValue );


//------------------------------------------------------------------------------

        /**
         * Whether a failed check also writes its message to the system log
         * ( syslog, identity "belfem", LOG_CRIT ). Defaults to true so a
         * crashed run leaves a trace in journalctl even when stderr is
         * lost. Deliberately independent of throw_on_error(): a SERIAL debug
         * build still throws -- the daily gdb and make check path -- and
         * gating syslog on the abort branch would silence exactly those runs. Test mains that
         * exercise error paths with EXPECT_THROW disable it alongside
         * set_throw_on_error( true ), so purpose-triggered failures do not
         * spam the log.
         */
        bool
        syslog_on_error();

        void
        set_syslog_on_error( const bool aValue );

        //! write one LOG_CRIT line for the failed check ( rank + location )
        //! and one per message line ( rank only ); LOG_PID identifies
        //! processes, not ranks.
        //! Normal control flow only — syslog() is not async-signal-safe,
        //! so this must never be called from a signal handler
        void
        log_to_syslog(
                const std::string & aLocation,
                const std::string & aCheck,
                const std::vector< std::string > & aMessage );

//------------------------------------------------------------------------------

        template< typename Exception >
        void
        error(
                const std::string    & aLocation,
                const std::string    & aTask,
                const std::string    & aCheck,
                const Exception & aException = Exception()
        )
        {

            std::istringstream tExceptionMessage( aException.what() );
            std::string  tExceptionLine ;


            std::vector< std::string  > tMessage ;
            while ( std::getline( tExceptionMessage, tExceptionLine) )
            {
                tMessage.push_back( tExceptionLine );
            }

            print_errorbox( aLocation, aTask, aCheck, tMessage );

            // the syslog hook sits before the reaction branch, and it must
            // stay there: a debug build still throws ( that is the daily gdb /
            // make check path, serial or parallel ), so an abort-only
            // placement would go silent exactly where a developer reads the
            // message
            if ( syslog_on_error() )
            {
                log_to_syslog( aLocation, aCheck, tMessage );
            }

            if ( throw_on_error() )
            {
                throw aException;
            }

            error_abort();
        }

//------------------------------------------------------------------------------

        std::string
        extract_function_name( const std::string & aPrettyFunction );

//------------------------------------------------------------------------------

        template < typename ... Args >
        void
        belfem_assert(
                const std::string & aFile,
                const std::size_t & aLine,
                const std::string & aFunction,
                const std::string & aCheck,
                const Args ...      aArgs )
        {
            std::stringstream tLocation;

            // get the basename (deliberately not including stringtools here)
            // find last entry of directory delimeter
            tLocation <<  aFile.substr( aFile.find_last_of("/\\") + 1 ) << " (line " << aLine << ")";

            std::stringstream tTask;
            tTask << "complete call to function " << extract_function_name( aFunction ) << "()";

            std::stringstream tReason;
            tReason << "Assertion " << aCheck << " failed.";

            // format output message
            std::string tMessage = belfem::sprint( aArgs ... );

            belfem::assert::error(
                    tLocation.str(),
                    tTask.str(),
                    tReason.str(),
                    std::runtime_error( tMessage.c_str() ) );
        }


//------------------------------------------------------------------------------
    } /* namespace assert */
} /* namespace belfem */

//------------------------------------------------------------------------------

#if BELFEM_ASSERTIONS_ACTIVE
#define BELFEM_ASSERT( aCheck, ... ) \
    do \
    { \
        if ( ! ( aCheck ) ) \
        { \
            belfem::assert::belfem_assert(  \
                __FILE__, \
                __LINE__, \
                __PRETTY_FUNCTION__, \
                #aCheck, \
                __VA_ARGS__ \
                ); \
        } \
    } while ( false )
#else
#define BELFEM_ASSERT( aCheck, ... )
#endif

//---------------------------------------------------------------------------

#define BELFEM_ERROR( aCheck, ... ) \
    do \
    { \
        if ( ! ( aCheck ) ) \
        { \
            belfem::assert::belfem_assert(  \
                __FILE__, \
                __LINE__, \
                __PRETTY_FUNCTION__, \
                #aCheck, \
                __VA_ARGS__ \
                ); \
        } \
    } while ( false )

//------------------------------------------------------------------------------

#endif //BELFEM_ASSERT_HPP
