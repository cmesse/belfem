/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * ASCII Art (c) 1996-2001   Joan G. Stark, used with permission.
 *
 */


#include <iomanip>
#include "assert.hpp"
#include "stringtools.hpp"

#include <syslog.h>

#ifdef BELFEM_MPI
#include "cl_Communicator.hpp"
extern belfem::Communicator gComm;
#endif


namespace belfem
{
    namespace assert
    {

//------------------------------------------------------------------------------

        std::string
        extract_function_name( const std::string & aPrettyFunction )
        {
            // Extract the function name from __PRETTY_FUNCTION__
            // Examples:
            //   "virtual void belfem::solver::PETSC::solve(belfem::SpMatrix&, ...)"
            //   -> "belfem::solver::PETSC::solve"
            //
            //   "void belfem::fem::DofManager::compute_jacobian()"
            //   -> "belfem::fem::DofManager::compute_jacobian"

            // Find the opening parenthesis
            size_t tParenPos = aPrettyFunction.find( '(' );
            if ( tParenPos == std::string::npos )
            {
                return aPrettyFunction;  // Fallback if no parenthesis found
            }

            // Work backwards from the parenthesis to find the start of the function name
            // Skip spaces before the parenthesis
            size_t tEnd = tParenPos;
            while ( tEnd > 0 && std::isspace( aPrettyFunction[ tEnd - 1 ] ) )
            {
                --tEnd;
            }

            // Find the start: look for the last space before the function name
            // This handles "virtual void Class::function" -> starts after "void "
            size_t tStart = aPrettyFunction.rfind( ' ', tEnd - 1 );
            if ( tStart == std::string::npos )
            {
                tStart = 0;
            }
            else
            {
                ++tStart;  // Skip the space
            }

            return aPrettyFunction.substr( tStart, tEnd - tStart );
        }

//------------------------------------------------------------------------------

        std::vector<std::string> wrap_lines(std::size_t aMaxWidth, const std::string& aLine)
        {
            std::vector<std::string> aWrappedLines;
            std::istringstream tStream(aLine);
            std::string tWord;
            std::string tLine;

            while (tStream >> tWord)
            {
                // Check if the word fits in the current line
                if (utf8_character_count(tLine) + utf8_character_count(tWord) + (tLine.empty() ? 0 : 1) <= aMaxWidth)
                {
                    if (!tLine.empty()) tLine += " ";
                    tLine += tWord;
                }
                else
                {
                    // If the word itself is longer than aMaxWidth, split it
                    if (utf8_character_count(tWord) > aMaxWidth)
                    {
                        // If there's any content in the current line, add it to aWrappedLines
                        if (!tLine.empty())
                        {
                            aWrappedLines.push_back(tLine);
                            tLine.clear();
                        }

                        // Split the long word into chunks of aMaxWidth characters
                        size_t tStart = 0;
                        size_t tBytePos = 0;
                        size_t tCharCount = 0;

                        while ( tBytePos < tWord.length() )
                        {
                            size_t tNextBytePos = tBytePos;
                            unsigned char c = tWord[tNextBytePos++];
                            if (c <= 0x7F) {}           // 1-byte (ASCII)
                            else if (c <= 0xDF) tNextBytePos++;
                            else if (c <= 0xEF) tNextBytePos += 2;
                            else tNextBytePos += 3;

                            if ( ++tCharCount > aMaxWidth )
                            {
                                aWrappedLines.push_back(tWord.substr(tStart, tBytePos - tStart));
                                tStart = tBytePos;
                                tCharCount = 1;
                            }
                            tBytePos = tNextBytePos;
                        }

                        if (tStart < tWord.length()) aWrappedLines.push_back(tWord.substr(tStart));
                    }
                    else
                    {
                        // Start a new line with the word
                        if (!tLine.empty()) aWrappedLines.push_back(tLine);
                        tLine = tWord;
                    }
                }
            }

            // Add the last line if it's not empty
            if (!tLine.empty()) aWrappedLines.push_back( tLine );

            return aWrappedLines;
        }

//------------------------------------------------------------------------------

        void
        hatch_dragon( std::vector< std::string > & aDragon )
        {
            aDragon.clear();
            aDragon.push_back("    ║                                " );
            aDragon.push_back("    ║                     ,          " );
            aDragon.push_back("    ║                    (\\\\__       " );
            aDragon.push_back("    ║                   C/`   `\\     " );
            aDragon.push_back("    ║                  C|    <oo'--  " );
            aDragon.push_back("    ║                   C\\    ,___\"/ ");
            aDragon.push_back("    ║                    C) \\___/ `  ");
            aDragon.push_back("    ║                   C/    |-.    ");
            aDragon.push_back("    ║                  C/    '-. \\   ");
            aDragon.push_back("    ║   .-=-.         C|   \\__  \\,)  ");
            aDragon.push_back("    ║     `) )       C/       \\,,)   ");
            aDragon.push_back("    ║    .' /       C' _    .==|     ");
            aDragon.push_back("    ║   (  (   _.-'`    \\  .==/      ");
            aDragon.push_back("    ║    \\  `'`     \\   | .==/\\__    ");
            aDragon.push_back("    ║ jgs `-.,,-'``-/  /__=;`   /    ");
            aDragon.push_back("    ║              (    `\\  \\_.'     ");
            aDragon.push_back("    ║               `\"\"\"\"\"`          ");
            aDragon.push_back("    ║                                ");
        }

//------------------------------------------------------------------------------

        void
        print_line(const std::vector<std::string> & aDragon, std::size_t & aCounter )
        {
            if ( aCounter < aDragon.size() )
            {
                std::cerr << aDragon[ aCounter++ ] <<  std::setw( 43 ) << "║" << std::endl ;
            }
            else
            {
                ++aCounter ;
            }
        }

//------------------------------------------------------------------------------

        void
        print_line( const std::vector<std::string> & aDragon, const string & aLine, std::size_t & aCounter )
        {
            std::vector< std::string > tWrappedLines = wrap_lines( 39, aLine );

            for ( const std::string & tLine : tWrappedLines )
            {
                if( aCounter < aDragon.size() )
                {
                    std::cerr << aDragon[ aCounter ] << "  " << tLine <<  std::setw( 38 - utf8_character_count( tLine ) ) << "" << "║" << std::endl ;
                }
                else
                {
                    std::cerr << aDragon[ 0 ] << "  " << tLine << std::setw( 38 - utf8_character_count( tLine ) ) << "" << "║" << std::endl ;
                }
                ++aCounter ;
            }
        }
//------------------------------------------------------------------------------

        void
        print_line( const string & aLine )
        {
            std::vector< std::string > tWrappedLines = wrap_lines( 68, aLine );

            for ( const std::string & tLine : tWrappedLines )
            {
                std::cerr << "    ║   " << tLine <<  std::setw( 69 - utf8_character_count( tLine ) ) << "" << "║" << std::endl ;
            }
        }

//------------------------------------------------------------------------------

        void
        get_lines( const string & aWhat, std::vector< std::string > & aLines )
        {
            aLines.clear();
            std::istringstream tStream( aWhat );

            string tLine ;

            while ( std::getline( tStream, tLine) )
            {
                aLines.push_back( tLine );
            }
        }

//------------------------------------------------------------------------------
        void
        print_errorbox(
            const std::string    & aLocation,
            const std::string    & aTask,
            const std::string    & aCheck,
            const std::vector< std::string  > & aMessage )
        {
            std::vector< std::string > tDragon ;
            hatch_dragon( tDragon );

            std::string tErrorLine   = aTask + "." ;
            std::string tReasonLine  = "Reason :    " + aCheck + "." ;
            std::string tWhereLine   = "Where  :     " + aLocation + "." ;
#ifdef BELFEM_MPI
            std::string tProcline    = "This error occured on proc " + std::to_string( gComm.rank() );
#endif
            std::cerr << std::endl ;
            std::cerr << "    ╔════════════════════════════════════════════════════════════════════════╗" << std::endl ;
            std::size_t tCount = 0 ;
            for( uint k=0; k<3; ++k )
            {
                print_line( tDragon,tCount );
            }
            print_line( tDragon, "ERROR:", tCount );
            print_line( tDragon,tCount );
            print_line( tDragon, "Unable to " + tErrorLine, tCount );
            print_line( tDragon, tCount );
            print_line( tDragon, "REASON:", tCount );
            print_line( tDragon, tCount );
            print_line( tDragon, aCheck, tCount );
            print_line( tDragon, tCount );
            print_line( tDragon, "WHERE:", tCount );
            print_line( tDragon, tCount );
            print_line( tDragon, aLocation, tCount );
            print_line( tDragon, tCount );
            while( tCount < tDragon.size() )
            {
                print_line( tDragon, tCount );
            }
            std::cerr << "    ║                                                                        ║" << std::endl ;
            std::cerr << "    ╠════════════════════════════════════════════════════════════════════════╣" <<  std::endl ;
            std::cerr << "    ║                                                                        ║" << std::endl ;
            for ( const std::string  & tLine : aMessage )
            {
                 print_line( tLine );
            }
#ifdef BELFEM_MPI
            if( gComm.size() != gNoOwner && gComm.size() > 1 )
            {
                std::cerr << "    ║                                                                        ║" << std::endl ;
                print_line( tProcline  );
            }
#endif
            std::cerr << "    ║                                                                        ║" << std::endl ;
            std::cerr << "    ╚════════════════════════════════════════════════════════════════════════╝" << std::endl ;
            std::cerr << std::endl ;
        }

        void
        error_abort()
        {
            // the MPI lifecycle contract that used to live here - the
            // MPI_COMM_WORLD-not-gComm.world() choice, the Initialized /
            // Finalized guards, the unconditional std::abort() fall-through -
            // moved verbatim into comm_abort ( declared in cl_Communicator.hpp,
            // defined in commtools.cpp ) on 2026-08-30, so that every caller
            // inherits it rather than only this one. L-21: vendor calls go
            // through the wrapper.
            comm_abort( 1 );
        }

//------------------------------------------------------------------------------

        // Reaction of a failed check; see throw_on_error() in the header for why
        // this lives here and not as a static in assert.hpp. Initialized to the
        // build's compile-time behaviour so nothing changes unless a test asks.
        static bool gThrowOnError = BELFEM_ASSERTIONS_ACTIVE ;

//------------------------------------------------------------------------------

        bool
        throw_on_error()
        {
            return gThrowOnError ;
        }

//------------------------------------------------------------------------------

        void
        set_throw_on_error( const bool aValue )
        {
            gThrowOnError = aValue ;
        }

//------------------------------------------------------------------------------

        // syslog reaction; independent of gThrowOnError — see the header
        static bool gSyslogOnError = true ;

//------------------------------------------------------------------------------

        bool
        syslog_on_error()
        {
            return gSyslogOnError ;
        }

//------------------------------------------------------------------------------

        void
        set_syslog_on_error( const bool aValue )
        {
            gSyslogOnError = aValue ;
        }

//------------------------------------------------------------------------------

        void
        log_to_syslog(
                const std::string & aLocation,
                const std::string & aCheck,
                const std::vector< std::string > & aMessage )
        {
            // per-event open/close: this path is terminal ( abort or a
            // throw that ends the run ), so persistent state buys nothing
            openlog( "belfem", LOG_PID, LOG_USER );

            // rank in the payload: LOG_PID identifies local processes,
            // not MPI ranks, and rank-local checks fire on ONE rank —
            // which may not be rank 0, so no rank gate here either.
            // Before gComm.init() the rank is the gNoOwner sentinel; an
            // error that early logs "rank ?" instead of 2147483647
#ifdef BELFEM_MPI
            const std::string tRank = gComm.rank() == gNoOwner ?
                std::string( "?" ) : std::to_string( gComm.rank() ) ;
#else
            const std::string tRank = "0" ;
#endif
            syslog( LOG_CRIT, "rank %s: %s: %s",
                tRank.c_str(),
                aLocation.c_str(),
                aCheck.c_str() );

            // one line per message line: syslog truncates around 1 KiB
            // per call, and the multi-line error boxes survive per-line
            for ( const std::string & tLine : aMessage )
            {
                if ( tLine.size() > 0 )
                {
                    syslog( LOG_CRIT, "rank %s: %s",
                        tRank.c_str(), tLine.c_str() );
                }
            }

            closelog();
        }

//------------------------------------------------------------------------------
    }
}
