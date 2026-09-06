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
 */

#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#ifdef  OMP
#include <omp.h>
#endif
#include "banner.hpp"
#include "belfem_version.hpp"
#include "cl_Communicator.hpp"
#include "assert.hpp"
#include "filetools.hpp"

// Externally Defined Global Communicator
extern belfem::Communicator gComm;

namespace belfem
{
//------------------------------------------------------------------------------

    std::string
    exec( const std::string & aCommand )
    {
        // output
        std::string aString;

        // size for buffer
        const int tBufferSize = 128;

        // temporary buffer
        char tBuffer[ tBufferSize ];

        // create pointer for stream object
        std::shared_ptr<FILE> tStream( popen( aCommand.c_str(), "r" ), pclose );

        if ( tStream )
        {
            // read result from command line
            while ( !feof( tStream.get() ) )
            {
                if ( fgets( tBuffer, tBufferSize, tStream.get() ) != nullptr )
                {
                    aString.append( tBuffer );
                }
            }

            // trim string
            auto tStart = aString.find_first_of(':') + 1;
            auto tEnd   = aString.find_last_not_of('\n');

            // return info
            return aString.substr( tStart, ( tEnd - tStart + 1 ) );

        }
        else
        {
            BELFEM_ASSERT( false, "could not execute command %s", aCommand.c_str() );

            return aString;
        }
    }
//------------------------------------------------------------------------------

    std::string
    uname()
    {
        return exec("uname");
    }

//------------------------------------------------------------------------------

    std::string
    cpu_info()
    {
        // get os type
        std::string tUname( uname() );

        if( tUname == "Darwin" )
        {
            return exec( "sysctl -n machdep.cpu.brand_string" );
        }
        else if ( tUname == "Linux" )
        {
            std::ifstream tProcCpuInfo("/proc/cpuinfo");

            // test if file exists
            if (tProcCpuInfo)
            {
                return exec( "cat /proc/cpuinfo | grep \"model name\" | head -n 1 2>&1" );
            } else
            {
                return "unknown";
            }
        }
        else
        {
            return "unknown";
        }
    }

//------------------------------------------------------------------------------

    std::string
    os_string()
    {
        // get os type
        std::string tUname( uname() );

        if( tUname == "Darwin" )
        {
            std::string tName = clean_string( exec( "sw_vers | grep ProductName | cut -d: -f2") );
            std::string tVersion = clean_string( exec( "sw_vers | grep ProductVersion | cut -d: -f2") );
            std::string tBuild = clean_string( exec( "sw_vers | grep BuildVersion | cut -d: -f2") );
            return tName + " " + tVersion + " " + tBuild ;
        }
        else if ( tUname == "Linux" )
        {
            if ( file_exists( "/etc/system-release" ))
            {
                std::string tLabel = exec( "cat /etc/system-release | cut -d'(' -f 1 " ) ;
                return clean_string( search_and_replace( tLabel, "release", "" ) );
            }
            else if ( exec( "which lsb_release" ).size() > 0 )
            {
                return clean_string( exec( "lsb_release -d | grep Description | cut -d: -f2" ) );
            }
            else
            {
                return tUname ;
            }
        }
        else
        {
            return tUname ;
        }
    }

//------------------------------------------------------------------------------

    std::string
    version()
    {
        return gVersionString;
    }

//------------------------------------------------------------------------------

    bool
    is_built_from_git()
    {
        return gGitAvailable;
    }

//------------------------------------------------------------------------------

    std::string
    git_commit_hash()
    {
        return gGitCommitHash;
    }

//------------------------------------------------------------------------------

    std::string
    git_commit_hash_short()
    {
        return gGitCommitHashShort;
    }

//------------------------------------------------------------------------------

    std::string
    git_branch()
    {
        return gGitBranch;
    }

//------------------------------------------------------------------------------

    bool
    git_is_dirty()
    {
        return gGitIsDirty;
    }

//-----------------------------------------------------------------------------

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#elif BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#endif

    void
    print_banner( const std::string aExecName )
    {
        // banner is only printed by first proc
        if( gComm.rank() == 0 )
        {
            // Get current time
            std::time_t tSystemTime = std::time(NULL);
            // Convert to local time structure
            std::tm* tCalendarTime = std::localtime(&tSystemTime);

            int tDay   = tCalendarTime->tm_mday ;
            int tMonth = tCalendarTime->tm_mon + 1 ;
            int tWday  = tCalendarTime->tm_wday ;
            int tYear  = tCalendarTime->tm_year + 1900 ;

            bool tJGS = false ;

            // Get system locale name (e.g., "en_US.UTF-8")
            std::string locale_name;
            try
            {
                std::locale user_locale("");
                locale_name = user_locale.name();
            }
            catch (...)
            {
                locale_name = "en_US.UTF-8";
            }

            // Extract country code (after '_', take 2 letters before '.' if present)
            std::string country_code = "US";  // Default to US
            size_t underscore_pos = locale_name.find('_');
            if (underscore_pos != std::string::npos && underscore_pos + 2 < locale_name.size()) {
                size_t dot_pos = locale_name.find('.', underscore_pos);
                country_code = locale_name.substr(underscore_pos + 1, (dot_pos != std::string::npos) ? (dot_pos - underscore_pos - 1) : 2);
                // Uppercase for consistency (e.g., "us" -> "US")
                std::transform(country_code.begin(), country_code.end(), country_code.begin(), ::toupper);
            }

            // Check for holidays and display themed banners
            switch ( tMonth )
            {
                case 3 :
                case 4:
                {
                    // St. Patrick's Day - March 17
                    if ( tDay == 17 && tMonth == 3 )
                    {
                        tJGS = banners::print_stpatrick() ;
                    }
                    else
                    {
                        // Easter calculation using the Anonymous Gregorian algorithm (Computus)
                        // This algorithm calculates the date of Easter Sunday for any Gregorian year
                        // Easter falls on the first Sunday after the first ecclesiastical full moon
                        // that occurs on or after March 21 (the ecclesiastical "vernal equinox")

                        int a = tYear % 19;            // Golden number - 1 (position in Metonic cycle)
                        int b = tYear / 100;           // Century
                        int c = tYear % 100;           // Year within century
                        int d = b / 4;                 // Number of leap centuries
                        int e = b % 4;                 // Remainder of leap centuries
                        int f = (b + 8) / 25;          // Correction for lunar orbit
                        int g = (b - f + 1) / 3;       // Correction for Gregorian calendar
                        int h = (19 * a + b - d - g + 15) % 30;  // Days from March 21 to Paschal full moon
                        int i = c / 4;                 // Number of leap years within century
                        int k = c % 4;                 // Remainder of leap years
                        int l = (32 + 2 * e + 2 * i - h - k) % 7;  // Number of days from Paschal full moon to next Sunday
                        int m = (a + 11 * h + 22 * l) / 451;  // Correction for dates that fall after April 25

                        // Calculate Easter Sunday date
                        std::tm tEtm = {};
                        tEtm.tm_mday = ((h + l - 7 * m + 114) % 31) + 1 ;  // Day of month (1-31)
                        tEtm.tm_mon = (h + l - 7 * m + 114) / 31 - 1 ;     // Month (0-based: 2=March, 3=April)
                        tEtm.tm_year =tCalendarTime->tm_year ;

                        // Calculate days from current date to Easter Sunday
                        int tDiff = static_cast<int>(std::difftime( std::mktime(&tEtm), std::mktime(tCalendarTime)) / (60 * 60 * 24));

                        // Display Easter bunny banner during Easter week (0-6 days before Easter Sunday)
                        if ( tDiff >= 0 && tDiff <= 6 )
                        {
                            tJGS = banners::print_easter();
                        }
                    }
                    if ( ! tJGS )
                    {
                        tJGS = banners::print_default();
                    }
                    break ;
                }
                case 7 :
                {
                    // Canada Day - July 1 (Canadian national day)
                    if ( tDay == 1 && country_code == "CA" )
                    {
                        tJGS = banners::print_canada();
                    }


                    // US Independence Day - July 4
                    else if ( tDay == 4 && country_code == "US" )
                    {
                        tJGS = banners::print_usa();
                        banners::print_default();
                    }
                    if ( ! tJGS )
                    {
                        tJGS = banners::print_default();
                    }
                    break ;
                }
                case 10 :
                case 11 :
                {
                    // non-sensical values to sppress special banner by default
                    int week_start = 100 ;
                    int week_end = 0 ;

                    // Canadian Thanksgiving - second Monday in October
                    if ( country_code == "CA" && tMonth == 10 )
                    {
                        // Algorithm to find the second Monday of October:
                        // 1. Back-calculate what day of week October 1 was
                        // 2. Find the first Monday
                        // 3. Add one week to get the second Monday
                        // 4. Display banner for the entire week (Sunday-Saturday)

                        // Compute weekday of Oct 1 by working backwards from current day
                        // tm_wday: 0=Sunday, 1=Monday, ..., 6=Saturday
                        int weekday_oct1 = (tWday - ((tDay - 1) % 7) + 7) % 7;

                        // Calculate how many days until the first Monday (day 1)
                        int days_to_first_mon = (1 - weekday_oct1 + 7) % 7;
                        int first_mon = 1 + days_to_first_mon;

                        // Second Monday is one week after first Monday (Canadian Thanksgiving)
                        int td = first_mon + 7;

                        // Calculate the full week range (Sunday before to Saturday after)
                        int weekday_td = (weekday_oct1 + (td - 1)) % 7;  // Weekday of second Monday (should be 1)
                        week_start = td - weekday_td;  // Sunday of Thanksgiving week
                        week_end = week_start + 6;     // Saturday of Thanksgiving week
                    }

                    // US Thanksgiving - fourth Thursday in November
                    else if ( country_code == "US" && tMonth == 11 )
                    {
                        // Algorithm to find the fourth Thursday of November:
                        // 1. Back-calculate what day of week November 1 was
                        // 2. Find the first Thursday (day 4 in tm_wday: 0=Sun, 1=Mon, ..., 4=Thu, ..., 6=Sat)
                        // 3. Add three weeks (21 days) to get the fourth Thursday
                        // 4. Display banner for the entire week (Sunday-Saturday)

                        // Compute weekday of Nov 1 by working backwards from current day
                        int weekday_nov1 = (tWday - ((tDay - 1) % 7) + 7) % 7;

                        // Calculate how many days until the first Thursday (day 4)
                        int days_to_first_thu = (4 - weekday_nov1 + 7) % 7;
                        int first_thu = 1 + days_to_first_thu;

                        // Fourth Thursday is three weeks after first Thursday (US Thanksgiving)
                        int td = first_thu + 21;  // +3 weeks = +21 days

                        // Calculate the full week range (Sunday before to Saturday after)
                        // Since Thursday is weekday 4, go back 4 days to Sunday and forward 2 days to Saturday
                        week_start = td - 4;  // Sunday of Thanksgiving week
                        week_end = td + 2;    // Saturday of Thanksgiving week
                    }

                    // Display turkey banner during Thanksgiving week
                    if ( tDay >= week_start && tDay <= week_end )
                    {
                        tJGS = banners::print_thanksgiving();
                    }
                    if ( ! tJGS )
                    {
                        tJGS = banners::print_default();
                    }
                    break;
                }
                case 12 :
                {
                    // December / Holiday season - display festive tree for entire month
                    tJGS = banners::print_christmas();
                    break ;
                }

                default:
                {
                    tJGS = banners::print_default();
                }
            }
            if( aExecName.length() > 0 )
            {
                std::fprintf( stdout, "\n    %s\n\n",aExecName.c_str() );
            }

            /*if ( is_built_from_git() )
            {
                std::fprintf( stdout, "    Version %s GIT '%s' %s ", version().c_str(), git_branch().c_str(), git_commit_hash_short().c_str() );
            }
            else
            {
                std::fprintf( stdout, "    Version %s", version().c_str() );
            }*/

            std::fprintf( stdout, "\n\n" );
            std::fprintf( stdout, "    Copyright (c) 2026 The Regents of the University of California,\n");
            std::fprintf( stdout, "    through Lawrence Berkeley National Laboratory (subject to receipt,\n");
            std::fprintf( stdout, "    of any required approvals from the U.S. Dept. of Energy).\n");
            std::fprintf( stdout, "    All rights reserved.\n\n");
            std::fprintf( stdout, "    Developers: Christian Messe, Gregory Giard\n\n");
            std::fprintf( stdout, "    \n\n");
            std::fprintf( stdout, "    See the top-level LICENSE file for the complete license and disclaimer.\n\n");
            std::fprintf( stdout, "    \n\n");
            if ( tJGS )
            {
                std::fprintf( stdout, "    ASCII Art (c) 1996-2001   Joan G. Stark, used with permission.\n\n");
                std::fprintf( stdout, "    \n\n");
            }
            // Parallel flags
#ifdef  BELFEM_MPI
            std::fprintf( stdout, "    using MPI\n\n" );
#endif
#ifdef  OMP
            std::fprintf( stdout, "    using OpenMP\n\n" );
#endif
#ifdef BELFEM_PROFILER
            std::fprintf( stdout, "    using Google Profiler\n\n" );
#endif
#if !defined(NDEBUG) || defined(DEBUG)
            std::fprintf( stdout, "    DEBUG flags are on.\n\n");
#endif


            // What version?
            std::fprintf( stdout, "    Version       : %s\n", version().c_str() );

            // When built?
            std::fprintf( stdout, "    Build Date    : %s at %s\n", __DATE__ , __TIME__ );

            if ( is_built_from_git() )
            {
                // What commit?
                if ( git_is_dirty() )
                {
                    std::fprintf( stdout, "    Git Commit    : %s (%s, with uncommitted changes)\n",
                                  git_commit_hash_short().c_str(), git_branch().c_str() );
                }
                else
                {
                    std::fprintf( stdout, "    Git Commit    : %s (%s)\n",
                                  git_commit_hash_short().c_str(), git_branch().c_str() );
                }
            }

            // insert blank line
            std::fprintf( stdout, "\n");


            // Who?
            std::fprintf( stdout, "    User/Host     : %s @ %s \n",
                          exec( "whoami").c_str(),
                          exec( "hostname").c_str() );

            // operating system
            std::fprintf( stdout, "    System        : %s ( %s ) \n",
                                  os_string().c_str(),
                          exec( "uname -m").c_str() );


            std::string tCpuInfo = cpu_info();

            std::fprintf( stdout, "    CPU Info      : %s\n", clean_string( tCpuInfo.c_str() ).c_str() );

#ifdef  BELFEM_MPI
            std::fprintf( stdout, "    Procs Used    : %i\n", ( int ) gComm.size() );
#endif
// NOT a mistake: OMP, not BELFEM_OMP. This number is what the THIRD-PARTY solvers
// will thread with; it stays meaningful on a default build, where BELFEM_OMP is OFF
// and BELFEM's own kernels are serial.
#ifdef  OMP
            std::fprintf( stdout, "    Threads Used  : %i\n", omp_get_max_threads() );
#endif

            // insert blank line
            std::fprintf( stdout, "\n");

            // What Compiler
#ifdef BELFEM_GCC
            std::fprintf( stdout, "    Compiler      : GNU Compiler Collection\n");
#elif BELFEM_CLANG
            std::fprintf( stdout, "    Compiler      : Apple Clang\n");
#elif BELFEM_INTEL
            std::fprintf( stdout, "    Compiler      : Intel oneAPI\n");
#elif BELFEM_PGI
            std::fprintf( stdout, "    Compiler      : NVidia PGI \n");
#endif

            // What LAPACK Lib
#ifdef BELFEM_NETLIB
            std::fprintf( stdout, "    BLAS & LAPACK : Netlib\n");
#elif BELFEM_ACCELLERATE
            std::fprintf( stdout, "    BLAS & LAPACK : Apple Accellerate Framework\n");
#elif BELFEM_MKL
            std::fprintf( stdout, "    BLAS & LAPACK : Intel Math Kernel Library (MKL)\n");
#endif


            // What Matrix lib?
#ifdef BELFEM_ARMADILLO
            std::fprintf( stdout, "    Matrix Lib    : Armadillo\n");
#elif  BELFEM_BLAZE
            std::fprintf( stdout, "    Matrix Lib    : Blaze\n");
#endif
            // What Solver libs

            string tSolverLibs = "";
#ifdef BELFEM_SUITESPARSE
            tSolverLibs += "Suite Sparse, ";
#endif
#ifdef BELFEM_SUPERLU
            tSolverLibs += "SuperLU, ";
#endif
#ifdef BELFEM_MUMPS
            tSolverLibs += "MUMPS, ";
#endif
#ifdef BELFEM_STRUMPACK
            tSolverLibs += "STRUMPACK, ";
#endif
#ifdef BELFEM_PARDISO
            tSolverLibs += "PARDISO, ";
#endif
#ifdef BELFEM_PETSC
            tSolverLibs += "PETSc, ";
#endif

            tSolverLibs = clean_string( tSolverLibs ).substr( 0, tSolverLibs.length() - 2 );
            std::fprintf( stdout, "    Solver Libs   : %s\n", tSolverLibs.c_str() );

            // insert blank line
            std::fprintf( stdout, "\n");

            // What?
            std::fprintf( stdout, "    Executable    : %s\n", gComm.exec_path().c_str() );
            std::fprintf( stdout, "    Arguments     : %s\n", gComm.argument_string().c_str()  );

            // Where?
            std::fprintf( stdout, "    Run Dir       : %s\n", gComm.workdir().c_str() );
            std::fprintf( stdout, "\n");

        }
    }

    bool
    banners::print_default()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "    .______    _______  __       _______  _______ .___  ___.\n");
        std::fprintf( stdout, "    |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
        std::fprintf( stdout, "    |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
        std::fprintf( stdout, "    |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
        std::fprintf( stdout, "    |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
        std::fprintf( stdout, "    |______/  |_______||_______||__|     |_______||__|  |__|\n\n");

        std::fprintf( stdout, "    %s\n\n", gLongName.c_str());
        std::fprintf( stdout, "    %s\n\n", gURL.c_str());

        return false ;
    }

    bool
    banners::print_easter()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "         _     _\n");
        std::fprintf( stdout, "        /\\`\\ /`/\\      .______    _______  __       _______  _______ .___  ___.\n");
        std::fprintf( stdout, "        \\/\\ V /\\/      |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
        std::fprintf( stdout, "          /6 6\\        |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
        std::fprintf( stdout, "         (= Y =)       |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
        std::fprintf( stdout, "         /`\"^\"`\\       |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
        std::fprintf( stdout, "        / /   \\ \\      |______/  |_______||_______||__|     |_______||__|  |__|\n");
        std::fprintf( stdout, "       (_/     \\_)\n");
        std::fprintf( stdout, "        /       \\o     %s\n", gLongName.c_str());
        std::fprintf( stdout, " jgs ___\\       /___\n");
        std::fprintf( stdout, "    (((____/^\\____)))  %s\n\n", gURL.c_str());
        return true ;
    }

    bool
    banners::print_stpatrick()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "       .-.-.     .______    _______  __       _______  _______ .___  ___.\n");
        std::fprintf( stdout, "      (     )    |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
        std::fprintf( stdout, "    .-.\\ : /.-.  |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
        std::fprintf( stdout, "   (   .`:`.   ) |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
        std::fprintf( stdout, "    (   /|\\   )  |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
        std::fprintf( stdout, " jgs `\"` | `\"`   |______/  |_______||_______||__|     |_______||__|  |__|\n\n");
        std::fprintf( stdout, "                 %s\n\n", gLongName.c_str());
        std::fprintf( stdout, "                 %s\n\n", gURL.c_str());

        return true ;
    }

    bool
    banners::print_usa()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "                                            o\n");
        std::fprintf( stdout, "                                           /\\\n");
        std::fprintf( stdout, "                                          /::\\\n");
        std::fprintf( stdout, "                                         /::::\\\n");
        std::fprintf( stdout, "                           ,a_a         /\\::::/\\\n");
        std::fprintf( stdout, "                          {/ ''\\_      /\\ \\::/\\ \\\n");
        std::fprintf( stdout, "                          {\\ ,_oo)    /\\ \\ \\/\\ \\ \\\n");
        std::fprintf( stdout, "                          {/  (_^____/  \\ \\ \\ \\ \\ \\\n");
        std::fprintf( stdout, "                .=.      {/ \\___)))*)    \\ \\ \\ \\ \\/\n");
        std::fprintf( stdout, "               (.=.`\\   {/   /=;  ~/      \\ \\ \\ \\/\n");
        std::fprintf( stdout, "                   \\ `\\{/(   \\/\\  /        \\ \\ \\/\n");
        std::fprintf( stdout, "                    \\  `. `\\  ) )           \\ \\/\n");
        std::fprintf( stdout, "                jgs  \\    // /_/_            \\/\n");
        std::fprintf( stdout, "                       '==''---)))");
        return true ;
    }

    bool
    banners::print_canada()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "                  .______    _______  __       _______  _______ .___  ___.\n");
        std::fprintf( stdout, "    . |`|/| .     |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
        std::fprintf( stdout, "    |\\|\\|'|/|     |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
        std::fprintf( stdout, " .--'-\\`|/-''--.  |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
        std::fprintf( stdout, "  \\`-._\\|./.-'/   |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
        std::fprintf( stdout, "   >`-._|/.-'<    |______/  |_______||_______||__|     |_______||__|  |__|\n");
        std::fprintf( stdout, "  '~|/~~|~~\\|~\n");

        std::fprintf( stdout, " jgs    |         %s\n\n", gLongName.c_str());
        std::fprintf( stdout, "                  %s\n\n", gURL.c_str());
        return true ;
    }

    bool
    banners::print_thanksgiving()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "\n");
        std::fprintf( stdout, "                   .______    _______  __       _______  _______ .___  ___.\n");
        std::fprintf( stdout, "            .-.    |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
        std::fprintf( stdout, "    .;;;;. ( ^_>   |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
        std::fprintf( stdout, "   <;<;  \\;>\\ !    |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
        std::fprintf( stdout, "  <;<;   '-.>) \\   |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
        std::fprintf( stdout, "   <;<; <'=.    |  |______/  |_______||_______||__|     |_______||__|  |__|\n");
        std::fprintf( stdout, "   <;<; '-     /\n");
        std::fprintf( stdout, "     <;,\\.\\--'`    %s\n", gLongName.c_str());
        std::fprintf( stdout, " jgs    `==`==\n");
        std::fprintf( stdout, "                   %s\n\n", gURL.c_str());

        return true ;
    }

    bool
    banners::print_christmas()
    {
        std::fprintf( stdout, "\n\n");
        std::fprintf( stdout, "               .______    _______  __       _______  _______ .___  ___.\n");
        std::fprintf( stdout, "      \\'/      |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
        std::fprintf( stdout, "    -= * =-    |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
        std::fprintf( stdout, "      {.}      |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
        std::fprintf( stdout, "     {.-'}     |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
        std::fprintf( stdout, "    {`_.-'}    |______/  |_______||_______||__|     |_______||__|  |__|\n");
        std::fprintf( stdout, "   {-` _.-'}\n");
        std::fprintf( stdout, "    `\":=:\"`    %s\n", gLongName.c_str());
        std::fprintf( stdout, " jgs `---`\n");
        std::fprintf( stdout, "               %s\n\n", gURL.c_str());

        return true ;
    }

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_CLANG
#pragma clang diagnostic pop
#endif



//-----------------------------------------------------------------------------
}
