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

#include <iostream>

#ifdef BELFEM_PROFILER
#include <gperftools/profiler.h>
#endif
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_Profiler.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Profiler::Profiler( const string aLogFile )
    {
        // generate the paths
        string tBasePath = aLogFile.substr( 0, aLogFile.find_last_of("."));
#ifdef BELFEM_MPI
        if( gComm.size() > 1 )
        {
            // parallel mode
            string tSuffix = "." + std::to_string( gComm.size() ) + "." + std::to_string( gComm.rank() );

            // parallel path of logfile
            mLogFile = tBasePath + tSuffix + aLogFile.substr( aLogFile.find_last_of("."), aLogFile.length() );
            mCallgrindFile = tBasePath + tSuffix + ".callgrind" ;

        }
        else
        {
#endif
            mLogFile = aLogFile ;
            mCallgrindFile = tBasePath + ".callgrind" ;
#ifdef BELFEM_MPI
        }
#endif
    }

//------------------------------------------------------------------------------

    void
    Profiler::start()
    {
#ifdef BELFEM_PROFILER
        message( InfoLevel::Default, "\n--------------------------------------------------------------------------------") ;
        message( InfoLevel::Default, "    starting profiler ...") ;
        message( InfoLevel::Default, "--------------------------------------------------------------------------------\n") ;
        ProfilerStart( mLogFile.c_str() );
#endif
    }

//------------------------------------------------------------------------------

    void
    Profiler::stop()
    {
#ifdef BELFEM_PROFILER

        ProfilerStop();

        message( InfoLevel::Default, "\n--------------------------------------------------------------------------------") ;
        message( InfoLevel::Default, "    ... stopped profiler") ;
        message( InfoLevel::Default, "--------------------------------------------------------------------------------\n") ;

        // get path to executable
        const string & tExecPath = basename( gComm.exec_path() ) ;

        // assemble command line
        string tCmd = "pprof --callgrind $(which " + tExecPath
                + ") " + mLogFile + " > " + mCallgrindFile ;

        message(  InfoLevel::Default, tCmd.c_str() );

        string tMessage = "creating callgrind file " + mCallgrindFile + "...";

        message(  InfoLevel::Default, tMessage.c_str() );

        system( tCmd.c_str() );

        message(  InfoLevel::Default, " ... done.") ;
#endif
    }

//------------------------------------------------------------------------------
}
