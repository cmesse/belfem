//
// Created by christian on 7/22/21.
//

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
            mLogFile = aLogFile ;
            mCallgrindFile = tBasePath + ".callgrind" ;
        }

    }

//------------------------------------------------------------------------------

    void
    Profiler::start()
    {
#ifdef BELFEM_PROFILER
        message( 1, "starting profiler ...") ;

        ProfilerStart( mLogFile.c_str() );
#endif
    }

//------------------------------------------------------------------------------

    void
    Profiler::stop()
    {
#ifdef BELFEM_PROFILER

        ProfilerStop();

        message( 1, " ... stopped profiler.") ;

        // get path to executable
        const string & tExecPath = gComm.exec_path() ;

        // assemble command line
        string tCmd = "pprof --callgrind " + tExecPath
                + " " + mLogFile + " > " + mCallgrindFile ;

        string tMessage = "creating callgrind file " + mCallgrindFile + "...";

        message( 1, tMessage.c_str() );

        system( tCmd.c_str() );

        message( 1, " ... done.") ;
#endif
    }

//------------------------------------------------------------------------------
}
