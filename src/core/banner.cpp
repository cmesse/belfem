//
// Created by Christian Messe on 2018-12-30.
//

#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

#include "banner.hpp"
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
            std::fprintf( stdout, "    .______    _______  __       _______  _______ .___  ___.\n");
            std::fprintf( stdout, "    |   _  \\  |   ____||  |     |   ____||   ____||   \\/   |\n");
            std::fprintf( stdout, "    |  |_)  | |  |__   |  |     |  |__   |  |__   |  \\  /  |\n");
            std::fprintf( stdout, "    |   _  <  |   __|  |  |     |   __|  |   __|  |  |\\/|  |\n");
            std::fprintf( stdout, "    |  |_)  | |  |____ |  `----.|  |     |  |____ |  |  |  |\n");
            std::fprintf( stdout, "    |______/  |_______||_______||__|     |_______||__|  |__|\n\n");

            std::fprintf( stdout, "    The Berkeley Lab Finite Element Framework\n\n");

            if( aExecName.length() > 0 )
            {
                std::fprintf( stdout, "\n    %s\n\n",aExecName.c_str() );
            }
            else
            {
                std::fprintf( stdout, "\n\n");
            }

            std::fprintf( stdout, "    Copyright (c) 2023 The Regents of the University of California,\n");
            std::fprintf( stdout, "    through Lawrence Berkeley National Laboratory (subject to receipt,\n");
            std::fprintf( stdout, "    of any required approvals from the U.S. Dept. of Energy).\n");
            std::fprintf( stdout, "    All rights reserved.\n\n");
            std::fprintf( stdout, "    Developed by Christian Messe  - cmesse@lbl.gov\n\n");

// Parallel flags
#ifdef  BELFEM_MPI
            std::fprintf( stdout, "    using MPI\n\n" );
#endif
#ifdef  BELFEM_OPENMP
            std::fprintf( stdout, "    using OpenMP\n\n" );
#endif

#if !defined(NDEBUG) || defined(DEBUG)
            std::fprintf( stdout, "    DEBUG flags are on.\n\n");
#endif

            // Who?
            std::fprintf( stdout, "    User/Host     : %s @ %s \n",
                          exec( "whoami").c_str(),
                          exec( "hostname").c_str() );

            // operating system
            std::fprintf( stdout, "    System        : %s ( %s ) \n",
                                  os_string().c_str(),
                          exec( "uname -m").c_str() );


            std::string tCpuInfo = cpu_info();

            std::fprintf( stdout, "    CPU Info      : %s \n", clean_string( tCpuInfo.c_str() ).c_str() );
            std::fprintf( stdout, "    Procs Used    : %i \n", ( int ) gComm.size() );

            // insert blank line
            std::fprintf( stdout, "\n");

            // When built?
            std::fprintf( stdout, "    Build Date    : %s at %s\n", __DATE__ , __TIME__ );


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
#ifdef BELFEM_SUITESPARSE
            std::fprintf( stdout, "    Solver Lib    : Suite Sparse\n");
#endif
#ifdef BELFEM_MUMPS
            std::fprintf( stdout, "    Solver Lib    : MUMPS\n");
#endif
#ifdef BELFEM_STRUMPACK
            std::fprintf( stdout, "    Solver Lib    : STRUMPACK\n");
#endif
#ifdef BELFEM_PARDISO
            std::fprintf( stdout, "    Solver Lib    : PARDISO\n");
#endif
#ifdef BELFEM_PETSC
            std::fprintf( stdout, "    Solver Lib    : PETSc\n");
#endif
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

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_CLANG
#pragma clang diagnostic pop
#endif



//-----------------------------------------------------------------------------
}