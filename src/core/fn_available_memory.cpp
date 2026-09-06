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

#include "fn_available_memory.hpp"

#if defined( __linux__ )
#include <fstream>
#include <string>
#include <limits>
#include <cstdlib>
#include <cerrno>
#elif defined( __APPLE__ )
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <mach/vm_statistics.h>
#include <unistd.h>
#endif

namespace belfem
{
#if defined( __linux__ )
    namespace
    {
        // the tightest cgroup limit is ~2^63 bytes on v1 and the literal
        // "max" on v2; both mean "no limit" and both are folded to this
        constexpr unsigned long long gNoLimit
            = std::numeric_limits< unsigned long long >::max();

        // a decimal token to a number without exceptions: strtoull sets
        // errno on overflow and leaves the end pointer at the start on a
        // non-number, and this probe must never throw ( audit finding )
        bool
        parse_number( const std::string & aToken, unsigned long long & aValue )
        {
            if ( aToken.empty() || aToken[ 0 ] < '0' || aToken[ 0 ] > '9' )
            {
                return false ;
            }
            errno = 0 ;
            char * tEnd = nullptr ;
            const unsigned long long tValue
                = std::strtoull( aToken.c_str(), &tEnd, 10 );
            if ( errno != 0 || tEnd == aToken.c_str() )
            {
                return false ;
            }
            aValue = tValue ;
            return true ;
        }

        // one number from a one-line file; false when the file is absent
        // or does not start with a number. "max" ( cgroup v2 ) reads as
        // gNoLimit
        bool
        read_number( const std::string & aPath, unsigned long long & aValue )
        {
            std::ifstream tFile( aPath );
            if ( ! tFile.is_open() )
            {
                return false ;
            }
            std::string tToken ;
            tFile >> tToken ;
            if ( tToken == "max" )
            {
                aValue = gNoLimit ;
                return true ;
            }
            return parse_number( tToken, aValue );
        }

        // true if a cgroup v1 controller list ( "cpu,cpuacct,memory" )
        // names the memory controller; a plain substring test misses the
        // co-mounted forms ( audit finding )
        bool
        names_memory_controller( const std::string & aList )
        {
            std::size_t tStart = 0 ;
            while ( tStart <= aList.size() )
            {
                std::size_t tComma = aList.find( ',', tStart );
                if ( tComma == std::string::npos )
                {
                    tComma = aList.size();
                }
                if ( aList.compare( tStart, tComma - tStart, "memory" ) == 0 )
                {
                    return true ;
                }
                tStart = tComma + 1 ;
            }
            return false ;
        }

        // MemAvailable in bytes; 0 if the kernel does not report it
        unsigned long long
        mem_available()
        {
            std::ifstream tFile( "/proc/meminfo" );
            std::string tLine ;
            while ( std::getline( tFile, tLine ) )
            {
                if ( tLine.rfind( "MemAvailable:", 0 ) == 0 )
                {
                    // "MemAvailable:   43477368 kB"
                    std::size_t tPos = tLine.find_first_of( "0123456789" );
                    unsigned long long tKiloBytes = 0 ;
                    if (    tPos == std::string::npos
                         || ! parse_number( tLine.substr( tPos, tLine.find( ' ', tPos ) - tPos ),
                                            tKiloBytes ) )
                    {
                        return 0 ;
                    }
                    return tKiloBytes * 1024ULL ;
                }
            }
            return 0 ;
        }

        // headroom under the tightest cgroup memory limit on the path from
        // this process's cgroup up to the root; gNoLimit when none is set.
        // A limit already exceeded reports 1 byte, not 0, so the caller
        // cannot mistake "over the limit" for "unknown".
        //
        // aResolved: /proc/self/cgroup named a memory hierarchy AND at least
        // one level of it was readable below the standard mount point
        // ( /sys/fs/cgroup ). A hierarchy that is declared but nowhere to
        // be found ( a nonstandard mount ) is reported as unresolved so the
        // caller can refuse to fall back to the host figure -- the host
        // figure overstates what a limited job may take ( audit finding )
        unsigned long long
        cgroup_headroom( bool & aResolved )
        {
            unsigned long long tHeadroom = gNoLimit ;
            bool tDeclared = false ;
            aResolved = false ;

            std::ifstream tCgroup( "/proc/self/cgroup" );
            std::string tLine ;
            while ( std::getline( tCgroup, tLine ) )
            {
                // v2: "0::/a/b/c"      -> /sys/fs/cgroup/a/b/c/memory.max
                // v1: "N:memory:/a/b"  -> /sys/fs/cgroup/memory/a/b/memory.limit_in_bytes
                std::string tRoot ;
                std::string tLimitFile ;
                std::string tUsageFile ;
                std::string tPath ;

                if ( tLine.rfind( "0::", 0 ) == 0 )
                {
                    tRoot      = "/sys/fs/cgroup" ;
                    tLimitFile = "memory.max" ;
                    tUsageFile = "memory.current" ;
                    tPath      = tLine.substr( 3 );
                }
                else
                {
                    // "N:controller-list:/path"
                    std::size_t tFirst  = tLine.find( ':' );
                    std::size_t tSecond = tFirst == std::string::npos ?
                        std::string::npos : tLine.find( ':', tFirst + 1 );
                    if (    tSecond == std::string::npos
                         || ! names_memory_controller(
                                tLine.substr( tFirst + 1, tSecond - tFirst - 1 ) ) )
                    {
                        continue ;
                    }
                    tRoot      = "/sys/fs/cgroup/memory" ;
                    tLimitFile = "memory.limit_in_bytes" ;
                    tUsageFile = "memory.usage_in_bytes" ;
                    tPath      = tLine.substr( tSecond + 1 );
                }
                tDeclared = true ;

                // walk from the leaf to the root: a scheduler puts the
                // limit on the job's cgroup, which is an ancestor of the
                // step's
                while ( true )
                {
                    unsigned long long tLimit = gNoLimit ;
                    unsigned long long tUsage = 0 ;

                    if (    read_number( tRoot + tPath + "/" + tLimitFile, tLimit )
                         && read_number( tRoot + tPath + "/" + tUsageFile, tUsage ) )
                    {
                        aResolved = true ;
                        // v1 reports "unlimited" as a value near 2^63,
                        // page-rounded; anything that large is no limit
                        if ( tLimit < ( 1ULL << 62 ) )
                        {
                            unsigned long long tRoom
                                = tLimit > tUsage ? tLimit - tUsage : 1ULL ;
                            if ( tRoom < tHeadroom )
                            {
                                tHeadroom = tRoom ;
                            }
                        }
                    }

                    if ( tPath.empty() || tPath == "/" )
                    {
                        break ;
                    }
                    std::size_t tSlash = tPath.find_last_of( '/' );
                    tPath = tSlash == std::string::npos ? "" : tPath.substr( 0, tSlash );
                }
            }
            // no memory hierarchy declared at all: nothing to resolve, and
            // the host figure is the truth
            if ( ! tDeclared )
            {
                aResolved = true ;
            }
            return tHeadroom ;
        }
    }
#endif

//------------------------------------------------------------------------------

    std::size_t
    available_memory()
    {
#if defined( __linux__ )
        unsigned long long tHost = mem_available();
        if ( tHost == 0 )
        {
            return 0 ;
        }
        bool tResolved = false ;
        unsigned long long tRoom = cgroup_headroom( tResolved );
        if ( ! tResolved )
        {
            return 0 ;
        }
        return ( std::size_t ) ( tRoom < tHost ? tRoom : tHost );

#elif defined( __APPLE__ )
        mach_port_t tHost = mach_host_self();
        vm_statistics64_data_t tStat ;
        mach_msg_type_number_t tCount = HOST_VM_INFO64_COUNT ;
        kern_return_t tStatus = host_statistics64(
                tHost, HOST_VM_INFO64, ( host_info64_t ) &tStat, &tCount );

        // the send right from mach_host_self() is ours to release; the
        // probe runs once per workspace failure and would leak one per call
        mach_port_deallocate( mach_task_self(), tHost );

        if ( tStatus != KERN_SUCCESS )
        {
            return 0 ;
        }
        long tPage = sysconf( _SC_PAGESIZE );
        if ( tPage <= 0 )
        {
            return 0 ;
        }
        // free plus inactive: inactive pages are reclaimable, and "free"
        // alone understates what an allocation can claim on a warm machine
        return ( std::size_t ) ( ( unsigned long long ) tStat.free_count
                                + ( unsigned long long ) tStat.inactive_count )
             * ( std::size_t ) tPage ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------
}
