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
#include <filesystem>
#include "cl_Communicator.hpp"
#include "globals.hpp"
#include "filetools.hpp"

// Externally Defined Global Communicator
extern belfem::Communicator gComm;

namespace belfem
{
//------------------------------------------------------------------------------

    bool
    file_exists( const std::string & aPath )
    {
        return std::filesystem::exists( std::filesystem::path ( aPath ) );
    }

//------------------------------------------------------------------------------

    std::string
    search_data_file( const std::string & aFile,
                      const std::string & aSubDirectory )
    {
        // the run directory wins over the shared database
        if( file_exists( aFile ) )
        {
            return aFile;
        }

        // without a root there is nothing to fall back to. Note this is the
        // global, not $BELFEM_DATA: an installed tree fills it from
        // BELFEM_INSTALL_DATADIR even when the environment is silent
        if( gBelfemDataPath.size() == 0 )
        {
            return aFile;
        }

        // the subdirectory carries no separator of its own, so the join is
        // explicit here. An empty subdirectory searches the root itself
        const std::string tBase = aSubDirectory.size() > 0 ?
                gBelfemDataPath + "/" + aSubDirectory : gBelfemDataPath;

        // same relative layout below the data directory
        std::string tCandidate = tBase + "/" + aFile;

        if( file_exists( tCandidate ) )
        {
            return tCandidate;
        }

        // by name alone, for a path that only exists on the machine the
        // input file came from. Guarded on the separator: without one,
        // the candidate above already was the name alone
        const std::size_t tSlash = aFile.find_last_of( '/' );

        if( tSlash != std::string::npos )
        {
            tCandidate = tBase + "/" + aFile.substr( tSlash + 1 );

            if( file_exists( tCandidate ) )
            {
                return tCandidate;
            }
        }

        // hand back the original, so the error names the file the user wrote
        return aFile;
    }

//------------------------------------------------------------------------------

    /**
     * this function takes a path and makes it parallel
     */
    std::string
    make_path_parallel( const std::string & aPath )
    {

        // test if running in parallel mode
        if ( gComm.size() > 1 )
        {
            // get file extesion
            std::string tFileExt = filetype( aPath );

            // get base path
            std::string tBasePath = aPath.substr( 0, aPath.find_last_of(".") );

            // add proc number to path
            std::string aParallelPath = tBasePath + "_"
                                        +  std::to_string( gComm.size() ) + "."
                                        +  std::to_string( gComm.rank() ) + "."
                                        + tFileExt;
            return aParallelPath;
        }
        else
        {
            // do not modify path
            return aPath;
        }
    }

//------------------------------------------------------------------------------
}