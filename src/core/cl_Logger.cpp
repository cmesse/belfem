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

#include "assert.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Logger::Logger( const uint aInfoLevel ) :
            mInfoLevel( aInfoLevel ),
            mStream( stdout )
    {

    }

    Logger::Logger( const InfoLevel aInfoLevel ) :
        mInfoLevel( static_cast< uint > ( aInfoLevel ) ),
        mStream( stdout )
    {

    }

//------------------------------------------------------------------------------

    Logger::Logger( const InfoLevel aInfoLevel, const std::string & aPath ) :
            mInfoLevel( static_cast< uint > ( aInfoLevel )  )
    {
        // create new file
        mStream = fopen( aPath.c_str(), "w+" );

        BELFEM_ERROR( mStream != nullptr,
            "Failed to open log file: %s", aPath.c_str() );

        // set flag that we use an ascii for the output
        mWriteToAscii = true;
    }

//------------------------------------------------------------------------------

    Logger::~Logger()
    {
        // check if logger runs in ASCII mode
        if( mWriteToAscii )
        {
            // close file
            fclose( mStream );
        }
    }

//------------------------------------------------------------------------------

    uint
    Logger::info_level() const
    {
        return mInfoLevel;
    }

//------------------------------------------------------------------------------

    void
    Logger::set_info_level( const uint aInfoLevel )
    {
        mInfoLevel = aInfoLevel;
    }

//------------------------------------------------------------------------------

    void
    Logger::set_info_level( const InfoLevel aInfoLevel )
    {
        mInfoLevel = static_cast< uint > ( aInfoLevel );
    }

//------------------------------------------------------------------------------
}