//
// Created by Christian Messe on 2019-01-04.
//



#include "cl_Logger.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Logger::Logger( const uint aInfoLevel ) :
        mInfoLevel( aInfoLevel ),
        mStream( stdout )
    {

    }

//------------------------------------------------------------------------------

    Logger::Logger( const uint aInfoLevel, const std::string & aPath ) :
            mInfoLevel( aInfoLevel )
    {
        // create new file
        mStream = fopen( aPath.c_str(), "w+" );

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
}