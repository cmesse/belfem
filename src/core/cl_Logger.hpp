//
// Created by Christian Messe on 19.11.18.
//

#ifndef BELFEM_CL_LOGGER_HPP
#define BELFEM_CL_LOGGER_HPP

#include <cstdio>
#include <fstream>

#include "stringtools.hpp"
#include "typedefs.hpp"

#ifdef BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#elif BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#endif

namespace belfem
{

//------------------------------------------------------------------------------

    class Logger
    {
//------------------------------------------------------------------------------

              uint       mInfoLevel;
              std::FILE* mStream;
              bool       mWriteToAscii = false;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * default constructor using stdout
         * todo : add description of levels
         */
        Logger( const uint aInfoLevel );

//------------------------------------------------------------------------------

        /**
         *  constructor using file
         */
        Logger( const uint aInfoLevel, const std::string & aPath );

//------------------------------------------------------------------------------

        /**
         * destructor
         */
         ~Logger();

//------------------------------------------------------------------------------

        /**
         * return the info level of the logger
         */
        uint
        info_level() const;

//------------------------------------------------------------------------------

        template < typename ... Args >
        void
        message(
                const uint          aInfoLevel,
                const std::string & aFormat,
                const Args ...      aArgs )
        {
            if( aInfoLevel <= mInfoLevel )
            {
                // format message and append a new line
                std::string tMessage = sprint(
                        aFormat.c_str(),
                        aArgs ... ) + "\n";

                std::fprintf( mStream, tMessage.c_str() );
            }
        }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
}
    // Externally Defined Global Communicator
    extern belfem::Logger gLog;

//------------------------------------------------------------------------------



    template < typename ... Args >
    void
    message(
            const belfem::uint     aInfoLevel,
            const std::string & aFormat,
            const Args ...      aArgs )
    {
        gLog.message( aInfoLevel, aFormat, aArgs ... );
    }

#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#elif BELFEM_GCC
#pragma GCC diagnostic pop
#endif

//------------------------------------------------------------------------------

#endif //BELFEM_CL_LOGGER_HPP
