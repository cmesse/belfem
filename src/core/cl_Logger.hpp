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
    enum class InfoLevel
    {
        Silent      = 0, // nothing
        Minimal     = 1, // minimal output
        Default     = 2, // default output
        Detailed    = 3, // a little more info about the mesh
        Verbose     = 4, // BELEM specific debugging info
        Everything  = 5  // also print debugging info from third party libraries
    };

//------------------------------------------------------------------------------

    /**
     * @brief Hierarchical logging with info levels.
     *
     * @ingroup grp_core
     * @see @ref core_core_usage_guide
     */
    class Logger
    {
//------------------------------------------------------------------------------

              uint       mInfoLevel = 0 ;
              std::FILE* mStream;
              bool       mWriteToAscii = false;

        // non-copyable, non-movable ( owns the FILE handle in ASCII mode )
        Logger( const Logger & ) = delete ;
        Logger & operator=( const Logger & ) = delete ;
        Logger( Logger && ) = delete ;
        Logger & operator=( Logger && ) = delete ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * legacy constructor using stdout
         */
        Logger( const uint aInfoLevel );

        /**
         * default constructor using stdout
         */
        Logger( const InfoLevel aInfoLevel );

//------------------------------------------------------------------------------

        /**
         *  constructor using file
         */
        Logger( const InfoLevel aInfoLevel, const std::string & aPath );

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

        /**
         * set the info level at runtime, e.g. from the command line
         */
        void
        set_info_level( const uint aInfoLevel );

        void
        set_info_level( const InfoLevel aInfoLevel );

//------------------------------------------------------------------------------

        template < typename ... Args >
        void
        message(
                const InfoLevel     aInfoLevel,
                const std::string & aFormat,
                const Args ...      aArgs )
        {
            if( static_cast< uint > ( aInfoLevel ) <= mInfoLevel )
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
    // Externally Defined Global Logger
    extern belfem::Logger gLog;

//------------------------------------------------------------------------------



    template < typename ... Args >
    void
    message(
            const belfem::InfoLevel     aInfoLevel,
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
