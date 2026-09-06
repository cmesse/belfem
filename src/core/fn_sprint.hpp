/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_SPRINT_HPP
#define BELFEM_FN_SPRINT_HPP


#include <string>
#include <iostream>
#include <memory>

#ifdef BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#elif BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#elif BELFEM_INTEL
#pragma warning push
#pragma warning disable 1595
#endif

namespace belfem
{
    //------------------------------------------------------------------------------

    /**
     * A format script similar to write( *,* ) in fortran
     *
     * @tparam Args    type of arguments to be passed
     * @param aFormat  format string
     * @param aArgs    arguments
     * @return
     */
    template < typename ... Args >
    std::string sprint( const char * aFormat, const Args ... aArgs )
    {
        // Determine size of string.
        auto tSize = std::snprintf( nullptr, 0, aFormat, aArgs ... );

        // create unique pointer with length tSize. Add Extra space for '\0'.
        std::unique_ptr< char[] > tBuffer( new char[ tSize + 1 ] );

        // write formatted string into buffer
        std::snprintf( tBuffer.get(), tSize + 1, aFormat, aArgs ... );

        // return formatted string
        return string( tBuffer.get(), tBuffer.get() + tSize );
    }
}

#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#elif BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_INTEL
#pragma warning pop
#endif

#endif //BELFEM_FN_SPRINT_HPP