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

#ifndef BELFEM_FN_GT_PARSE_HPP
#define BELFEM_FN_GT_PARSE_HPP

#include <charconv>

#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        // strip leading/trailing whitespace and a leading '+'
        // ( from_chars accepts neither )
        inline std::string_view
        trim_field( const std::string_view aField )
        {
            size_t tBegin = aField.find_first_not_of( " \t\r\n" );
            if ( tBegin == std::string_view::npos )
            {
                return std::string_view();
            }
            size_t tEnd = aField.find_last_not_of( " \t\r\n" );
            std::string_view aView = aField.substr( tBegin, tEnd - tBegin + 1 );
            if ( ! aView.empty() && aView.front() == '+' )
            {
                aView.remove_prefix( 1 );
            }
            return aView;
        }

//------------------------------------------------------------------------------

        /**
         * locale-independent, non-throwing replacement for std::stod in the
         * fixed-column parsers; aborts with context on a malformed field
         */
        inline real
        parse_real( const std::string_view aField, const char * aContext )
        {
            std::string_view tView = trim_field( aField );

            real tValue;
            std::from_chars_result tResult =
                    std::from_chars( tView.data(), tView.data() + tView.size(), tValue );

            BELFEM_ERROR( tResult.ec == std::errc() && tResult.ptr == tView.data() + tView.size(),
                    "Malformed number '%.*s' while reading %s",
                    ( int ) aField.size(), aField.data(), aContext );

            return tValue;
        }

//------------------------------------------------------------------------------

        /**
         * like parse_real, but an empty or whitespace-only field
         * returns the given fallback instead of aborting
         */
        inline real
        parse_real_optional( const std::string_view aField, const real aFallback,
                             const char * aContext )
        {
            std::string_view tView = trim_field( aField );

            if ( tView.empty() )
            {
                return aFallback;
            }

            real tValue;
            std::from_chars_result tResult =
                    std::from_chars( tView.data(), tView.data() + tView.size(), tValue );

            BELFEM_ERROR( tResult.ec == std::errc() && tResult.ptr == tView.data() + tView.size(),
                    "Malformed number '%.*s' while reading %s",
                    ( int ) aField.size(), aField.data(), aContext );

            return tValue;
        }

//------------------------------------------------------------------------------

        /**
         * locale-independent, non-throwing replacement for std::stoi
         */
        inline int
        parse_int( const std::string_view aField, const char * aContext )
        {
            std::string_view tView = trim_field( aField );

            int aValue = 0;
            std::from_chars_result tResult =
                    std::from_chars( tView.data(), tView.data() + tView.size(), aValue );

            BELFEM_ERROR( tResult.ec == std::errc() && tResult.ptr == tView.data() + tView.size(),
                    "Malformed integer '%.*s' while reading %s",
                    ( int ) aField.size(), aField.data(), aContext );

            return aValue;
        }

//------------------------------------------------------------------------------

        /**
         * bounds-checked substring for fixed-column records; aborts with
         * context if the line is shorter than the requested window start
         */
        inline std::string
        parse_field( const std::string & aLine, const size_t aPos,
                     const size_t aCount, const char * aContext )
        {
            BELFEM_ERROR( aPos < aLine.length(),
                    "Line too short while reading %s : '%s'",
                    aContext, aLine.c_str() );

            return aLine.substr( aPos, aCount );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GT_PARSE_HPP
