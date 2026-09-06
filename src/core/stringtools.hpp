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

#ifndef BELFEM_STRINGTOOLS_HPP
#define BELFEM_STRINGTOOLS_HPP

#include <memory>
#include <iostream>
#include <string>
#include <cstdio>
#include <regex>
#include <algorithm> // for std::all_of
#include <cctype>    // for isdigit

#include "typedefs.hpp"
#include "fn_sprint.hpp"

#include "cl_Cell.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

// temporarily disable warnings


#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#pragma clang diagnostic ignored "-Wformat"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"

#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#elif defined(BELFEM_INTEL)
    // Intel compiler diagnostics
#pragma warning(push)
    // disable format warnings
#pragma warning(disable: 1011)
    // disable unused variable warnings
#pragma warning(disable: 177)
    // keep your original one
#pragma warning(disable: 1595)

#endif

//------------------------------------------------------------------------------

    template < typename T >
    inline std::string
    datatype_string()
    {
        return "unknown";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<bool>()
    {
        return "bool";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<int>()
    {
        return "int";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<unsigned int>()
    {
        return "unsigned int";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<long unsigned int>()
    {
        return "long unsigned int";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<double>()
    {
        return "double";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<std::string>()
    {
        return "string";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    datatype_string<std::complex<double>>()
    {
        return "complex<double>";
    }

//------------------------------------------------------------------------------

    /**
     * returns the basename of a file
     */
    std::string
    basename( const std::string & aFilePath );

//------------------------------------------------------------------------------

    /**
     * returns the directory name of a file
     */
    std::string
    dirname( const std::string & aFilePath );

//------------------------------------------------------------------------------

    /**
     * returns the filetype of a file
     */
    std::string
    filetype( const std::string & aFilePath );

//------------------------------------------------------------------------------

    /**
     * returns the name of a file without the path
     */
    std::string
    filename( const std::string & aFilePath );

//------------------------------------------------------------------------------

    /**
     * tidy up the string
     */
    std::string
    clean_string( const std::string & aString );

//------------------------------------------------------------------------------

    /**
     * return the first word of the string
     */
    string
    first_word( const std::string & aString, const char aDelimiter = ' ');

//------------------------------------------------------------------------------
    /**
     * create a cell of words from a string
     */
    Cell< string >
    string_to_words( const std::string & aString, const char aDelimiter = ' ' );

//------------------------------------------------------------------------------

    std::string
    search_and_replace(
            const std::string & aString,
            const std::string & aSearch,
            const std::string & aReplace );

//------------------------------------------------------------------------------

    /**
     * convert the string to lower case
     */
    std::string
    string_to_lower( const std::string & aString  );

//------------------------------------------------------------------------------

    /**
     * convert the string to upper case
     */
    std::string
    string_to_upper( const std::string & aString  );

//------------------------------------------------------------------------------

    /**
     * return true if string is either 1, on, true or yes
     */
    bool
    string_to_bool( const std::string & aString );

//------------------------------------------------------------------------------

    /**
     * convert string to real, return NAN if it is not real
     */
     real
     to_real( const std::string & aString );

//------------------------------------------------------------------------------

     value
     unit_to_si( const string & aString );

//------------------------------------------------------------------------------

    string
    format_with_leading_zeros( const uint aNumber );

//------------------------------------------------------------------------------

    bool
    is_integer( const string & aString );

//------------------------------------------------------------------------------

    size_t
    utf8_character_count(const string & aString ) ;

//------------------------------------------------------------------------------

    template < typename T >
    void
    string_to_cell( const string & aString, Cell< T > & aValues )
    {

        Cell< string > tWords = string_to_words(
                search_and_replace(
                        search_and_replace(
                                search_and_replace(
                                        search_and_replace( search_and_replace(
                                                clean_string( aString ),
                                                "}","" ), "{",""),
                                          " ", "" ), ";", " ; " ),
                ","," " ) );


        // count memory
        for( string tWord : tWords )
        {
            if( tWord.find(":") < tWord.size() )
            {
                Cell< string > tNumbers = string_to_words( search_and_replace( tWord, ":", " " ) );
                T tA = ( T ) std::stoll( tNumbers(0) );
                T tB = ( T ) std::stoll( tNumbers(1) );

                if( tB > tA )
                {
                    T tN = tB - tA + 1 ;

                    T tC = tA ;
                    for( T k=0; k<tN; ++k )
                    {
                        aValues.push( tC );
                        tC += 1 ;
                    }
                }
                else
                {
                    T tN = tA - tB + 1 ;

                    T tC = tA ;
                    for( T k=0; k<tN; ++k )
                    {
                        aValues.push( tC );
                        tC -= 1 ;
                    }
                }
            }
            else if( tWord == ";" )
            {
                break;
            }
            else
            {

                aValues.push( ( T ) std::stoll( tWord ) );

            }
        }
    }


    template < typename T >
    void
    to_pair( const std::string & aString, Cell< std::pair< std::string, T >> & aResult )
    {
        aResult.clear() ;

        // Matches one or more letters ([A-Za-z]+) followed immediately by one or more digits (\d+)
        // The outer parentheses create the full match; inner ones capture the symbol and the number separately.
        std::regex re(R"(([A-Za-z]+)(\d+))");

        // std::sregex_iterator walks through ALL non-overlapping matches in the string
        for (std::sregex_iterator it(aString.begin(), aString.end(), re);
             it != std::sregex_iterator(); ++it)
        {
            const std::smatch& match = *it;

            std::string label = match[1].str();   // group 1 = letters (e.g. "Pb", "Sn", "RRR")
            T           value = std::stoi(match[2].str()); // group 2 = digits converted to int

            aResult.push( { std::move( label ), value } );
        }
    }

#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#elif BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_INTEL
#pragma warning pop
#endif

//------------------------------------------------------------------------------
}
#endif //BELFEM_STRINGTOOLS_HPP
