//
// Created by Christian Messe on 18.11.18.
//

#ifndef BELFEM_STRINGTOOLS_HPP
#define BELFEM_STRINGTOOLS_HPP

#include <memory>
#include <iostream>
#include <string>
#include <cstdio>

#include "typedefs.hpp"
#include "cl_Cell.hpp"


namespace belfem
{
//------------------------------------------------------------------------------

// temporarily disable warnings


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
    string sprint( const char * aFormat, const Args ... aArgs )
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


#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#elif BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_INTEL
#pragma warning pop
#endif

//------------------------------------------------------------------------------

    template < typename T >
    inline std::string
    get_datatype_string( const T & aSample )
    {
        return "unknown";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const bool & aSample )
    {
        return "bool";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const int & aSample )
    {
        return "int";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const unsigned int & aSample )
    {
        return "unsigned int";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const long unsigned int & aSample )
    {
        return "long unsigned int";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const double & aSample )
    {
        return "double";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const std::string & aSample )
    {
        return "string";
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline std::string
    get_datatype_string( const std::complex<double> & aSample )
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
    first_word( const std::string & aString );

//------------------------------------------------------------------------------
    /**
     * create a cell of words from a string
     */
    Cell< string >
    string_to_words( const std::string & aString );

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
}
#endif //BELFEM_STRINGTOOLS_HPP