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

#include "stringtools.hpp"
#include "assert.hpp"
#include "units.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    std::string
    basename( const std::string & aFilePath )
    {
        // find last entry of directory delimeter
        std::size_t tFound = aFilePath.find_last_of("/\\");

        // return basename
        return aFilePath.substr( tFound + 1 );
    }

//------------------------------------------------------------------------------

    std::string
    dirname( const std::string & aFilePath )
    {
        // find last entry of directory delimeter
        std::size_t tFound = aFilePath.find_last_of("/\\");

        // return dirname
        if( tFound < aFilePath.length() )
        {
            return aFilePath.substr( 0, tFound );
        }
        else
        {
            return "";
        }
    }

//------------------------------------------------------------------------------

    /**
     * returns the filetype of a file
     */
    std::string
    filetype( const std::string & aFilePath )
    {
        // find last entry of directory delimeter
        std::size_t tFound = aFilePath.find_last_of(".");

        // return basename
        return aFilePath.substr( tFound + 1 );
    }

//------------------------------------------------------------------------------

    /**
     * returns the name of a file without the path
     */
    std::string
    filename( const std::string & aFilePath )
    {
        // find last entry of directory delimeter
        std::size_t tFound = aFilePath.find_last_of("/");

        // return basename
        return aFilePath.substr( tFound + 1, aFilePath.size() - tFound - 1 );
    }


//------------------------------------------------------------------------------

    /**
     * tidy up the string
     */
    std::string
    clean_string( const std::string & aString )
    {
        // Quick check if cleaning is needed
        bool needs_cleaning = false;
        for (char c : aString) {
            if (c <= 32 || c == '"' || c == '#') {
                needs_cleaning = true;
                break;
            }
        }
        if (!needs_cleaning && aString.front() != ' ' && aString.back() != ' ') {
            return aString;
        }

        // length of string
        size_t tLength = aString.size();

        if( tLength > 1 )
        {
            // special characters
            const char tSpace = 32;
            const char tTab   =  9;
            const char tCr    = 13;
            const char tLf    = 10;
            const char tQuote = 34;
            const char tHash  = 35;

            // grab array of instring
            const char * tChars = aString.c_str();

            // create new string
            string aOutString;
            aOutString.reserve( tLength );

            bool tIsSpace    = true;
            bool tIsSpaceOld;

            bool tInQuote    = false;

            // loop over all chars
            for( size_t c=0; c < tLength; ++c )
            {
                // shift space flag
                tIsSpaceOld = tIsSpace;

                // generate new space flag
                tIsSpace = (    tChars[ c ] == tSpace
                             || tChars[ c ] == tTab
                             || tChars[ c ] == tCr
                             || tChars[ c ] == tLf );


                // test if flag is space
                if ( tIsSpace )
                {
                    // test if is in quote
                    if( tInQuote )
                    {
                        aOutString.append( aString.substr( c, 1 ) );
                    }
                    else if( ! tIsSpaceOld )
                    {
                        // add space to string
                        aOutString.append( " " );
                    }
                }
                else if ( tChars[ c ] == tQuote )
                {
                    // flip quote switch
                    tInQuote = ! tInQuote;

                    // add character to outstring
                    aOutString.append( aString.substr( c, 1 ) );
                }
                else if ( tChars[ c ] == tHash )
                {
                    break;
                }
                else
                {
                    // add character to outstring
                    aOutString.append( aString.substr( c, 1 ) );
                }
            }

            // treat last string
            if ( tIsSpace && aOutString.size() > 0 )
            {
                // erase last character
                aOutString.erase( aOutString.end() - 1 );
            }

            // output result
            return aOutString;
        }
        else
        {
            return aString;
        }
    }

//------------------------------------------------------------------------------

    string
    first_word( const std::string & aString, char aDelimiter )
    {
        // tidy up string
        string tString = clean_string( aString );

        return tString.substr( 0, tString.find_first_of( aDelimiter, 0 ) );

    }

//------------------------------------------------------------------------------

    Cell<string> string_to_words(const std::string& aString, char aDelimiter) {
        Cell<string> words;

        // Skip cleaning if delimiter is not space
        const std::string& str = (aDelimiter == ' ') ? clean_string(aString) : aString;

        if (str.empty()) return words;

        // Pre-count words for reserve
        size_t word_count = 1;
        for (char c : str) {
            if (c == aDelimiter) word_count++;
        }
        words.reserve(word_count);

        // Use string_view in C++17 or find_first_not_of for better performance
        size_t start = 0;
        size_t end = 0;

        while ((end = str.find(aDelimiter, start)) != std::string::npos) {
            if (end > start) {
                words.push(str.substr(start, end - start));
            }
            start = end + 1;
        }

        if (start < str.length()) {
            words.push(str.substr(start));
        }

        return words;
    }

//------------------------------------------------------------------------------

    std::string search_and_replace(
     const std::string& aString,
     const std::string& aSearch,
     const std::string& aReplace) {

        if (aSearch.empty()) return aString;

        // Quick check if search string exists
        if (aString.find(aSearch) == std::string::npos) {
            return aString;
        }

        std::string result;
        result.reserve(aString.length() * 1.5); // Reasonable estimate

        size_t lastPos = 0;
        size_t findPos = 0;

        while ((findPos = aString.find(aSearch, lastPos)) != std::string::npos) {
            result.append(aString, lastPos, findPos - lastPos);
            result.append(aReplace);
            lastPos = findPos + aSearch.length();
        }

        result.append(aString, lastPos, std::string::npos);
        return result;
    }

//------------------------------------------------------------------------------

    std::string
    string_to_lower( const std::string & aString  )
    {
        std::string aOutString( aString );

        std::transform( aOutString.begin(), aOutString.end(), aOutString.begin(),
            [](unsigned char c)
            {
                return std::tolower( c );
            }
        );

        return aOutString;
    }

//------------------------------------------------------------------------------

    std::string
    string_to_upper( const std::string & aString  )
    {
        std::string aOutString( aString );

        std::transform( aOutString.begin(), aOutString.end(), aOutString.begin(),
                        [](unsigned char c)
                        {
                            return std::toupper( c );
                        }
        );

        return aOutString;
    }

//------------------------------------------------------------------------------

    bool
    string_to_bool( const std::string & aString )
    {
        // lower string of aString
        std::string tLowerString( string_to_lower( aString ) );

        return (    tLowerString == "true"
                 || tLowerString == "on"
                 || tLowerString == "yes"
                 || tLowerString == "1" ) ;
    }


//------------------------------------------------------------------------------

    real
    to_real( const std::string & aString )
    {
        char * tEnd ;
        real aValue = std::strtod( aString.c_str(), & tEnd );

        // check if this value works
        if( tEnd == aString.c_str() )
        {
            return BELFEM_QUIET_NAN ;
        }
        else
        {
            return aValue ;
        }
    }

//------------------------------------------------------------------------------

    value
    unit_to_si( const string & aString )
    {
        value aValue;
        aValue.first = 1.0;
        aValue.second = { 0, 0, 0, 0, 0, 0, 0 };


        real & tLength = aValue.second[ 0 ];
        real & tMass = aValue.second[ 1 ];
        real & tTime = aValue.second[ 2 ];
        real & tCurrent = aValue.second[ 3 ];
        real & tTemperature = aValue.second[ 4 ];
        real & tSubstance = aValue.second[ 5 ];
        real & tBrightness = aValue.second[ 6 ];

        // remove spaces
        string tString = search_and_replace( aString, " ", "" );

        // no unit
        if ( tString == "-" || tString == "" )
        {
            return aValue;
        }

        // remove brackets
        tString = search_and_replace( tString, "(", "" );
        tString = search_and_replace( tString, ")", "" );

        // change multiplication with space
        tString = search_and_replace( tString, "*", " " );

        // mu symbol
        tString = search_and_replace( tString, "µ", "mu" );

        // exponents
        tString = search_and_replace( tString, "²", "^2" );
        tString = search_and_replace( tString, "³", "^3" );

        // numerator string
        size_t tSplit = tString.find( "/", 0 );

        Cell< string > tNumerators = string_to_words( tString.substr( 0, tSplit ));


        Cell< string > tDenominators;

        if ( tSplit < tString.length())
        {
            tDenominators = string_to_words( tString.substr( tSplit + 1, tString.length()));
        }


        for ( uint s = 0; s < 2; ++s )
        {
            Cell< string > & tUnits = s == 0 ? tNumerators : tDenominators;

            for ( string & tUnit: tUnits )
            {
                // real tValue ;
                int tPower = s == 0 ? 1.0 : -1.0;

                // get power
                tSplit = tUnit.find( "^", 0 );
                if ( tSplit < tUnit.length())
                {
                    tPower *= std::stod( tUnit.substr( tSplit + 1, tUnit.length()));
                    tUnit = tUnit.substr( 0, tSplit );
                }

                // scale factor for this unit
                real tScale = BELFEM_QUIET_NAN;

                if ( tUnit == "nm" )
                {
                    tScale = 1e-9;
                    tLength += tPower;
                }
                else if ( tUnit == "mum" )
                {
                    tScale = 1e-6;
                    tLength += tPower;
                }
                else if ( tUnit == "mm" )
                {
                    tScale = 1e-3;
                    tLength += tPower;
                }
                else if ( tUnit == "cm" )
                {
                    tScale = 1e-2;
                    tLength += tPower;
                }
                else if ( tUnit == "dm" )
                {
                    tScale = 1e-1;
                    tLength += tPower;
                }
                else if ( tUnit == "m" )
                {
                    tScale = 1.0;
                    tLength += tPower;
                }
                else if ( tUnit == "km" )
                {
                    tScale = 1e3;
                    tLength += tPower;
                }
                else if ( tUnit == "Mm" )
                {
                    tScale = 1e6;
                    tLength += tPower;
                }
                else if ( tUnit == "in" || tUnit == "inch" )
                {
                    tScale = constant::in;
                    tLength += tPower;
                }
                else if ( tUnit == "ft" )
                {
                    tScale = constant::ft;
                    tLength += tPower;
                }
                else if ( tUnit == "kft" )
                {
                    tScale = 1000 * constant::ft;
                    tLength += tPower;
                }
                else if ( tUnit == "mi" )
                {
                    tScale = constant::mi;
                    tLength += tPower;
                }
                else if ( tUnit == "gal" )
                {
                    tScale = constant::gal;
                    tLength += tPower * 3;
                }
                else if ( tUnit == "oz" )
                {
                    tScale = constant::oz;
                    tLength += tPower * 3;
                }
                else if ( tUnit == "fur" )
                {
                    tScale = 660 * constant::ft;
                    tLength += tPower;
                }

                // mass
                else if ( tUnit == "mg" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                }
                else if ( tUnit == "g" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                }
                else if ( tUnit == "kg" )
                {
                    tScale = 1;
                    tMass += tPower;
                }
                else if ( tUnit == "t" )
                {
                    tScale = 1e3;
                    tMass += tPower;

                }
                else if ( tUnit == "lb" )
                {
                    tScale = constant::lb;
                    tMass += tPower;
                }
                else if ( tUnit == "slug" )
                {
                    tScale = constant::lbf / constant::ft;
                    tMass += tPower;
                }
                else if ( tUnit == "fir" )
                {
                    tScale = 90 * constant::lb;
                    tMass += tPower;
                }

                // time
                else if ( tUnit == "ns" )
                {
                    tScale = 1e-9;
                    tTime += tPower;
                }
                else if ( tUnit == "mus" )
                {
                    tScale = 1e-6;
                    tTime += tPower;
                }
                else if ( tUnit == "ms" )
                {
                    tScale = 1e-3;
                    tTime += tPower;
                }
                else if ( tUnit == "s" )
                {
                    tScale = 1.0;
                    tTime += tPower;
                }
                else if ( tUnit == "min" )
                {
                    tScale = 60.0;
                    tTime += tPower;
                }
                else if ( tUnit == "h" )
                {
                    tScale = 3600.0;
                    tTime += tPower;
                }
                else if ( tUnit == "Hz" )
                {
                    tScale = 1.0;
                    tTime -= tPower;
                }
                else if ( tUnit == "kHz" )
                {
                    tScale = 1e3 ;
                    tTime -= tPower;
                }
                else if ( tUnit == "MHz" )
                {
                    tScale = 1e6 ;
                    tTime -= tPower;
                }
                else if ( tUnit == "ftn" )
                {
                    tScale = 1209600;
                    tTime += tPower;
                }

                // force
                else if ( tUnit == "nN" )
                {
                    tScale = 1e-9;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "muN" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "mN" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "N" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kN" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "MN" )
                {
                    tScale = 1e6;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "mlbf" )
                {
                    tScale = constant::lbf * 1e-3;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "lbf" )
                {
                    tScale = constant::lbf;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "klbf" )
                {
                    tScale = constant::lbf * 1e3;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "Mlbf" )
                {
                    tScale = constant::lbf * 1e6;
                    tMass += tPower;
                    tLength += tPower;
                    tTime -= tPower * 2;
                }

                // energy
                else if ( tUnit == "mJ" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "J" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kJ" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "MJ" )
                {
                    tScale = 1e6;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "cal" )
                {
                    tScale = constant::calTh;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kcal" )
                {
                    tScale = constant::calTh * 1e3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 2;
                }

                // power
                else if ( tUnit == "nW" )
                {
                    tScale = 1e-9;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "muW" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "mW" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "W" )
                {
                    tScale = 1;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "kW" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "MW" )
                {
                    tScale = 1e6;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "GW" ) // what the hell is a jigawatt?
                {
                    tScale = 1e9;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tTime -= tPower * 3;
                }

                // pressure
                else if ( tUnit == "muPa" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "mPa" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "Pa" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kPa" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "MPa" )
                {
                    tScale = 1e6;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "GPa" )
                {
                    tScale = 1e9;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "mbar" )
                {
                    tScale = 1e2;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kbar" )
                {
                    tScale = 1e8;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "atm" )
                {
                    tScale = 101325;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }

                else if ( tUnit == "psi" )
                {
                    tScale = constant::psi;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kpsi" )
                {
                    tScale = constant::psi * 1e3;
                    tMass += tPower;
                    tLength -= tPower;
                    tTime -= tPower * 2;
                }

                // current
                else if ( tUnit == "mA" )
                {
                    tScale = 1e-3;
                    tCurrent += tPower;
                }
                else if ( tUnit == "A" )
                {
                    tScale = 1.0;
                    tCurrent += tPower;
                }
                else if ( tUnit == "kA" )
                {
                    tScale = 1e3;
                    tCurrent += tPower;
                }


                // voltage
                else if ( tUnit == "muV" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "mV" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "V" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "kV" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }

                // magnetic density: tesla is V*s/m^2 = kg/( s^2 * A ), so mass
                // +1, current -1, time -2 and NO length term. Until 2026-08-30
                // every entry here carried the VOLT exponents instead ( mass +1,
                // length +2, current -1, time -3 ) -- the "V" block above copied
                // with only the string changed. check_unit compares nothing but
                // these seven exponents, so a voltage was accepted wherever a
                // flux density was required and a correctly written V*s/m^2 was
                // refused. The scale factors were never wrong.
                else if ( tUnit == "G" )
                {
                    tScale = 1e-4;
                    tMass += tPower;
                    tCurrent -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "muT" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tCurrent -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "mT" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tCurrent -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "T" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tCurrent -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "kT" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tCurrent -= tPower;
                    tTime -= tPower * 2;
                }
                else if ( tUnit == "MT" )
                {
                    tScale = 1e6;
                    tMass += tPower;
                    tCurrent -= tPower;
                    tTime -= tPower * 2;
                }

                // temperature
                else if ( tUnit == "K" || tUnit == "°K" || tUnit == "C" || tUnit == "°C" )
                {
                    tScale = 1.0;
                    tTemperature += tPower;
                }
                else if (  tUnit == "°F" || tUnit == "R" || tUnit == "°R" )
                {
                    tScale = 5.0 / 9.0;
                    tTemperature += tPower;
                }

                // substance
                else if ( tUnit == "mol" )
                {
                    tScale = 1.0;
                    tSubstance += tPower;
                }
                else if ( tUnit == "kmol" )
                {
                    tScale = 1e3;
                    tSubstance += tPower;
                }

                // brightness
                else if ( tUnit == "cd" )
                {
                    tScale = 1.0;
                    tBrightness += tPower;
                }

                // angle
                else if ( tUnit == "rad" )
                {
                    tScale = 1.0 ;
                }
                else if ( tUnit == "°" || tUnit == "deg" )
                {
                    tScale = constant::pi / 180.0 ;
                }

                // brightness
                else if ( tUnit == "cd" )
                {
                    tScale = 1.0;
                    tBrightness += tPower;
                }

                // resistance
                else if ( tUnit == "nOhm" || tUnit == "nΩ" )
                {
                    tScale = 1e-9;
                    tMass += tPower;
                    tLength += tPower * 2.0;
                    tTime -= tPower * 3;
                    tCurrent -= tPower * 2;
                }
                else if ( tUnit == "muOhm" || tUnit == "muΩ" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tLength += tPower * 2.0;
                    tTime -= tPower * 3;
                    tCurrent -= tPower * 2;
                }
                else if ( tUnit == "mOhm" || tUnit == "mΩ" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength += tPower * 2.0;
                    tTime -= tPower * 3;
                    tCurrent -= tPower * 2;
                }
                else if ( tUnit == "Ohm" || tUnit == "Ω" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tLength += tPower * 2.0;
                    tTime -= tPower * 3;
                    tCurrent -= tPower * 2;
                }
                else if ( tUnit == "kOhm" || tUnit == "kΩ" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength += tPower * 2.0;
                    tTime -= tPower * 3;
                    tCurrent -= tPower * 2;
                }
                //capacitance
                else if ( tUnit == "nF" )
                {
                    tScale = 1e-9;
                    tCurrent += tPower * 2 ;
                    tTime += tPower * 4 ;
                    tLength -= tPower * 2 ;
                    tMass -= tPower ;
                }
                else if ( tUnit == "muF" )
                {
                    tScale = 1e-6;
                    tCurrent += tPower * 2 ;
                    tTime += tPower * 4 ;
                    tLength -= tPower * 2 ;
                    tMass -= tPower ;
                }
                else if ( tUnit == "mF" )
                {
                    tScale = 1e-3;
                    tCurrent += tPower * 2 ;
                    tTime += tPower * 4 ;
                    tLength -= tPower * 2 ;
                    tMass -= tPower ;
                }
                else if ( tUnit == "F" )
                {
                    tScale = 1.0;
                    tCurrent += tPower * 2 ;
                    tTime += tPower * 4 ;
                    tLength -= tPower * 2 ;
                    tMass -= tPower ;
                }
                // inductance
                else if ( tUnit == "nH" )
                {
                    tScale = 1e-9;
                    tCurrent -= tPower * 2 ;
                    tTime -= tPower * 2 ;
                    tLength += tPower * 2 ;
                    tMass += tPower ;
                }
                else if ( tUnit == "muH" )
                {
                    tScale = 1e-6;
                    tCurrent -= tPower * 2 ;
                    tTime -= tPower * 2 ;
                    tLength += tPower * 2 ;
                    tMass += tPower ;
                }
                else if ( tUnit == "mH" )
                {
                    tScale = 1e-3;
                    tCurrent -= tPower * 2 ;
                    tTime -= tPower * 2 ;
                    tLength += tPower * 2 ;
                    tMass += tPower ;
                }
                else if ( tUnit == "H" )
                {
                    tScale = 1.0;
                    tCurrent -= tPower * 2 ;
                    tTime -= tPower * 2 ;
                    tLength += tPower * 2 ;
                    tMass += tPower ;
                }
                // magfield
                else if ( tUnit == "muOe" )
                {
                    tScale = 250.e-6 / constant::pi ;
                    tCurrent += tPower ;
                    tLength  -= tPower ;
                }
                else if ( tUnit == "mOe" )
                {
                    tScale = 250.e-3 / constant::pi ;
                    tCurrent += tPower ;
                    tLength  -= tPower ;
                }
                else if ( tUnit == "Oe" )
                {
                    tScale = 250.0 / constant::pi ;
                    tCurrent += tPower ;
                    tLength  -= tPower ;
                }
                else if ( tUnit == "kOe" )
                {
                    tScale = 250e3 / constant::pi ;
                    tCurrent += tPower ;
                    tLength  -= tPower ;
                }
                else if ( tUnit == "MOe" )
                {
                    tScale = 250e6 / constant::pi ;
                    tCurrent += tPower ;
                    tLength  -= tPower ;
                }
                // unknown
                else
                {
                    BELFEM_ERROR( false, "Unknown unit : %s", tString.c_str());
                    aValue.first = BELFEM_QUIET_NAN;
                    return aValue;
                }

                aValue.first *= std::pow( tScale, tPower );
            }

        }

        return aValue;
    }
//------------------------------------------------------------------------------

    string
    format_with_leading_zeros(const uint aNumber)
    {
        if (aNumber < 10)
        {
            return "%01u";
        }
        else if (aNumber < 100)
        {
            return "%02u";
        }
        else if (aNumber < 1000)
        {
            return "%03u";
        }
        else if (aNumber < 10000)
        {
            return "%04u";
        }
        else if (aNumber < 100000)
        {
            return "%05u";
        }
        else if (aNumber < 1000000)
        {
            return "%06u";
        }
        else if (aNumber < 10000000)
        {
            return "%07u";
        }
        else if (aNumber < 100000000)
        {
            return "%08u";
        }
        else if (aNumber < 1000000000)
        {
            return "%09u";
        }
        else
        {
            return "%10u";
        }
    }

//------------------------------------------------------------------------------

    bool is_integer( const string & aString )
    {
        // empty string is not an integer
        if ( aString.empty() ) return false;

        // check for optional sign and adjust starting position
        std::size_t tStart = aString[0] == '+' || aString[0] == '-' ? 1 : 0;

        // If the string is just a sign with no digits, it's not an integer
        if ( tStart == aString.size() ) return false;

        // check that all remaining characters are digits
        return std::all_of( aString.begin() + tStart, aString.end(), ::isdigit);
    }

//------------------------------------------------------------------------------

    size_t
    utf8_character_count(const string & aString )
    {
        size_t count = 0;
        for (size_t i = 0; i < aString.size();)
        {
            unsigned char c = aString[i];
            if (c <= 0x7F) i += 1;        // 1-byte (ASCII)
            else if (c <= 0xDF) i += 2;   // 2-byte
            else if (c <= 0xEF) i += 3;   // 3-byte
            else i += 4;                  // 4-byte
            count++;
        }
        return count;
    }

//------------------------------------------------------------------------------
} /* namespace belfem */
