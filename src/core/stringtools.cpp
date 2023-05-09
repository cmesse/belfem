//
// Created by Christian Messe on 2018-12-23.
//

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
    first_word( const std::string & aString )
    {
        // tidy up string
        string tString = clean_string( aString );

        return tString.substr( 0, tString.find_first_of( " ", 0 ) );

    }

//------------------------------------------------------------------------------

    Cell< string >
    string_to_words( const std::string & aString )
    {
        // Cell of words
        Cell< string > aWords;

        // tidy up string
        string tString = clean_string( aString );

        // get length of cleaned string
        size_t tLength = tString.size();

        // extract words from string
        if( tLength > 1 )
        {
            size_t tEnd = 0;
            size_t tStart =  0;

            while( tEnd < tLength )
            {
                tEnd = tString.find_first_of( " ", tStart );
                if( tEnd > tLength )
                {
                    aWords.push(tString.substr(tStart, tLength - tStart));
                    break;
                }
                else
                {
                    aWords.push(tString.substr(tStart, tEnd - tStart));
                    tStart = tEnd + 1;
                }
            }
        }
        else if ( tLength == 1 )
        {
            if( tString != " " )
            {
                aWords.push( tString );
            }
        }

        return aWords;
    }

//------------------------------------------------------------------------------

    std::string
    search_and_replace(
            const std::string & aString,
            const std::string & aSearch,
            const std::string & aReplace )

    {
        std::string aOutString( aString );

        std::size_t tPos = 0;

        while( ( tPos = aOutString.find( aSearch, tPos) ) != std::string::npos )
        {
            aOutString.replace( tPos, aSearch.length(), aReplace );
            tPos += aReplace.length();
        }

        return aOutString;
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
        if( tEnd == aString.c_str() && *tEnd != '\0' )
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
                    tScale = 1e3;
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
                    tScale = 1e-3;
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

                // magnetic density
                else if ( tUnit == "G" )
                {
                    tScale = 1e-4;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "muT" )
                {
                    tScale = 1e-6;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "mT" )
                {
                    tScale = 1e-3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "T" )
                {
                    tScale = 1.0;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "kT" )
                {
                    tScale = 1e3;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }
                else if ( tUnit == "MT" )
                {
                    tScale = 1e6;
                    tMass += tPower;
                    tLength += tPower * 2;
                    tCurrent -= tPower;
                    tTime -= tPower * 3;
                }

                // temperature
                else if ( tUnit == "K" || tUnit == "°K" || tUnit == "C" || tUnit == "°C" )
                {
                    tScale = 1.0;
                    tTemperature += tPower;
                }
                else if ( tUnit == "F" || tUnit == "°F" || tUnit == "R" || tUnit == "°R" )
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
} /* namespace belfem */