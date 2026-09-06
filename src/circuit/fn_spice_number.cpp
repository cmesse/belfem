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

#include <cctype>
#include <charconv>
#include <cmath>
#include <limits>

#include "fn_spice_number.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        real
        spice_number_to_si( const string & aToken )
        {
            const size_t tLength = aToken.length();
            size_t tPos = 0;

            // optional sign
            if ( tPos < tLength &&
                 ( aToken[ tPos ] == '+' || aToken[ tPos ] == '-' ) )
            {
                ++tPos;
            }

            // digits before and after the decimal point; the grammar demands
            // at least one digit in total ( "1", "1.", ".5" are all legal )
            size_t tNumDigits = 0;
            while ( tPos < tLength &&
                    std::isdigit( static_cast< unsigned char >( aToken[ tPos ] ) ) )
            {
                ++tPos;
                ++tNumDigits;
            }
            if ( tPos < tLength && aToken[ tPos ] == '.' )
            {
                ++tPos;
                while ( tPos < tLength &&
                        std::isdigit( static_cast< unsigned char >( aToken[ tPos ] ) ) )
                {
                    ++tPos;
                    ++tNumDigits;
                }
            }

            BELFEM_ERROR( tNumDigits > 0,
                          "not a SPICE number: '%s'",
                          aToken.c_str() );

            // optional integer exponent, consumed only when complete --
            // in "1e" or "1exact" the 'e' is a unit letter, not an exponent
            if ( tPos < tLength &&
                 ( aToken[ tPos ] == 'e' || aToken[ tPos ] == 'E' ) )
            {
                size_t tExpPos = tPos + 1;
                if ( tExpPos < tLength &&
                     ( aToken[ tExpPos ] == '+' || aToken[ tExpPos ] == '-' ) )
                {
                    ++tExpPos;
                }
                size_t tExpDigits = 0;
                while ( tExpPos < tLength &&
                        std::isdigit( static_cast< unsigned char >( aToken[ tExpPos ] ) ) )
                {
                    ++tExpPos;
                    ++tExpDigits;
                }
                if ( tExpDigits > 0 )
                {
                    tPos = tExpPos;
                }
            }

            // convert only the validated prefix with the locale-independent
            // std::from_chars ( precedent: gastables fn_GT_parse.hpp ), so
            // neither a non-C LC_NUMERIC nor strtod's extra grammars ( inf,
            // nan, hex floats ) can be reached. from_chars accepts a leading
            // '-' but not a leading '+', which the scanner already consumed
            const size_t tNumBegin = ( aToken[ 0 ] == '+' ) ? 1 : 0;
            real tValue = 0.0;
            const std::from_chars_result tParse =
                    std::from_chars( aToken.data() + tNumBegin,
                                     aToken.data() + tPos,
                                     tValue );

            BELFEM_ERROR( tParse.ec != std::errc::result_out_of_range,
                          "SPICE number out of range: '%s'",
                          aToken.c_str() );

            // policy, pinned toolchain-independently: accepted magnitudes
            // are zero or the normal double range -- subnormals are
            // rejected loudly ( libstdc++ 11 already reports them as
            // out-of-range, newer versions would accept them silently )
            BELFEM_ERROR( tValue == 0.0 ||
                          std::fabs( tValue ) >= std::numeric_limits< real >::min(),
                          "SPICE number out of range: '%s'",
                          aToken.c_str() );

            // the scanner validated the grammar, so the only consume
            // mismatch from_chars may leave behind is a bare trailing
            // decimal point ( "1." )
            const char * tPrefixEnd = aToken.data() + tPos;
            BELFEM_ERROR( tParse.ec == std::errc() &&
                          ( tParse.ptr == tPrefixEnd ||
                            ( tPrefixEnd - tParse.ptr == 1 && *tParse.ptr == '.' ) ),
                          "not a SPICE number: '%s'",
                          aToken.c_str() );

            // scale factor, case-insensitive, longest match first
            // ( ngspice manual, "Ngspice scale factors" )
            real tScale = 1.0;

            auto tFoldedMatch = [ &aToken, tLength ]( const size_t aPos,
                                                      const char * aFactor,
                                                      const size_t aFactorLength ) -> bool
            {
                if ( aPos + aFactorLength > tLength )
                {
                    return false;
                }
                for ( size_t k = 0; k < aFactorLength; ++k )
                {
                    if ( std::tolower( static_cast< unsigned char >( aToken[ aPos + k ] ) )
                         != aFactor[ k ] )
                    {
                        return false;
                    }
                }
                return true;
            };

            if ( tFoldedMatch( tPos, "meg", 3 ) )
            {
                tScale = 1.0e6;
                tPos += 3;
            }
            else if ( tFoldedMatch( tPos, "mil", 3 ) )
            {
                tScale = 25.4e-6;
                tPos += 3;
            }
            else if ( tPos < tLength )
            {
                switch ( std::tolower( static_cast< unsigned char >( aToken[ tPos ] ) ) )
                {
                    case 't' : { tScale = 1.0e12;  ++tPos; break; }
                    case 'g' : { tScale = 1.0e9;   ++tPos; break; }
                    case 'k' : { tScale = 1.0e3;   ++tPos; break; }
                    case 'm' : { tScale = 1.0e-3;  ++tPos; break; }
                    case 'u' : { tScale = 1.0e-6;  ++tPos; break; }
                    case 'n' : { tScale = 1.0e-9;  ++tPos; break; }
                    case 'p' : { tScale = 1.0e-12; ++tPos; break; }
                    case 'f' : { tScale = 1.0e-15; ++tPos; break; }
                    case 'a' : { tScale = 1.0e-18; ++tPos; break; }
                    default  : { break; }
                }
            }

            // letters after the number or the scale factor are ignored
            // ( 10Volts, 1kOhm, 1MSec ); anything else is a hard error, which
            // also rejects the old-style embedded suffix "2k5"
            while ( tPos < tLength &&
                    std::isalpha( static_cast< unsigned char >( aToken[ tPos ] ) ) )
            {
                ++tPos;
            }

            BELFEM_ERROR( tPos == tLength,
                          "unexpected character '%c' after SPICE number: '%s' "
                          "( an embedded suffix like '2k5' is not supported, write '2.5k' )",
                          aToken[ tPos ],
                          aToken.c_str() );

            // the scale factor can overflow a finite mantissa ( "1e300T" )
            // or shrink it into the subnormal range ( "1e-300f" ) -- the
            // same magnitude policy applies to the product
            const real tResult = tValue * tScale;

            BELFEM_ERROR( std::isfinite( tResult ) &&
                          ( tResult == 0.0 ||
                            std::fabs( tResult ) >= std::numeric_limits< real >::min() ),
                          "SPICE number out of range: '%s'",
                          aToken.c_str() );

            return tResult;
        }

//-----------------------------------------------------------------------------
    }
}
