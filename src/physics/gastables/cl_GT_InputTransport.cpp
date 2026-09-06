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

#include "assert.hpp"
#include "cl_GT_InputTransport.hpp"
#include "fn_GT_fix_label.hpp"
#include "fn_GT_parse.hpp"

namespace belfem
{
    namespace gastables
    {
//----------------------------------------------------------------------------

        InputTransport::InputTransport( const string & aPath )
                : Ascii( aPath, FileMode::OPEN_RDONLY )
        {
            uint tLineCount = 0;

            for ( const std::string & tLine : mBuffer )
            {
                // get first character of line
                char tChar = tLine.substr( 0, 1 ).c_str()[ 0 ];

                // first character must be neither space, exclamation mark nor dash
                if ( tChar != 32 && tChar != 33 && tChar != 45 && tLine.size() > 32 )
                {
                    // get label of this gas
                    string tFirst = fix_label( clean_string( tLine.substr( 0, 15 ) ) );

                    // get name of interaction parameter, if exists
                    string tSecond = fix_label( clean_string( tLine.substr( 16, 15 ) ) );

                    if ( tSecond.size() == 0 )
                    {

                        mMap[ tFirst ] = tLineCount;
                    }
                    else
                    {
                        // combine names
                        string tCombine1 = tFirst + "@" + tSecond;
                        string tCombine2 = tSecond + "@" + tFirst;

                        mMap[ tCombine1 ] = tLineCount;
                        mMap[ tCombine2 ] = tLineCount;
                    }
                }

                // increment counter
                ++tLineCount;
            }
        }
//------------------------------------------------------------------------------

        uint
        InputTransport::entry_exists( RefGas * aRefgas  )
        {
            if( mMap.key_exists( aRefgas->label() ) )
            {
                return  mMap( aRefgas->label() );
            }
            else
            {
                // check if there is an entry for the unionized species
                string tLabel = search_and_replace(
                        search_and_replace( aRefgas->label(), "+","" ),
                            "-","" );

                if( mMap.key_exists( tLabel ) )
                {
                    return mMap( tLabel );
                }
                else
                {
                    return BELFEM_UINT_MAX;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        InputTransport::read_data( RefGas * aRefgas )
        {
            // find line in file
            uint tLineCount = this->entry_exists( aRefgas );

            // test if entry exists
            if ( tLineCount != BELFEM_UINT_MAX )
            {
                BELFEM_ERROR( tLineCount < mBuffer.size(),
                        "trans.inp truncated at record of %s",
                        aRefgas->label().c_str() );

                // get line
                const string & tLine = mBuffer( tLineCount++ );

                BELFEM_ERROR( tLine.length() >= 38,
                        "Record line of %s in trans.inp is too short",
                        aRefgas->label().c_str() );

                // get number of viscosity polynomials
                uint tNumberOfViscosities = parse_int( tLine.substr( 35, 1 ),
                        "viscosity count in trans.inp" );

                // read viscosities
                for( uint k=0; k<tNumberOfViscosities; ++k )
                {
                    aRefgas->add_transport_poly( this->read_polynomial( tLineCount ) );
                }

                // get number of conductivity polynomials
                uint tNumberOfConductivities = parse_int( tLine.substr( 37, 1 ),
                        "conductivity count in trans.inp" );

                // read conductivities
                for( uint k=0; k<tNumberOfConductivities; ++k )
                {
                    aRefgas->add_transport_poly( this->read_polynomial( tLineCount ) );
                }
            }
        }

//------------------------------------------------------------------------------

        bool
        InputTransport::interaction_parameter_exists(
                const string & aLabelA,
                const string & aLabelB )
        {
            return mMap.key_exists( aLabelA + "@" + aLabelB );
        }

//------------------------------------------------------------------------------

        TransportPoly *
        InputTransport::read_polynomial( uint & aLineCount )
        {
            BELFEM_ERROR( aLineCount < mBuffer.size(),
                    "trans.inp truncated inside a polynomial block" );

            // get line
            const std::string & tLine = mBuffer( aLineCount++ );

            // fields reach to column 79
            BELFEM_ERROR( tLine.length() >= 80,
                    "Polynomial line in trans.inp is too short: '%s'",
                    tLine.c_str() );

            // read minimum temperature ( columns 2-10 )
            real tTmin = parse_real( tLine.substr( 2, 9 ), "Tmin in trans.inp" );

            // read maximum temperature ( columns 11-19 )
            real tTmax = parse_real( tLine.substr( 11, 9 ), "Tmax in trans.inp" );

            // create the coefficient vector
            Vector< real > tCoefficients( 4 );

            // populate the coefficient vector
            tCoefficients( 0 ) = this->word_to_real( tLine.substr( 20, 15 ) );
            tCoefficients( 1 ) = this->word_to_real( tLine.substr( 35, 15 ) );
            tCoefficients( 2 ) = this->word_to_real( tLine.substr( 50, 15 ) );
            tCoefficients( 3 ) = this->word_to_real( tLine.substr( 65, 15 ) );

            // the type of the polynomial
            auto tType = TransportPolyType::UNDEFINED;

            if( tLine.substr( 1, 1 ) == "C" )
            {
                tType = TransportPolyType::CONDUCTIVITY;
            }
            else if ( tLine.substr( 1, 1 ) == "V" )
            {
                tType = TransportPolyType::VISCOSITY;
            }
            else
            {
                BELFEM_ERROR( false, "Something went wrong while loading transport coefficients." );
            }

            // create the new object
            return new TransportPoly( tType, tTmin, tTmax, tCoefficients ) ;
        }


//------------------------------------------------------------------------------

        real
        InputTransport::word_to_real( const string & aWord )
        {
            string tWord( aWord );

            BELFEM_ERROR( tWord.length() >= 13,
                    "Coefficient field '%s' in trans.inp is too short",
                    aWord.c_str() );

            BELFEM_ERROR( tWord[ 11 ] == 'E' || tWord[ 11 ] == 'D'
                       || tWord[ 11 ] == 'e' || tWord[ 11 ] == 'd',
                    "Coefficient field '%s' in trans.inp is misaligned",
                    aWord.c_str() );

            // fix D letter
            tWord[ 11 ] = 'e';

            // fix plus sign
            if ( tWord[ 12 ] == '+' || tWord[ 12 ] == ' ' )
            {
                tWord[ 12 ] = '0';
            }

            return parse_real( tWord, "coefficient in trans.inp" );
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */