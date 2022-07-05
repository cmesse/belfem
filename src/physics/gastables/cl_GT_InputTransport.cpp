//
//
// Created by Christian Messe on 25.08.19.
//


#include "assert.hpp"
#include "cl_GT_InputTransport.hpp"
#include "fn_GT_fix_label.hpp"

namespace belfem
{
    namespace gastables
    {
//----------------------------------------------------------------------------

        InputTransport::InputTransport( const string & aPath )
                : Ascii( aPath, FileMode::OPEN_RDONLY )
        {
            uint tLineCount = 0;

            for ( std::string tLine : mBuffer )
            {
                // get first character of line
                char tChar = tLine.substr( 0, 1 ).c_str()[ 0 ];

                // first character must neither be space or exclamation mark
                if ( tChar != 32 && tChar != 45 && tLine.size() > 32 )
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
                    return 0;
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
            if ( tLineCount > 0 )
            {


                // get line
                const string & tLine = mBuffer( tLineCount++ );

                // get number of viscosity polynomials
                uint tNumberOfViscosities = std::stoi( tLine.substr( 35, 1 ) );

                // read viscosities
                for( uint k=0; k<tNumberOfViscosities; ++k )
                {
                    aRefgas->add_transport_poly( this->read_polynomial( tLineCount ) );
                }

                // get number of conductivity polynomials
                uint tNumberOfConductivities = std::stoi( tLine.substr( 37, 1 ) );

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
            // get line
            const std::string & tLine = mBuffer( aLineCount++ );

            // read minimum temperature
            real tTmin =  std::stod( tLine.substr( 2, 7 ) );

            // read maximum temperature
            real tTmax =  std::stod( tLine.substr( 11, 7 ) );

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

            // fix D letter
            tWord[ 11 ] = 'e';

            // fix plus sign
            if ( tWord[ 12 ] == '+' || tWord[ 12 ] == ' ' )
            {
                tWord[ 12 ] = '0';
            }

            return std::stod( tWord );
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */