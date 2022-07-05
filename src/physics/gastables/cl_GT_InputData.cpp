//
// Created by Christian Messe on 26.08.19.
//


#include <cmath>
#include "stringtools.hpp"
#include "cl_GT_InputData.hpp"


namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        InputData::InputData( const string & aPath )
            : Ascii(aPath, FileMode::OPEN_RDONLY )
        {
            // init counter
            uint tLineCount = 0;

            // loop over all lines
            for ( std::string tLine : mBuffer )
            {
                // test if this line is commented out
                if( tLine.at( 0 ) != '!' )
                {
                    // position of first space
                    uint tFirst = tLine.find_first_of( ' ' );

                    if ( tFirst > 0 && tFirst < tLine.length())
                    {
                        // get first word of file
                        std::string tLabel = tLine.substr( 0, tFirst );

                        // add label to map
                        mMap[ tLabel ] = tLineCount;
                    }
                }
                // increment line counter
                tLineCount++;
            }
        }

//------------------------------------------------------------------------------

        void
        InputData::read_data( GasData * aData )
        {
            if ( this->entry_exists( aData ) )
            {
                // get line
                std::string & tLine = mBuffer( mMap( aData->label() ) );

                Vector< uint > tStart = { 92, 102, 112, 121, 129, 138 };
                Vector< uint > tLength = { 10, 9, 8, 8, 7, 0 };

                // last entry
                tLength( tLength.length()-1  )
                        = tLength.length() - tStart( tLength.length()-1 );

                Vector< real > tValues( 6, BELFEM_QUIET_NAN );

                // read values
                for( uint k=0; k<6; ++k )
                {
                    if( tStart( k ) < tLine.length() )
                    {
                        std::string tValue = tLine.substr( tStart( k ), tLength( k ));

                        // cleanup string
                        tValue = search_and_replace( tValue, "(", " " );
                        tValue = search_and_replace( tValue, ")", " " );
                        tValue = clean_string( tValue );

                        if ( tValue.length() > 0 )
                        {
                            tValues( k ) = std::stod( tValue );
                        }
                    }
                    else
                    {
                        break;
                    }
                }

                // copy data to data struct ( with right case )
                //aData->set_label( clean_string( tLine.substr( 0, 12 ) ) );

                // copy the name
                aData->set_name( clean_string( tLine.substr( 12, 60 ) ) );

                // the cas number
                aData->set_cas_number( clean_string( tLine.substr( 75, 16 ) ) );

                // the molar mass in kg / Mol
                aData->set_molar_mass(  tValues( 0 ) * 0.001 );

                // critical data
                aData->set_t_crit( tValues( 1 ) );
                aData->set_p_crit( tValues( 2 ) * 1e5 );
                aData->set_z_crit( tValues( 3 ) );

                // acentric factor
                aData->set_acentric_factor( tValues( 4 ) );

                // dipole moment
                aData->set_dipole_moment( tValues( 5 ) );
            }
        }

//------------------------------------------------------------------------------

        bool
        InputData::entry_exists( GasData * aData )
        {
            return mMap.key_exists( aData->label() );
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */