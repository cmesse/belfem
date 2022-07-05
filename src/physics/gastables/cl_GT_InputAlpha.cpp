//
// Created by Christian Messe on 26.08.19.
//

#include "cl_GT_InputAlpha.hpp"

namespace belfem
{
    namespace gastables
    {
//----------------------------------------------------------------------------

        InputAlpha::InputAlpha( const string & aPath )
                : Ascii( aPath, FileMode::OPEN_RDONLY )
        {
            // finit counter
            uint tLineCount = 0;

            // loop over all lines
            for ( std::string tLine : mBuffer )
            {
                // test if this line is commented out
                if ( tLine.at( 0 ) != '!' )
                {
                    // convert string to words
                    Cell< string > tWords = string_to_words( tLine );

                    // add cas number to map
                    mMap[ tWords( 1 ) ] = tLineCount;
                }
                // increment line counter
                tLineCount++;
            }
        }

//----------------------------------------------------------------------------

        bool
        InputAlpha::entry_exists( GasData * aData )
        {
            return mMap.key_exists( aData->cas() );
        }

//----------------------------------------------------------------------------

        void
        InputAlpha::read_data( GasData * aData  )
        {

            Vector< real > & tSRK = aData->srk();
            Vector< real > & tPR  = aData->pr();

            tSRK.set_size( 3, 0.0 );
            tPR.set_size( 3, 0.0 );

            const real tOmega = aData->acentric();

            if ( this->entry_exists( aData ) )
            {
                // data container
                Vector< real > tData( 12 );

                Cell< string > tWords = string_to_words( mBuffer( mMap( aData->cas() ) ) );

                // read coefficients from file
                for( uint k=0; k<12; ++k )
                {
                    tData( k ) = std::stod( tWords( k+2 ) );
                }

                // overwrite critical temperature and pressure to be consistent with paper
                aData->set_t_crit( tData( 0 ) );
                aData->set_p_crit( tData( 1 )*1e5 );

                // get Z_crit
                real tZcrit = aData->Z_crit();

                // overwrite Z_crit in order to get rho_crit consistent with T_crit and p_crit
                aData->set_z_crit( tZcrit );

                // data for Soave-Redlich-Kwong
                tSRK( 0 ) =  2.0*tData( 4 );
                tSRK( 1 ) = -std::pow( tData( 5 ), 2 );
                tSRK( 2 ) =  2.0/3.0*std::pow( tData( 6 ), 3 );

                // data for Peng-Robinson
                tPR( 0 ) =  2.0*tData( 9 );
                tPR( 1 ) =  -std::pow( tData( 10 ), 2 );
                tPR( 2 ) =  2.0/3.0*std::pow( tData( 11 ), 3 );

                aData->set_cubic_flag();
            }
            else if (  ! std::isnan( tOmega ) ) // test if acentric factor exists
            {
                // fixme: these data are never used. This is OK

                // see 10.1023/B:IJOT.0000022331.46865.2f
                tSRK( 0 ) = tOmega * ( 0.1441 * tOmega + 1.3838 ) + 0.387;
                tSRK( 1 ) = tOmega * ( 0.6939 - 2.5214 * tOmega ) + 0.0325;
                tSRK( 2 ) = 0.6225 * tOmega + 0.2236;

                // data for Peng-Robinson
                tPR( 0 ) = tOmega * ( 1.3569  * tOmega + 0.9957 ) + 0.4077;
                tPR( 1 ) = tOmega * ( 3.5590 - 11.2986 * tOmega ) - 0.1146;
                tPR( 2 ) = tOmega * ( 11.7802 * tOmega - 3.8901 ) + 0.5033;
            }
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
