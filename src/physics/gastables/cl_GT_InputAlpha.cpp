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

#include <cmath>

#include "cl_GT_InputAlpha.hpp"
#include "fn_GT_parse.hpp"

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
            for ( const std::string & tLine : mBuffer )
            {
                // test if this line is empty or commented out
                if ( tLine.length() > 0 && tLine.at( 0 ) != '!' )
                {
                    // convert string to words
                    Cell< string > tWords = string_to_words( tLine );

                    // add cas number to map
                    if ( tWords.size() > 1 )
                    {
                        mMap[ tWords( 1 ) ] = tLineCount;
                    }
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

            if ( this->entry_exists( aData ) )
            {
                // data container
                Vector< real > tData( 12 );

                Cell< string > tWords = string_to_words( mBuffer( mMap( aData->cas() ) ) );

                BELFEM_ERROR( tWords.size() >= 14,
                        "Entry of %s in cubicalpha.inp has too few fields",
                        aData->label().c_str() );

                // read coefficients from file
                for( uint k=0; k<12; ++k )
                {
                    tData( k ) = parse_real( tWords( k+2 ), "coefficient in cubicalpha.inp" );
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
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
