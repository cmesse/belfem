//
// Created by Christian Messe on 02.12.19.
//

#include "cl_AtmosphereModel_ISA1976.hpp"
#include "constants.hpp"

namespace belfem
{
    namespace atmosphere
    {
//------------------------------------------------------------------------------
        AtmosphereModel_ISA1976::AtmosphereModel_ISA1976() :
                mH(
                        {
                           0.0e0,  //  Troposphere
                           11.0e3, // Tropopause
                           20.0e3, // Stratosphere
                           32.0e3, // Stratosphere
                           47.0e3, // Stratopause
                           51.0e3, // Mesosphere
                           71.0e3, // Mesosphere
                           84.85e3,
                           90.69e3,
                           constant::R_earth * 100.0e3 /
                           ( constant::R_earth + 100.0e3 )
                   } ),
                mL( {
                            -6.5e-3, // Troposphere
                             0.0,    // Tropopause
                             1.0e-3, // Stratosphere
                             2.8e-3, // Stratosphere
                             0.0e00, // Stratopause
                            -2.8e-3, // Mesosphere
                            -2.0e-3,
                             0.0,
                             1.03e-3
                    } ),
              mMaxAltitude( 1000.0e3 )
        {
            // compute temperature and pressure at 100 km
            this->compute_T_and_p( 100.0e3, mT100, mP100 );

            mTheta100 = mT100 / mT0 + 2.85;
        }

//------------------------------------------------------------------------------

        real
        AtmosphereModel_ISA1976::max_altitude() const
        {
            return mMaxAltitude;
        }


//------------------------------------------------------------------------------

        void
        AtmosphereModel_ISA1976::compute_T_and_p(
                const real & aAltitude,
                      real & aTemperature,
                      real & aPressure ) const
        {
            if( aAltitude <= 100e3 )
            {
                // geopotential height
                real tH = constant::R_earth * aAltitude /
                          ( constant::R_earth + aAltitude );

                // set initial values
                aTemperature = mT0;
                aPressure = mP0;

                uint j = mH.length() - 1;

                real tT0;

                for ( uint k = 1; k < mH.length(); ++k )
                {

                    if ( mH( k ) >= tH )
                    {
                        j = k - 1;
                        break;
                    }
                    else
                    {
                        if ( mL( k - 1 ) != 0 )
                        {
                            tT0 = aTemperature;
                            aTemperature += mL( k - 1 ) * ( mH( k ) - mH( k - 1 ));
                            aPressure *= std::pow( aTemperature / tT0, -constant::g0 / ( mR * mL( k - 1 )));
                        }
                        else
                        {
                            aPressure *= std::exp( constant::g0 * ( mH( k - 1 ) - mH( k )) /
                                                   ( mR * aTemperature ));
                        }
                    }
                }

                if ( mL( j ) != 0 )
                {
                    tT0 = aTemperature;
                    aTemperature += mL( j ) * ( tH - mH( j ));

                    aPressure *= std::pow( aTemperature / tT0, -constant::g0 / ( mR * mL( j )));
                }
                else
                {
                    aPressure *= std::exp( constant::g0 * ( mH( j ) - tH ) / ( mR * aTemperature ));
                }
            }
            else
            {
                // eq. 2.13
                aTemperature = ( mTheta100 - 2.85 * std::exp( ( 100e3-aAltitude ) / 70e3 ) ) * mT0 ;

                // help magnitude
                real tX = ( aAltitude - 100e3 ) / 100e3;

                // interpolated from table
                aPressure = mP100 / std::exp( std::sqrt( ( ( 0.9373 - 0.1933*tX ) * tX + 33.0483 ) * tX ) );
            }
        }


//------------------------------------------------------------------------------
    }
}