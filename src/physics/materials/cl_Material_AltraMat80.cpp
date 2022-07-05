//
// Created by Christian Messe on 02.07.20.
//

#include "cl_Material_AltraMat80.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        AltraMat80::AltraMat80() :
                IsotropicMaterial( MaterialType::Altramat80 )
        {

            // set maximum temperature
            mTmax = 1870.0;

            // polynomial from 200 K to 800 K ;
            mSpecificHeatPoly1 = {
                    2.69167849E-14,
                    -9.97699292E-11,
                    1.44731055E-07,
                    -1.00988631E-04,
                    3.13840963E-02,
                    -1.01556976E+00,
                    0.0 };

            // polynomial from 1200 K to 3000 K ;
            mSpecificHeatPoly3 = {
                    2.44097624E-08,
                    -1.67032085E-04,
                    4.42266451E-01,
                    9.58744394E+02
            };

            // create polynomial for cryogenic region
            create_beam_poly(
                    0.0, // T = 0 K
                    0.0, // cp @ T=0K
                    0.0, // dcpdT @0K
                    mSwitchCpT0, // T = 200 K
                    polyval( mSpecificHeatPoly1, mSwitchCpT0 ), // cp @T=200 K
                    dpolyval( mSpecificHeatPoly1, mSwitchCpT0 ), // dcpdT @T=200 K,
                    mSpecificHeatPoly0
            );

            // create connecting polynomial
            create_beam_poly(
                    mSwitchCpT1, // T = 800 K
                    polyval( mSpecificHeatPoly1, mSwitchCpT1 ), // cp @T=800 K
                    dpolyval( mSpecificHeatPoly1, mSwitchCpT1 ), // dcpdT @T=800 K,
                    mSwitchCpT2, // T = 1200 K
                    polyval( mSpecificHeatPoly3, mSwitchCpT2 ), // cp @T=200 K
                    dpolyval( mSpecificHeatPoly3, mSwitchCpT2 ), // dcpdT @T=200 K,
                    mSpecificHeatPoly2 );

            // set density
            mDensityPoly.set_size( 1, 80.0 );

            // thermal conductivity, taken form datasheet
            mThermalConductivityPoly = {
                    1.52057787E-13,
                    -1.50099379E-10,
                    9.88094187E-08,
                    9.44512445E-05,
                    0.0
            };

            mHasThermal = true;

        }

//----------------------------------------------------------------------------

        real
        AltraMat80::c( const real aT ) const
        {
            if ( aT < mSwitchCpT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if ( aT < mSwitchCpT1 )
            {
                return polyval( mSpecificHeatPoly1, aT );
            }
            else if ( aT < mSwitchCpT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly3, aT );
            }
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */