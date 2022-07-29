//
// Created by christian on 7/28/22.
//

#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_Material_SAE314.hpp"
#include "nist_functions.hpp"
#include "fn_create_beam_poly.hpp"

namespace belfem
{
    namespace material
    {


//----------------------------------------------------------------------------

        SAE314::SAE314() :
                IsotropicMaterial( MaterialType::SAE314 )
        {
            mNumber = "1.4841";

            mTmax = 1273.15;
            this->create_expansion_poly() ;
            this->create_density_poly( 7800.0 );
            this->create_mech_polys() ;

            this->create_specific_heat_polys() ;
            this->create_conductivity_polys();
            mHasThermal = true;
        }

//----------------------------------------------------------------------------

        void
        SAE314::create_expansion_poly()
        {
            // based on data from Sirius Datasheet
            mThermalExpansionPoly = {
                    -2.859307e-13,
                    1.222925e-08,
                    5.956721e-06,
                    0.000000e+00 } ;
        }

//----------------------------------------------------------------------------

        void
        SAE314::create_mech_polys()
        {
            // based on datasheet
            mYoungPoly = { -5.617653e3, -7.310931e7, 2.216255e11} ;
            mShearPoly = { -6.750954e3, -2.241774e7, 8.216975e10 } ;
            mHasMechanical = true;
        }

//----------------------------------------------------------------------------

        void
        SAE314::create_specific_heat_polys()
        {
            // cp at 200 K
            real tCpRef = 416.4285148 ;

            // scaled data, 10 < T < 40
            mSpecificHeatPoly1 = {
                     3.394447978e-1,
                    -1.645011015e-1,
                     2.534585565e-2,
                    -5.307666665e-5,
                    -2.754876756e-5,
                     4.012000063e-6,
                    -3.558628323e-8 } ;

            mSpecificHeatPoly1 *= tCpRef ;

            // scaled data, 60 < T < 200
            mSpecificHeatPoly3 = {
                     9.683999191e2,
                    -2.069186178e1,
                    -1.022043858e0,
                     3.722615066e-2,
                    -2.793368253e-4,
                     9.994039027e-7,
                    -1.353454483e-9 } ;

            mSpecificHeatPoly3 *= tCpRef ;

            // compute deep cryo polynomial
            create_beam_poly(
                    0,
                    0,
                    0,
                    mSwitchCT0,
                    nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                    nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                    mSpecificHeatPoly0 );

            // compute connecting polynomial
            create_beam_poly(
                    mSwitchCT1,
                    nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                    nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                    mSwitchCT2,
                    nist::cp_janaf(mSpecificHeatPoly3, mSwitchCT2 ),
                    nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT2 ),
                    mSpecificHeatPoly2 );

            // linear extrapolation, so that cp @ 298.15 K = 500.0
            mSpecificHeatPoly4.set_size( 2 );

            mSpecificHeatPoly4( 0 )
                = nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT3 );
            mSpecificHeatPoly4( 1 )
                = nist::cp_janaf( mSpecificHeatPoly3, mSwitchCT3 )
                  - mSpecificHeatPoly4( 0 ) * mSwitchCT3 ;
        }

//----------------------------------------------------------------------------

        void
        SAE314::create_conductivity_polys()
        {
            // lambda at 200 K
            real tLambdaRef = 12.63269588 ;

            // T < 20 K
            mThermalConductivityPoly0
                = { -6.50800451e-6,
                     3.42504617e-4,
                     4.34442084e-3,
                     0.0 };
            mThermalConductivityPoly0 *= tLambdaRef ;

            // 40 < T < 80
            mThermalConductivityPoly2
                    = {  2.18252745e-7,
                        -8.93276489e-5,
                         1.50890960e-2,
                        -1.04825514e-1 };
            mThermalConductivityPoly2 *= tLambdaRef ;

            // 100 < T < 150
            mThermalConductivityPoly4
                    = {   8.98306626e-8,
                        - 4.68392476e-5,
                          1.05185318e-2,
                          5.66961918e-2 };
            mThermalConductivityPoly4*= tLambdaRef ;

            // T > 250
            mThermalConductivityPoly6 = { 2.650406e-2, 7.357599e0 };

            // connecting polynomials
            create_beam_poly(
                    mSwitchLT0,
                    polyval( mThermalConductivityPoly0, mSwitchLT0 ),
                    dpolyval( mThermalConductivityPoly0, mSwitchLT0 ),
                    mSwitchLT1,
                    polyval( mThermalConductivityPoly2, mSwitchLT1 ),
                    dpolyval( mThermalConductivityPoly2, mSwitchLT1 ),
                    mThermalConductivityPoly1 );

            create_beam_poly(
                    mSwitchLT2,
                    polyval( mThermalConductivityPoly2, mSwitchLT2 ),
                    dpolyval( mThermalConductivityPoly2, mSwitchLT2 ),
                    mSwitchLT3,
                    polyval( mThermalConductivityPoly4, mSwitchLT3 ),
                    dpolyval( mThermalConductivityPoly4, mSwitchLT3 ),
                    mThermalConductivityPoly3 );

            create_beam_poly(
                    mSwitchLT4,
                    polyval( mThermalConductivityPoly4, mSwitchLT4 ),
                    dpolyval( mThermalConductivityPoly4, mSwitchLT4 ),
                    mSwitchLT5,
                    polyval( mThermalConductivityPoly6, mSwitchLT5 ),
                    dpolyval( mThermalConductivityPoly6, mSwitchLT5 ),
                    mThermalConductivityPoly5 );

        }

//----------------------------------------------------------------------------

        real
        SAE314::c( const real aT ) const
        {
            if ( aT < mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if ( aT < mSwitchCT1 )
            {
                return nist::cp_janaf( mSpecificHeatPoly1, aT );
            }
            else if ( aT < mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if ( aT < mSwitchCT3 )
            {
                return nist::cp_janaf( mSpecificHeatPoly3, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly4, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        SAE314::lambda( const real aT ) const
        {
            if ( aT < mSwitchLT0 )
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
            else if ( aT < mSwitchLT1 )
            {
                return  polyval( mThermalConductivityPoly1, aT );
            }
            else if ( aT < mSwitchLT2 )
            {
                return  polyval( mThermalConductivityPoly2, aT );
            }
            else if ( aT < mSwitchLT3 )
            {
                return  polyval( mThermalConductivityPoly3, aT );
            }
            else if ( aT < mSwitchLT4 )
            {
                return  polyval( mThermalConductivityPoly4, aT );
            }
            else if ( aT < mSwitchLT5 )
            {
                return  polyval( mThermalConductivityPoly5, aT );
            }
            else
            {
                return  polyval( mThermalConductivityPoly6, aT );
            }
        }

//----------------------------------------------------------------------------
    }
}

//----------------------------------------------------------------------------