//
// Created by christian on 7/28/22.
//

#include "fn_create_beam_poly.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_Material_SAE301.hpp"
#include "nist_functions.hpp"
#include "fn_create_beam_poly.hpp"
#include "units.hpp"

namespace belfem
{
    namespace material
    {


//----------------------------------------------------------------------------

        SAE301::SAE301() :
                IsotropicMaterial( MaterialType::SAE301 )
        {
            mNumber = "1.4310";

            mTmax = 1173.15 ;
            this->create_expansion_poly() ;
            this->create_density_poly();
            this->create_mech_polys() ;

            this->create_specific_heat_polys() ;
            this->create_conductivity_polys();
            mHasThermal = true;
        }

//----------------------------------------------------------------------------

        void
        SAE301::create_expansion_poly()
        {
            // based on data from 10.1520/STP45025S
            mThermalExpansionPoly = {
                    -8.773394e-14,
                    4.620967e-09,
                    1.328919e-05,
                    0.000000e+00 } ;
        }

//----------------------------------------------------------------------------

        void
        SAE301::create_density_poly()
        {
            // based on data from 10.1520/STP45025S
            Vector< real > tT = {
                    116.483,
                    144.261,
                    199.817,
                    293.150,
                    366.483,
                    477.594,
                    588.706,
                    699.817,
                    810.928,
                    922.039,
                    1033.150,
                    1144.261,
                    1255.372 };

            Vector< real > tRho = {
                    7985.646,
                    7977.342,
                    7955.198,
                    7913.678,
                    7883.230,
                    7838.943,
                    7789.119,
                    7742.063,
                    7692.239,
                    7639.647,
                    7587.056,
                    7537.232,
                    7487.408
             };

            // first fit
            polyfit( tT, tRho, 2, mDensityPoly );

            // make consistent with thermal expansion
            IsotropicMaterial::create_density_poly( polyval( mDensityPoly, BELFEM_TREF ) );

        }

//----------------------------------------------------------------------------

        void
        SAE301::create_mech_polys()
        {
            // based on data from NASA TN D-2045, poisson guessed
            mYoungPoly = { -1.371312E+8, 2.331146E+11 };
            mShearPoly = { -5.402895E+07,  9.178941E+10 };

            mHasMechanical = true;
        }

//----------------------------------------------------------------------------

        void
        SAE301::create_specific_heat_polys()
        {
            // from 10.1520/STP45025S
            Vector< real > tT = {
                    116.483,
                    144.261,
                    199.817,
                    293.150,
                    366.483,
                    477.594,
                    588.706,
                    699.817,
                    810.928,
                    922.039,
                    1033.150,
                    1144.261,
                    1255.372 };

            Vector< real > tCp = {
                    334.944,
                    360.065,
                    401.933,
                    456.362,
                    489.856,
                    531.724,
                    556.845,
                    573.592,
                    586.152,
                    598.713,
                    619.647,
                    644.768,
                    678.262
                };

            polyfit( tT, tCp, 3, mSpecificHeatPoly5 );

            // cp at 200 K
            real tCpRef = polyval( mSpecificHeatPoly5, 200.0 );

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

            // compute connecting polynomial
            create_beam_poly(
                    mSwitchCT3,
                    nist::cp_janaf( mSpecificHeatPoly3, mSwitchCT3 ),
                    nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT3 ),
                    mSwitchCT4,
                    polyval(mSpecificHeatPoly5, mSwitchCT4 ),
                    dpolyval( mSpecificHeatPoly5, mSwitchCT4 ),
                    mSpecificHeatPoly4 );

        }

//----------------------------------------------------------------------------

        void
        SAE301::create_conductivity_polys()
        {
            // from 10.1520/STP45025S
            Vector< real > tT = {
                    //116.483,
                    //144.261,
                    //199.817,
                    //293.150,
                    366.483,
                    477.594,
                    588.706,
                    699.817,
                    810.928,
                    922.039,
                    1033.150,
                    1144.261,
                    1255.372 };

            Vector< real > tLambda = {
                    //11.077,
                    //11.769,
                    //12.981,
                    //14.884,
                    16.269,
                    18.000,
                    19.557,
                    21.115,
                    22.673,
                    24.057,
                    25.269,
                    26.480,
                    27.692 };

            polyfit( tT, tLambda, 1, mThermalConductivityPoly6 );

            // lambda at 200 K ( if all values are used and interpolation is quadratic
            real tLambdaRef = 12.946968 ;

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
                    = { -1.27655318e-5,
                        6.24092009e-3,
                        2.34353290e-1
                            };
            mThermalConductivityPoly4*= tLambdaRef ;


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
        SAE301::c( const real aT ) const
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
                return  nist::cp_janaf( mSpecificHeatPoly3, aT );
            }
            else if ( aT < mSwitchCT4 )
            {
                return  polyval(mSpecificHeatPoly4, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly5, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        SAE301::lambda( const real aT ) const
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