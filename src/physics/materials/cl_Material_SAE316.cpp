//
// Created by christian on 7/28/22.
//

#include "fn_create_beam_poly.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_Material_SAE316.hpp"
#include "nist_functions.hpp"
#include "fn_create_beam_poly.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"

namespace belfem
{
    namespace material
    {


//----------------------------------------------------------------------------

        SAE316::SAE316() :
                IsotropicMaterial( MaterialType::SAE316 )
        {
            mNumber = "1.4401";

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
        SAE316::create_expansion_poly()
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
        SAE316::create_mech_polys()
        {


            // create vandermonde matrix

            real tT1 = mSwitchM1 ;
            real tT2 = 100 ;
            real tT3 = 200 ;
            real tT4 = 300 ;
            Matrix< real > tA( 4, 4 );

            tA( 0,0 ) = 3.0 * tT1 * tT1 ;
            tA( 1,0 ) = tT2 * tT2 * tT2 ;
            tA( 2,0 ) = tT3 * tT3 * tT3 ;
            tA( 3,0 ) = tT4 * tT4 * tT4  ;

            tA( 0,1 ) = 2.0 * tT1 ;
            tA( 1,1 ) = tT2 * tT2 ;
            tA( 2,1 ) = tT3 * tT3 ;
            tA( 3,1 ) = tT4 * tT4 ;

            tA( 0,2 ) = 1.0 ;
            tA( 1,2 ) = tT2 ;
            tA( 2,2 ) = tT3 ;
            tA( 3,2 ) = tT4 ;

            tA( 0,3 ) = 0.0 ;
            tA( 1,3 ) = 1.0 ;
            tA( 2,3 ) = 1.0 ;
            tA( 3,3 ) = 1.0 ;

            Matrix< real > tB = tA ;
            Matrix< real > tC = tA ;

            Vector< int >  tWork(4 );

            // data based from NBSIR 80-1627
            mYoungPoly2 = { 0.0,  208.387751e9, 201.885361e9, 194.228713e9 };
            mShearPoly2 = { 0.0,   81.190236e9,  78.314787e9,  75.010562e9 };

            gesv( tA, mYoungPoly2, tWork );
            gesv( tB, mShearPoly2, tWork );

            // low temperatuere polynomials
            tT1 = 0.0 ;
            tT2 = 25.0 ;
            tT3 = mSwitchM1 ;

            tA( 0,0 ) = 3.0 * tT1 * tT1 ;
            tA( 1,0 ) = tT2 * tT2 * tT2 ;
            tA( 2,0 ) = tT3 * tT3 * tT3 ;
            tA( 3,0 ) = 3.0 * tT3 * tT3 ;

            tA( 0,1 ) = 2.0 * tT1 ;
            tA( 1,1 ) = tT2 * tT2 ;
            tA( 2,1 ) = tT3 * tT3 ;
            tA( 3,1 ) = 2.0 * tT3 ;

            tA( 0,2 ) = 1.0 ;
            tA( 1,2 ) = tT2 ;
            tA( 2,2 ) = tT3 ;
            tA( 3,2 ) = 1.0 ;

            tA( 0,3 ) = 0.0 ;
            tA( 1,3 ) = 1.0 ;
            tA( 2,3 ) = 1.0 ;
            tA( 3,3 ) = 0.0 ;

            // backup matrix
            tB = tA ;

            mYoungPoly1.set_size( 4, 0.0 );
            mYoungPoly1( 1 ) = 208.12312848e9 ;
            mYoungPoly1( 2 ) = polyval( mYoungPoly2, tT3 );
            gesv( tA, mYoungPoly1, tWork );


            mShearPoly1.set_size( 4, 0.0 );
            mShearPoly1( 1 ) =  81.2433079e9 ;
            mShearPoly1( 2 ) =  polyval( mShearPoly2, tT3 );

            gesv( tB, mShearPoly1, tWork );

            // linear extrapolation
            mYoungPoly3.set_size( 2 );
            mYoungPoly3( 0 ) = dpolyval( mYoungPoly2, mSwitchM2 );
            mYoungPoly3( 1 ) =  polyval( mYoungPoly2, mSwitchM2 )
                    - mYoungPoly3( 0 ) * mSwitchM2 ;

            // based on datasheet, but use only data for very high temperature
            Vector< real > tE = { -5.617653e3, -7.310931e7, 2.216255e11 } ;
            Vector< real > tG = { -6.750954e3, -2.241774e7, 8.216975e10 } ;

            tT1 = mSwitchM2 ;
            tT2 = 800.0 ;
            tT3 = 1200.0 ;

            tA( 0,0 ) = 3.0 * tT1 * tT1 ;
            tA( 1,0 ) = tT1 * tT1 * tT1 ;
            tA( 2,0 ) = tT2 * tT2 * tT2 ;
            tA( 3,0 ) = tT3 * tT3 * tT3 ;

            tA( 0,1 ) = 2.0 * tT1 ;
            tA( 1,1 ) = tT1 * tT1 ;
            tA( 2,1 ) = tT2 * tT2 ;
            tA( 3,1 ) = tT3 * tT2 ;

            tA( 0,2 ) = 1.0 ;
            tA( 1,2 ) = tT1 ;
            tA( 2,2 ) = tT2 ;
            tA( 3,2 ) = tT3 ;

            tA( 0,3 ) = 0.0 ;
            tA( 1,3 ) = 1.0 ;
            tA( 2,3 ) = 1.0 ;
            tA( 3,3 ) = 1.0 ;

            tB = tA ;

            mYoungPoly3.set_size( 4 );
            mYoungPoly3( 0 ) = dpolyval( mYoungPoly2, tT1 );
            mYoungPoly3( 1 ) = polyval( mYoungPoly2, tT1 );
            mYoungPoly3( 2 ) = polyval( tE, tT2 );
            mYoungPoly3( 3 ) = polyval( tE, tT3 );
            gesv( tA, mYoungPoly3, tWork );

            mShearPoly3.set_size( 4 );
            mShearPoly3( 0 ) = dpolyval( mShearPoly2, tT1 );
            mShearPoly3( 1 ) = polyval( mShearPoly2, tT1 );
            mShearPoly3( 2 ) = polyval( tG, tT2 );
            mShearPoly3( 3 ) = polyval( tG, tT3 );
            gesv( tB, mShearPoly3, tWork );

            mHasMechanical = true;

        }

//----------------------------------------------------------------------------

        void
        SAE316::create_specific_heat_polys()
        {
            // data from https://doi.org/10.1520/STP45025S
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
                    1255.372 } ;

            Vector< real > tCp = {
                    330.757,
                    355.878,
                    393.560,
                    452.175,
                    485.669,
                    527.537,
                    548.471,
                    565.218,
                    573.592,
                    586.152,
                    615.460,
                    648.955,
                    690.823 };

            // cp at 200 K, newer source
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

            // scaled data, 60 < T < 150
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


            polyfit( tT, tCp, 3, mSpecificHeatPoly5 );

            // compute turning point, d2cpDt2 == 0
            // let cp = a*T^3 + b*T^2 + c*T + d
            mSwitchCT4 = - mSpecificHeatPoly5( 1 ) /
                    ( 3.0 * mSpecificHeatPoly5( 0 ) );

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
        SAE316::create_conductivity_polys()
        {
            // data from Watson: Thermal Conductivity of a Sample of Type 316 Stainless Steel, 1963
            Vector< real > tT = {
                    363.15,
                    403.15,
                    443.15,
                    483.15,
                    523.15,
                    563.15,
                    603.15,
                    643.15,
                    683.15,
                    723.15,
                    763.15,
                    803.15,
                    843.15,
                    883.15,
                    923.15,
                    963.15,
                    1003.15,
                    1043.15,
                    1083.15,
                    1123.15 };

            Vector< real > tK = {
                    14.9,
                    15.5,
                    16.2,
                    16.8,
                    17.4,
                    18.1,
                    18.7,
                    19.3,
                    19.9,
                    20.5,
                    21.1,
                    21.8,
                    22.4,
                    23.0,
                    23.6,
                    24.3,
                    24.9,
                    25.6,
                    26.2,
                    26.9 };
            polyfit( tT, tK, 1, mThermalConductivityPoly6 );

            // lambda at 200 K
            real tLambdaRef = 12.63269588;

            // see G. Ventura and M. Perfetti, Thermal Properties of Solids at Room and Cryogenic Temperatures (Dordrecht, the
            //Netherlands: Springer Netherlands, 2014).

            // T < 20 K
            mThermalConductivityPoly0
                    = { -6.50800451e-6,
                        3.42504617e-4,
                        4.34442084e-3,
                        0.0 };
            mThermalConductivityPoly0 *= tLambdaRef;

            // 40 < T < 80
            mThermalConductivityPoly2
                    = { 2.18252745e-7,
                        -8.93276489e-5,
                        1.50890960e-2,
                        -1.04825514e-1 };
            mThermalConductivityPoly2 *= tLambdaRef;

            // 100 < T < 150
            mThermalConductivityPoly4
                    = { -1.27655318e-5,
                         6.24092009e-3,
                         2.34353290e-1
                        };

            mThermalConductivityPoly4 *= tLambdaRef;

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

            // this is a bit tricky because we want this to be smooth
            real tSwitchLT5 = 0.0;
            mSwitchLT5 = 300.0;

            uint tCount = 0 ;

            while( std::abs( tSwitchLT5 - mSwitchLT5 ) > 1e-12 )
            {
                create_beam_poly(
                        mSwitchLT4,
                        polyval( mThermalConductivityPoly4, mSwitchLT4 ),
                        dpolyval( mThermalConductivityPoly4, mSwitchLT4 ),
                        mSwitchLT5,
                        polyval( mThermalConductivityPoly6, mSwitchLT5 ),
                        dpolyval( mThermalConductivityPoly6, mSwitchLT5 ),
                        mThermalConductivityPoly5 );
                // shift
                tSwitchLT5 = mSwitchLT5;
                mSwitchLT5 *= 0.1 ;
                mSwitchLT5 += -0.9 * mThermalConductivityPoly5( 1 ) / ( 3.0 * mThermalConductivityPoly5( 0 ));
                if( std::abs( tSwitchLT5 - mSwitchLT5 ) < 1e-12 )
                {
                    break ;
                }
                if( tCount++ > 100 )
                {
                    BELFEM_ERROR( false, "too many iterations");
                }
            }
        }

//----------------------------------------------------------------------------

        real
        SAE316::c( const real aT ) const
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
            else if ( aT < mSwitchCT4 )
            {
                return polyval( mSpecificHeatPoly4, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly5, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        SAE316::lambda( const real aT ) const
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

        real
        SAE316::E( const real aT ) const
        {
            if( aT < mSwitchM1 )
            {
                return polyval( mYoungPoly1, aT );
            }
            else if ( aT < mSwitchM2 )
            {
                return polyval( mYoungPoly2, aT );
            }
            else
            {
                return polyval( mYoungPoly3, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        SAE316::G( const real aT ) const
        {
            if( aT < mSwitchM1 )
            {
                return polyval( mShearPoly1, aT );
            }
            else if ( aT < mSwitchM2 )
            {
                return polyval( mShearPoly2, aT );
            }
            else
            {
                return polyval( mShearPoly3, aT );
            }
        }

//----------------------------------------------------------------------------
    }
}

//----------------------------------------------------------------------------