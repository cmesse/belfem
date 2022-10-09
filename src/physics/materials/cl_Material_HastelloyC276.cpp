//
// Created by Christian Messe on 18.04.22.
//

#include "cl_Material_HastelloyC276.hpp"

#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        HastelloyC276::HastelloyC276() :
                IsotropicMaterial( MaterialType::HastelloyC276 )
        {
            // set maximum temperature
            mTmax = 1594;
            mNumber = "2.4819";

            this->create_resistivity_polys() ;
            this->create_specific_heat_polys();
            this->create_conductivity_polys();
            // this->create_mech_polys();
            this->create_expansion_poly();
            this->create_density_poly( 8890.0, 295.372 );

            mHasThermal = true;
            mHasMechanical = false ;
            mHasExpansion = true;
            mHasResistivity = false ;
        }

//----------------------------------------------------------------------------

        void
        HastelloyC276::create_resistivity_polys()
        {
            // from DOI: 10.1063/1.2899058
            Vector< real > tT1 = { 3.996, 4.098, 4.611, 5.020, 5.584, 6.045,
                                   6.660, 7.121, 7.633, 8.094, 8.504, 9.016,
                                   9.477, 10.195, 10.605, 11.219, 11.885,
                                   12.346, 13.012, 13.678, 14.344, 15.113,
                                   15.932, 16.650, 17.572, 18.391, 19.365,
                                   20.338, 21.260, 22.439, 23.668, 25.000,
                                   26.383, 27.971, 29.713, 31.609, 33.658,
                                   36.014, 38.576, 41.598, 44.980, 48.719 };

            Vector< real > tRho1 = { 1.2299028, 1.2296877, 1.2292958, 1.2290051,
                                     1.2287019, 1.2284998, 1.2282599, 1.2280832,
                                     1.2279698, 1.2278816, 1.2277808, 1.2276927,
                                     1.2276045, 1.2276053, 1.2275171, 1.2274671,
                                     1.2274804, 1.2276201, 1.2276208, 1.2276214,
                                     1.2276095, 1.2275976, 1.2277123, 1.2277131,
                                     1.2277899, 1.2279047, 1.2279816, 1.2281092,
                                     1.2282114, 1.2284278, 1.2284923, 1.2287089,
                                     1.2289381, 1.2291169, 1.2293086, 1.2296269,
                                     1.2299328, 1.2301124, 1.2304062, 1.2309282,
                                     1.231438, 1.2320367 };
            tRho1 *= 1E-6;
            Vector< real > tT3 = { 50.310, 53.360, 57.940, 64.050, 69.650, 77.290,
                                   84.420, 93.080, 102.240, 114.460, 125.660, 140.430,
                                   154.180, 168.940, 182.690, 197.960, 212.730, 226.480,
                                   240.730, 254.990, 269.760, 284.010, 299.290, 313.540,
                                   327.290, 341.040, 354.790, 369.550, 382.280, 393.990, 404.680 };

            Vector< real > tRho3 = { 1.231363, 1.232045, 1.232424, 1.232954, 1.233939,
                                     1.234848, 1.235681, 1.236742, 1.237954, 1.239014,
                                     1.24053, 1.24212, 1.243787, 1.245529, 1.24712,
                                     1.248711, 1.250074, 1.251665, 1.253407, 1.255757,
                                     1.257954, 1.259393, 1.260756, 1.261968, 1.263483,
                                     1.264694, 1.266209, 1.267572, 1.269012, 1.270148,
                                     1.271056 };
            tRho3 *= 1e-6;

            polyfit( tT1, tRho1, 5, mElectricResistivityPoly0 );
            polyfit( tT3, tRho3, 1, mElectricResistivityPoly2 );

            // connecting polynomial
            create_beam_poly(
                    mSwitchRhoT0, // T = 25 K
                    polyval( mElectricResistivityPoly0, mSwitchRhoT0 ), // rho @T=25K
                    dpolyval( mElectricResistivityPoly0, mSwitchRhoT0 ), // drhodT @T=25 K,
                    mSwitchRhoT1,// T = 60 K
                    polyval( mElectricResistivityPoly2, mSwitchRhoT1 ), // rho @T=60 K
                    dpolyval( mElectricResistivityPoly2, mSwitchRhoT1 ), // drhodT @T=60 K
                    mElectricResistivityPoly1 );

        }

//----------------------------------------------------------------------------

        void
        HastelloyC276::create_specific_heat_polys()
        {
            // from DOI: 10.1063/1.2899058

            Vector< real > tT1 = { 2.158, 2.447, 2.732, 3.074, 3.431, 3.830,
                                   3.891, 4.276, 4.849, 5.498, 6.186, 6.906,
                                   7.769, 8.673, 9.758, 10.892, 12.158, 13.785,
                                   15.388, 17.311, 19.320, 21.733 };

            Vector< real > tC1 = { 0.8360, 0.9273, 1.0094, 1.0989, 1.2076,
                                   1.2779, 1.3270, 1.4176, 1.4861, 1.5877,
                                   1.6961, 1.8291, 2.0099, 2.2087, 2.4501,
                                   2.7695, 3.1602, 3.7448, 4.4790, 5.4591,
                                   7.0402, 9.2527 };

            Vector< real > tT3 = { 24.640, 27.717, 30.932, 35.068, 39.138,
                                   43.680, 48.750, 54.839, 62.178, 69.949,
                                   78.080, 88.536, 98.827 };

            Vector< real > tC3 = { 12.3911, 16.5939, 22.4335, 30.6142,
                                   41.7780, 55.4214, 72.8389, 93.0578,
                                   116.6692, 140.8706, 168.4901, 195.9291,
                                   221.4763 };

            Vector< real > tT5 = { 112.061, 123.137, 139.638, 155.881,
                                   175.379, 200.457, 220.287, 245.912,
                                   283.303 };

            Vector< real > tC5 = { 250.3546, 275.0980, 302.3103, 316.9341,
                                   341.7811, 358.3143, 379.1899, 390.1094,
                                   405.1284 };

            polyfit( tT1, tC1, 3, mSpecificHeatPoly1 );
            polyfit( tT3, tC3, 3, mSpecificHeatPoly3 );
            polyfit( tT5, tC5, 3, mSpecificHeatPoly5 );

            // connecting polynomials
            create_beam_poly(
                    mSwitchCT1, // T = 20 K
                    polyval( mSpecificHeatPoly1, mSwitchCT1 ), // lambda @T=20K
                    dpolyval( mSpecificHeatPoly1, mSwitchCT1 ), // dlambdadT @T=20 K,
                    mSwitchCT2,// T = 40 K
                    polyval( mSpecificHeatPoly3, mSwitchCT2 ), // lambda @T=40 K
                    dpolyval( mSpecificHeatPoly3, mSwitchCT2 ), // dlambdadT @T=40 K
                    mSpecificHeatPoly2 );

            create_beam_poly(
                    mSwitchCT3, // T = 90 K
                    polyval( mSpecificHeatPoly3, mSwitchCT3 ), // lambda @T=90K
                    dpolyval( mSpecificHeatPoly3, mSwitchCT3 ), // dlambdadT @T=90 K,
                    mSwitchCT4,// T = 110 K
                    polyval( mSpecificHeatPoly5, mSwitchCT4 ), // lambda @T=110 K
                    dpolyval( mSpecificHeatPoly5, mSwitchCT4 ), // dlambdadT @T=110 K
                    mSpecificHeatPoly4 );

            // the last curve is just a linear extrapolation
            mSpecificHeatPoly6.set_size( 2 );
            mSpecificHeatPoly6( 0 ) = dpolyval( mSpecificHeatPoly5, mSwitchCT5 );
            mSpecificHeatPoly6( 1 ) = polyval( mSpecificHeatPoly5, mSwitchCT5 )
                                      - mSpecificHeatPoly6( 0 ) * mSwitchCT5;

        }

//----------------------------------------------------------------------------

        real
        HastelloyC276::rho_el( const real aJ, const real aT, const real aB, const real aAngle ) const
        {
            if ( aT < mSwitchRhoT0 )
            {
                return polyval( mElectricResistivityPoly0, aT );
            }
            else if ( aT < mSwitchRhoT1 )
            {
                return polyval( mElectricResistivityPoly1, aT );
            }
            else
            {
                return polyval( mElectricResistivityPoly2, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        HastelloyC276::c( const real aT ) const
        {
            if ( aT < mSwitchCT1 )
            {
                return polyval( mSpecificHeatPoly1, aT );
            }
            else if ( aT < mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if ( aT < mSwitchCT3 )
            {
                return polyval( mSpecificHeatPoly3, aT );
            }
            else if ( aT < mSwitchCT4 )
            {
                return polyval( mSpecificHeatPoly4, aT );
            }
            else if ( aT < mSwitchCT5 )
            {
                return polyval( mSpecificHeatPoly5, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly6, aT );
            }
        }

//----------------------------------------------------------------------------

        void
        HastelloyC276::create_conductivity_polys()
        {
            Vector< real > tT1 = { 14.725, 15.365, 16.389, 17.798, 18.822,
                                   19.462, 20.871, 21.639, 23.175, 24.2,
                                   25.864, 27.529, 29.577, 30.986, 32.522,
                                   36.108, 39.181, 41.485, 44.558, 48.656,
                                   53.009 };

            Vector< real > tLambda1 = { 3.4671, 3.5711, 3.7559, 3.9522, 4.0675,
                                        4.2641, 4.4488, 4.5758, 4.6679, 4.8296,
                                        5.1068, 5.2103, 5.5337, 5.6027, 5.7989,
                                        5.9828, 6.2595, 6.4671, 6.6049, 6.9276,
                                        7.0186 };

            Vector< real > tT3 = { 57.618, 63.252, 70.166, 76.825, 85.787,
                                   95.262, 105.506, 117.542, 131.114, 144.43,
                                   158.259, 171.831, 185.659, 199.232 };

            Vector< real > tLambda3 = { 7.1791, 7.3854, 7.8229, 7.728, 7.7248,
                                        8.1845, 8.3893, 8.5703, 8.8434, 9.1628,
                                        9.4821, 9.732, 10.1439, 10.5558 };

            polyfit( tT1, tLambda1, 2, mThermalConductivityPoly0 );
            polyfit( tT3, tLambda3, 1, mThermalConductivityPoly2 );

            create_beam_poly(
                    mSwitchLambdaT0, // T = 50 K
                    polyval( mThermalConductivityPoly0, mSwitchLambdaT0 ), // lambda @T=50K
                    dpolyval( mThermalConductivityPoly0, mSwitchLambdaT0 ), // dlambdadT @T=50 K,
                    mSwitchLambdaT1,// T = 100 K
                    polyval( mThermalConductivityPoly2, mSwitchLambdaT1 ), // lambda @T=100 K
                    dpolyval( mThermalConductivityPoly2, mSwitchLambdaT1 ), // dlambdadT @T=100 K
                    mThermalConductivityPoly1 );

        }

//----------------------------------------------------------------------------

        void
        HastelloyC276::create_expansion_poly()
        {
            // from DOI: 10.1063/1.2899058
            mThermalExpansionPoly = {  2.316147e-08,
                                       4.595386e-06,
                                       0.000000e+00 };
            mHasExpansion = true ;
        }


//----------------------------------------------------------------------------

        real
        HastelloyC276::lambda( const real aT ) const
        {
            if ( aT < mSwitchLambdaT0 )
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
            else if ( aT < mSwitchLambdaT1 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        HastelloyC276::mu( const real aT,   const real aTref ) const
        {
            return   polyval( mIntAlphaPoly, aT )
                   - polyval( mIntAlphaPoly, aTref ) - 1.0 ;
        }

//----------------------------------------------------------------------------
    }
}