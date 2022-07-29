//
// Created by Christian Messe on 27.07.20.
//

#include "cl_Material_Aluminum.hpp"
#include "nist_functions.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_linspace.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Aluminum::Aluminum() :
                IsotropicMaterial( MaterialType::Aluminum )
        {
            // set maximum temperature
            mTmax = 933.45;

            // sources:
            // Sutton,                 1953,  10.1103/PhysRev.91.816
            // Kamm and Alers,         1964 , 10.1063/1.1713309
            // Gerlich and Fischer,    1969,  10.1016/0022-3697(69)90377-1
            // Tallonand A. Wolfenden, 1979,  10.1016/0022-3697(79)90037-4
            // McLellan and Isikawa,   1987,  10.1016/0022-3697(87)90147-8
            mYoungPoly = { 8.85028e-4, -6.74178e0, -9.41290e3, -2.77173e7, 7.96165e10 };
            mShearPoly = { 4.70082e-3, -1.39404e1, 3.74788e3, -1.16598e7, 2.96729e10 };

            this->create_specific_heat_polys();
            this->create_conductivity_polys();

            // The cryogenic thermal conductivity varies with purity.


            // Touloukian from 100 K to melting point
            mThermalExpansionPoly = { -6.99557e-12, 2.08484e-8, 1.14374e-5, 0.0 };

            this->create_density_poly( 2700.0 );

            // assumption for oxidized surface
            mEpsilonPoly.set_size( 1, 0.4 );

            mHasThermal = true;
            mHasMechanical = true;
            mHasExpansion = true;

        }

//----------------------------------------------------------------------------

        real
        Aluminum::lambda( const real aT ) const
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

        void
        Aluminum::create_specific_heat_polys()
        {
            // cryogenic polynomial, taken from
            // https://trc.nist.gov/cryogenics/materials/6061%20Aluminum/6061_T6Aluminum_rev.htm
            // see also doi: 10.1007/0-306-47112-4_84

            Vector <real> tC = { 0.07918, 1.09570, -0.07277, 0.08084,
                                 0.02803, -0.09464, 0.04179, -0.00571 };

            // populate low cryogenic data

            // extrapolated from above data using Excel for T < 50 K
            mSpecificHeatPoly0 = { -8.387573E-05, -1.369301E-06, +1.444508E+00, 0.0 };

            // populate cryogenic data
            uint tN2 = std::floor( mSwitchCT2 - mSwitchCT1 ) + 1;
            Vector <real> tT2;
            linspace( mSwitchCT1, mSwitchCT2, tN2, tT2 );
            Vector <real> tC2( tN2 );
            for ( uint k = 0; k < tN2; ++k )
            {
                tC2( k ) = nist::proppoly( tC, tT2( k ));
            }

            polyfit( tT2, tC2, 3, mSpecificHeatPoly2 );

            // Touloukian and JANAF
            mSpecificHeatPoly4 = { 1.041606e-11, -3.165283e-8, 3.833785e-5, -2.297527e-2, 7.221121e0, 0.0 };

            create_beam_poly(
                    mSwitchCT0,
                    polyval( mSpecificHeatPoly0, mSwitchCT0 ),
                    dpolyval( mSpecificHeatPoly0, mSwitchCT0 ),
                    mSwitchCT1,
                    polyval( mSpecificHeatPoly2, mSwitchCT1 ),
                    dpolyval( mSpecificHeatPoly2, mSwitchCT1 ),
                    mSpecificHeatPoly1 );

            create_beam_poly(
                    mSwitchCT2,
                    polyval( mSpecificHeatPoly2, mSwitchCT2 ),
                    dpolyval( mSpecificHeatPoly2, mSwitchCT2 ),
                    mSwitchCT3,
                    polyval( mSpecificHeatPoly4, mSwitchCT3 ),
                    dpolyval( mSpecificHeatPoly4, mSwitchCT3 ),
                    mSpecificHeatPoly3 );
        }

//----------------------------------------------------------------------------

        void
        Aluminum::create_conductivity_polys()
        {
            // Reference data from Powel, Ho, Liley: Thermal Conductivity of Selected Materials, NSDRS, 1966
            Vector <real> tT0 = {
                    0.0, // << -- leave 0K in!
                    92.847,
                    93.603,
                    93.748,
                    96.117,
                    106.381,
                    106.417,
                    106.58,
                    115.695,
                    116.27,
                    117.22,
                    117.817,
                    122.875,
                    122.995,
                    125.368,
                    125.76,
                    136.838,
                    136.847,
                    137.143 };

            Vector <real> tL0 = {
                    0.0,
                    316.0892,
                    312.117,
                    312.4894,
                    306.0897,
                    278.2291,
                    280.1324,
                    281.4564,
                    266.0914,
                    267.2498,
                    266.2843,
                    261.2502,
                    258.1604,
                    257.1535,
                    258.3808,
                    257.1533,
                    250.615,
                    251.1115,
                    252.1597 };

            // Touloukian
            Vector <real> tT2 = {
                    297.73,
                    373.25,
                    473.43,
                    572.85,
                    673.02,
                    773.18,
                    872.59,
                    934.63 };

            Vector <real> tL2 = {
                    236.888,
                    239.762,
                    237.965,
                    232.792,
                    228.138,
                    221.926,
                    214.935,
                    210.793 };

            // lower polynomial
            polyfit( tT0, tL0, 4, mThermalConductivityPoly0 );

            // upper polynomial
            polyfit( tT2, tL2, 2, mThermalConductivityPoly2 );

            // connecting polynomial
            create_beam_poly(
                    mSwitchLambdaT0, // T = 130 K
                    polyval( mThermalConductivityPoly0, mSwitchLambdaT0 ), // lambda @T=130K
                    dpolyval( mThermalConductivityPoly0, mSwitchLambdaT0 ), // dlambdadT @T=130 K,
                    mSwitchLambdaT1, // T = 300 K
                    polyval( mThermalConductivityPoly2, mSwitchLambdaT1 ), // lambda @T=300 K
                    dpolyval( mThermalConductivityPoly2, mSwitchLambdaT1 ), // dlambdadT @T=300 K
                    mThermalConductivityPoly1 );

        }

//----------------------------------------------------------------------------

        real
        Aluminum::c( const real aT ) const
        {
            if ( aT > mSwitchCT3 )
            {
                return polyval( mSpecificHeatPoly4, aT );
            }
            else if ( aT > mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly3, aT );
            }
            else if ( aT > mSwitchCT1 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if ( aT > mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly1, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */