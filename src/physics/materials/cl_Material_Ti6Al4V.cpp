//
// Created by Christian Messe on 02.08.20.
//

#include "cl_Material_Ti6Al4V.hpp"
#include "nist_functions.hpp"
#include "cl_Vector.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_linspace.hpp"
#include "cl_Material.hpp"
#include "fn_gesv.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Ti6Al4V::Ti6Al4V()  :
                IsotropicMaterial( MaterialType::TI6AL4V )
        {
            // set maximum temperature
            mTmax = 1941.0;
            mNumber = "3.7165";

            this->create_mech_polys();
            this->create_specific_heat_polys();
            this->create_conductivity_polys();

            // not very accurate.
            // cryogenic data fron NIST
            // hot temperature data based on 10.1007/pl00021868
            mThermalExpansionPoly = { -5.67231914556140E-24, 4.39449998265093E-20,
                                      -1.33020127385256E-16, 1.96525261089848E-13,
                                      -1.44541515205732E-10, 5.59126298885355E-08,
                                      -2.03662130605073E-06, 0.00000000000000E+00 };


            this->create_density_poly();

            mHasMechanical = true;
            mHasThermal = true;
            mHasExpansion = true;


        }

//----------------------------------------------------------------------------

        void
        Ti6Al4V::create_density_poly()
        {
            /*
             * McCormick, A; Brooks, R F: MTS Programme on processability:
             * Thermophysical property data for commercial alloys measured
             * in PMPl9 2 and 3 (Apr 93-Mar 96), NPL (1996) Chapter 1.
             */
            Vector <real> tT = { 295.75, 371.92, 473.06, 574.21, 672.86, 774.01,
                                 872.65, 973.80, 1072.45, 1173.60, 1267.25, 1370.89,
                                 1470.79, 1570.68, 1670.58, 1771.74, 1871.62 };

            Vector <real> tR = { 4419.90, 4407.35, 4394.82, 4382.29, 4367.47,
                                 4351.50, 4337.83, 4325.29, 4310.47, 4294.51,
                                 4283.12, 4267.15, 4252.33, 4240.94, 4226.12,
                                 4205.58, 4198.77 };


            polyfit( tT, tR, 1, mDensityPoly );
            IsotropicMaterial::create_density_poly(
                    polyval( mDensityPoly, BELFEM_TREF ) );
        }

//----------------------------------------------------------------------------

        void
        Ti6Al4V::create_specific_heat_polys()
        {
            // these data have been created from mass averaging of
            // JANAF data. We use them to determine the low cryogenic polynomial
            Vector <real> tT0 = {
                    100.0,
                    200.0,
                    298.2,
                    300.0,
                    400.0,
                    500.0 };

            Vector <real> tC0 = {
                    308.7124774,
                    485.1740871,
                    547.9106219,
                    548.7347306,
                    582.9982398,
                    604.993376 };

            /**
             * Bros, H; Michel, M L; Castanet, R: J. Thermal Analysis 41 (1994) 7/24.
             * Richardson, M J; Hayes, D L; Day, A P; Mills, K C: as in ref 2, Chapter 3.
             */

            // data for alpha-configuration
            Vector <real> tT3 = { 292.54, 292.81, 295.11, 369.28, 375.51, 425.04,
                                  473.13, 474.42, 567.13, 572.01, 580.75, 672.24,
                                  674.63, 725.39, 773.63, 873.67, 873.79, 971.47,
                                  1076.59, 1171.80, 1268.27 };

            Vector <real> tC3 = { 560.09, 522.04, 546.07, 565.62, 560.11, 561.63,
                                  590.68, 583.17, 607.23, 617.24, 606.73, 628.28,
                                  640.80, 643.32, 650.84, 688.42, 672.40, 693.45,
                                  713.51, 734.06, 752.61 };

            // data for beta-configuration
            Vector <real> tT5 = { 1266.58, 1370.47, 1471.89, 1570.84, 1669.78, 1771.21, 1871.39 };
            Vector <real> tC5 = { 640.45, 659.50, 677.55, 694.60, 713.15, 730.70, 749.26 };


            // create temporary polynomial
            Vector <real> tPoly;
            polyfit( tT0, tC0, 4, tPoly );

            // alpha-polynomial
            polyfit( tT3, tC3, 2, mSpecificHeatPoly3 );

            // beta polynomial
            polyfit( tT5, tC5, 1, mSpecificHeatPoly5 );

            real tCa = polyval( tPoly, mSwitchCT0 );
            real dCa = dpolyval( tPoly, mSwitchCT0 );

            real tCb = polyval( tPoly, mSwitchCT1 );
            real dCb = dpolyval( tPoly, mSwitchCT1 );

            // low cryogenic data
            create_beam_poly(
                    0.0,
                    0.0,
                    0.0,
                    mSwitchCT0,
                    tCa,
                    dCa,
                    mSpecificHeatPoly0 );

            // high cryogenic data
            create_beam_poly(
                    mSwitchCT0,
                    tCa,
                    dCa,
                    mSwitchCT1,
                    tCb,
                    dCb,
                    mSpecificHeatPoly1 );

            // connecting polynomial 1
            create_beam_poly(
                    mSwitchCT1,
                    tCb,
                    dCb,
                    mSwitchCT2,
                    polyval( mSpecificHeatPoly3, mSwitchCT2 ),
                    dpolyval( mSpecificHeatPoly3, mSwitchCT2 ),
                    mSpecificHeatPoly2 );

            // connecting polynomial 2
            create_beam_poly(
                    mSwitchCT3,
                    polyval( mSpecificHeatPoly3, mSwitchCT3 ),
                    dpolyval( mSpecificHeatPoly3, mSwitchCT3 ),
                    mSwitchCT4,
                    polyval( mSpecificHeatPoly5, mSwitchCT4 ),
                    dpolyval( mSpecificHeatPoly5, mSwitchCT4 ),
                    mSpecificHeatPoly4 );

        }

//----------------------------------------------------------------------------

        void
        Ti6Al4V::create_conductivity_polys()
        {
            // data from NIST
            // https://trc.nist.gov/cryogenics/materials/Ti6Al4V/Ti6Al4V_rev.htm
            // see also 10.1007/0-306-47112-4_84

            mThermalConductivityPoly0 = { -5107.8774, 19240.422, -30789.064, 27134.756,
                                          -14226.379, 4438.2154, -763.07767, 55.796592 };

            uint tN2 = 176;
            Vector <real> tT2;
            linspace( mSwitchLambdaT1 - 25.0, mSwitchLambdaT2 + 25.0, tN2, tT2 );
            Vector <real> tL2( tN2 );

            // populate data
            for ( uint k = 0; k < tN2; ++k )
            {
                tL2( k ) = nist::proppoly( mThermalConductivityPoly0, tT2( k ));
            }

            // create cryo poly
            polyfit( tT2, tL2, 6, mThermalConductivityPoly2 );

            // create connecting poly
            create_beam_poly(
                    mSwitchLambdaT0,
                    nist::proppoly( mThermalConductivityPoly0, mSwitchLambdaT0 ),
                    nist::dproppoly( mThermalConductivityPoly0, mSwitchLambdaT0 ),
                    mSwitchLambdaT1,
                    polyval( mThermalConductivityPoly2, mSwitchLambdaT1 ),
                    dpolyval( mThermalConductivityPoly2, mSwitchLambdaT1 ),
                    mThermalConductivityPoly1 );

            // high temperature data
            // Polev, V F; Zinovyev, V E; Korshunov, IG: High Temperatures, 23 (1985) 704/706.
            Vector <real> tT4 = { 297.09, 374.95, 472.93, 572.18, 672.66, 773.16, 872.41, 974.16, 1072.18, 1171.46,
                                  1265.73 };
            Vector <real> tL4 = { 7.0508, 7.4765, 8.7991, 10.2073, 11.4017, 12.6389, 14.218, 15.5406, 17.8034, 20.109,
                                  22.7138 };
            Vector <real> tT6 = { 1265.59, 1369.86, 1470.38, 1570.85, 1670.08, 1770.57, 1869.8 };
            Vector <real> tL6 = { 19.295, 20.9167, 22.7949, 23.6475, 24.5855, 25.7372, 26.9317 };

            // create alpha poly
            polyfit( tT4, tL4, 2, mThermalConductivityPoly4 );

            // create beta poly
            polyfit( tT6, tL6, 2, mThermalConductivityPoly6 );

            // connecting polys
            create_beam_poly(
                    mSwitchLambdaT2,
                    polyval( mThermalConductivityPoly2, mSwitchLambdaT2 ),
                    dpolyval( mThermalConductivityPoly2, mSwitchLambdaT2 ),
                    mSwitchLambdaT3,
                    polyval( mThermalConductivityPoly4, mSwitchLambdaT3 ),
                    dpolyval( mThermalConductivityPoly4, mSwitchLambdaT3 ),
                    mThermalConductivityPoly3 );

            create_beam_poly(
                    mSwitchLambdaT4,
                    polyval( mThermalConductivityPoly4, mSwitchLambdaT4 ),
                    dpolyval( mThermalConductivityPoly4, mSwitchLambdaT4 ),
                    mSwitchLambdaT5,
                    polyval( mThermalConductivityPoly6, mSwitchLambdaT5 ),
                    dpolyval( mThermalConductivityPoly6, mSwitchLambdaT5 ),
                    mThermalConductivityPoly5 );
        }

//----------------------------------------------------------------------------

        void
        Ti6Al4V::create_mech_polys()
        {
            // data for alpha configuration
            // based on doi 10.1007/bf00420541

            Vector <real> tTE0 = { 296.3, 308.85, 323.47, 348.19, 354.29, 365.02,
                                   381.12, 418.43, 437.42, 448.1, 462.8, 474.66,
                                   484.34, 493.64, 501.17, 511.21, 518.37, 534.48,
                                   549.19, 557.79, 568.18, 577.86, 585.36, 601.13,
                                   614.06, 632.68, 643.78, 652.37, 663.12, 671.72,
                                   686.03, 689.65, 706.47, 720.49, 730.17, 745.55,
                                   757.37, 767.76, 781.37, 804, 821.14, 834., 836.21,
                                   853.38, 863.81, 887.07, 896.73, 906.04, 916.07,
                                   930.04, 982.3, 986.66, 995.27, 1012.8, 1019.98,
                                   1025.01, 1037.21, 1057.27, 1070.17, 1078.75,
                                   1087.35, 1096.68, 1107.8 };

            Vector <real> tE0 = { 103.342, 102.447, 103.86, 102.621, 101.922, 101.929,
                                  101.838, 98.7, 97.608, 99.169, 97.974, 96.426, 95.73,
                                  95.735, 95.238, 94.442, 94.446, 93.803, 92.508, 92.162,
                                  91.767, 91.171, 91.225, 90.532, 88.985, 88.344, 88.,
                                  87.854, 87.158, 86.762, 86.62, 85.368, 85.027, 83.18,
                                  82.583, 82.542, 81.948, 81.452, 80.858, 78.213,
                                  78.976, 79.335, 77.48, 77.44, 75.741, 75.404, 75.409,
                                  75.064, 74.719, 74.175, 73.003, 71.049, 70.201, 70.212,
                                  69.313, 68.814, 67.568, 66.777, 66.183, 66.088, 65.691,
                                  64.794, 63.998 };
            tE0 *= 1.e9;

            Vector <real> tTG0 = { 298.96, 309.73, 321.52, 327.57, 347.64, 356.58, 369.45,
                                   380.54, 412.06, 420.65, 433.87, 442.48, 452.11, 459.64,
                                   475.41, 484.72, 497.59, 507.59, 516.57, 525.15, 538.39,
                                   550.56, 561.3, 571.67, 584.92, 597.45, 605.68, 616.79,
                                   637.53, 652.22, 664.38, 676.18, 687.97, 704.81, 723.08,
                                   735.6, 746.34, 760.65, 771.74, 780.71, 806.12, 823.29,
                                   836.52, 849.79, 863.38, 887.36, 898.45, 909.53, 927.42,
                                   975.72, 983.97, 997.58, 1011.54, 1025.86, 1037.31, 1059.51,
                                   1074.18, 1084.91 };

            Vector <real> tG0 = { 39.18, 37.982, 38.391, 39.247, 38.106, 38.161, 38.119,
                                  38.075, 36.89, 36.745, 37.104, 36.206, 37.015, 36.568,
                                  35.674, 35.329, 35.387, 35.694, 34.696, 34.751, 34.659,
                                  34.315, 34.121, 33.976, 33.533, 33.289, 33.094, 32.649,
                                  32.561, 32.068, 31.975, 31.932, 32.391, 31.448, 30.656,
                                  30.664, 30.219, 30.127, 30.284, 29.286, 28.75, 28.71,
                                  28.868, 27.823, 27.981, 27.394, 27.451, 27.457, 27.468,
                                  26.895, 26.298, 25.754, 25.211, 24.818, 24.775, 23.986,
                                  23.944, 24.001 };

            tG0 *= 1E9;

            /* The following data are not usable, and only saved for the sake of completeness
            Vector< real > tTE1 = { 1122.84, 1136.43, 1147.16, 1158.6, 1168.63, 1174.75, 1186.25,
                                    1199.6, 1218.67 } ;
            Vector< real > tE1 = { 63.455, 63.463, 63.52, 63.527, 63.031, 62.082, 60.534, 57.231,
                                   53.931 } ;
            Vector< real > tTG1 = { 1101.02, 1111.04, 1118.56, 1125.37, 1138.24, 1148.61, 1160.77,
                                    1171.51, 1181.17, 1189.42, 1198.75, 1217.4 };
            Vector< real > tG1 = {23.459, 23.314, 23.068, 22.621, 22.729, 22.835, 22.843, 22.699,
                                  22.554, 21.756, 20.959, 19.566 }; */

            Vector <real> tTE2 = { 1233.53, 1247.51, 1251.15, 1263.04, 1275.91,
                                   1292.33, 1306.61, 1312.63, 1323.66 };
            Vector <real> tE2 = { 48.472, 37.744, 35.79, 33.489, 33.547, 34.46, 35.221,
                                  36.78, 38.342 };
            tE2 *= 1e9;
            Vector <real> tTG2 = { 1231.42, 1246.22, 1254.48, 1262.37, 1278.11, 1293.11,
                                   1310.28, 1319.54, 1324.53, 1324.53 };
            Vector <real> tG2 = { 17.417, 13.563, 12.615, 12.017,
                                  12.077, 12.538, 12.698, 13.707, 14.112, 14.112 };
            tG2 *= 1e9;

            // first, we create the alpha polynomial
            polyfit( tTE0, tE0, 1, mYoungPoly0 );
            polyfit( tTG0, tG0, 1, mShearPoly0 );


            // create temporary polynomials
            polyfit( tTE2, tE2, 3, mYoungPoly2 );
            polyfit( tTG2, tG2, 3, mShearPoly2 );

            // reference point, where the beta begins
            real tTa = 1283.0;
            real tEa = polyval( mYoungPoly2, tTa );
            real tGa = polyval( mShearPoly2, tTa );

            // poisson number at this point
            real tNua = tEa / ( 2.0 * tGa ) - 1.0;

            // assumption for poisson derivative based on polynomial 0

            real tdNua = ( polyval( mShearPoly0, tTa )
                           * dpolyval( mYoungPoly0, tTa )
                           - polyval( mYoungPoly0, tTa )
                             * dpolyval( mShearPoly0, tTa )) /
                         ( 2.0 * std::pow( polyval( mShearPoly0, tTa ), 2 ));

            // assumption for polynomial
            real tTb = mTmax + 200.0;
            real tNub = 0.5;

            // create young poly
            mYoungPoly2.set_size( 2 );
            mYoungPoly2( 0 ) = mYoungPoly0( 0 );
            mYoungPoly2( 1 ) = tEa - mYoungPoly2( 0 ) * tTa;

            // temporary polynomial for poisson number
            Matrix <real> tV( 3, 3 );
            tV( 0, 0 ) = tTa * tTa;
            tV( 1, 0 ) = 2.0 * tTa;
            tV( 2, 0 ) = tTb * tTb;
            tV( 0, 1 ) = tTa;
            tV( 1, 1 ) = 1.0;
            tV( 2, 1 ) = tTb;
            tV( 0, 2 ) = 1.0;
            tV( 1, 2 ) = 0.0;
            tV( 2, 2 ) = 1.0;

            Vector < real > tNuPoly( 3 );
            Vector < int > tPivot( 3 );
            tNuPoly( 0 ) = tNua;
            tNuPoly( 1 ) = tdNua;
            tNuPoly( 2 ) = tNub;
            gesv( tV, tNuPoly, tPivot );            

            // create samples for shrear modulus
            uint tN = 100;

            Vector <real> tT;
            linspace( tTa, tTb, tN, tT );
            Vector <real> tG( tN );

            for ( uint k = 0; k < tN; ++k )
            {
                tG( k ) = polyval( mYoungPoly2, tT( k )) /
                          ( 2.0 * ( 1.0 + polyval( tNuPoly, tT( k ))));
            }

            // create shear polynomial for beta phase
            polyfit( tT, tG, 3, mShearPoly2 );

            // beta configuration educated guess
            // assuming same inclination for young and
            // poisson at 0.5 for 2000 K


            create_beam_poly(
                    mSwitchMechT0,
                    polyval( mYoungPoly0, mSwitchMechT0 ),
                    dpolyval( mYoungPoly0, mSwitchMechT0 ),
                    mSwitchMechT1,
                    polyval( mYoungPoly2, mSwitchMechT1 ),
                    dpolyval( mYoungPoly2, mSwitchMechT1 ),
                    mYoungPoly1 );

            create_beam_poly(
                    mSwitchMechT0,
                    polyval( mShearPoly0, mSwitchMechT0 ),
                    dpolyval( mShearPoly0, mSwitchMechT0 ),
                    mSwitchMechT1,
                    polyval( mShearPoly2, mSwitchMechT1 ),
                    dpolyval( mShearPoly2, mSwitchMechT1 ),
                    mShearPoly1 );

        }

//----------------------------------------------------------------------------

        real
        Ti6Al4V::E( const real aT ) const
        {
            if ( aT < mSwitchMechT0 )
            {
                return polyval( mYoungPoly0, aT );
            }
            else if ( aT < mSwitchMechT1 )
            {
                return polyval( mYoungPoly1, aT );
            }
            else
            {
                return polyval( mYoungPoly2, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Ti6Al4V::G( const real aT ) const
        {
            if ( aT < mSwitchMechT0 )
            {
                return polyval( mShearPoly0, aT );
            }
            else if ( aT < mSwitchMechT1 )
            {
                return polyval( mShearPoly1, aT );
            }
            else
            {
                return polyval( mShearPoly2, aT );
            }
        }

//----------------------------------------------------------------------------
        real
        Ti6Al4V::c( const real aT ) const
        {
            if ( aT < mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if ( aT < mSwitchCT1 )
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
            else
            {
                return polyval( mSpecificHeatPoly5, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Ti6Al4V::lambda( const real aT ) const
        {
            if ( aT < mSwitchLambdaT0 )
            {
                return nist::proppoly( mThermalConductivityPoly0, aT );
            }
            else if ( aT < mSwitchLambdaT1 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else if ( aT < mSwitchLambdaT2 )
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
            else if ( aT < mSwitchLambdaT3 )
            {
                return polyval( mThermalConductivityPoly3, aT );
            }
            else if ( aT < mSwitchLambdaT4 )
            {
                return polyval( mThermalConductivityPoly4, aT );
            }
            else if ( aT < mSwitchLambdaT5 )
            {
                return polyval( mThermalConductivityPoly5, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly6, aT );
            }
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
