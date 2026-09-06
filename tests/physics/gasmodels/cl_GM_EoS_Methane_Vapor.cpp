//
// Created by Christian Messe on 18.08.20.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Gas.hpp"

#define private public
#define protected public
#include "cl_GM_EoS_Methane.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

TEST( GASMODELS, Methane_Vapor )
{
    // crate reference mixture
    Gas tRef( "CH4" );

    // create EoS
    EoS_Methane tCH4( tRef );

//----------------------------------------------------------------------------
//  Vapor Pressure
//----------------------------------------------------------------------------

    /*
     * Dortmund Databank. Based on
     *
     * 1	Kleinrahm R.; Duschek W.; Wagner W.:
     *      (Pressure, density, temperature) measurements in the critical
     *      region of methane. J.Chem.Thermodyn. 18 (1986) 1103-1114
     *
     * 2	Hunter M.A.:
     *      The Molecular Aggregation of Liquefied Gases.
     *      J.Phys.Chem. 10 (1906) 330-360
     *
     * 3	Van Itterbeek A.; Staes K.; Verbeke O.; Theeuwes F.:
     *      Vapour pressure of saturated liquid methane.
     *      Physica (Amsterdam) 30 (1964) 1896-1900
     *
     * 4	Stock A.; Henning F.; Kuß E.:
     *      Dampfdrucktafeln fuer Temperaturbestimmungen zwischen +25°C und -185°C.
     *      Ber.Dtsch.Chem.Ges./B:Abhandl. 54 (1921) 1119-1129
     *
     * 5	Karwat E.:
     *      Der Dampfdruck des festen Chlorwasserstoffs, Methans und Ammoniaks.
     *      Z.Phys.Chem.(Leipzig) 112 (1924) 486-490
     *
     * 6	Eucken A.; Berger W.: Das I-T-Diagramm des Methans.
     *      Z.Ges.Kälteindustrie 41 (1934) 145-152
     */
    Vector< real > tTvap = { 76.89, 79.75, 81.74, 83.82, 85.42, 87.25, 92.15,
                             93.15, 94.15, 95.15, 96.15, 96.30, 97.15, 97.20,
                             97.60, 98.15, 98.60, 99.15, 99.40, 100.15, 101.15,
                             102.15, 103.15, 104.15, 105.15, 106.15, 107.15,
                             108.15, 109.15, 109.20, 109.80, 110.15, 112.35,
                             117.82, 120.48, 122.71, 124.90, 126.96, 128.33,
                             130.25, 131.41, 131.77, 132.13, 135.60, 136.14,
                             137.91, 139.86, 140.23, 140.67, 142.21, 144.20,
                             145.41, 146.49, 146.70, 147.84, 148.40, 149.80,
                             150.00, 151.75, 152.41, 154.07, 155.46, 155.82,
                             156.65, 157.27, 158.50, 160.04, 161.15, 161.61,
                             164.31, 165.84, 167.40, 169.01, 170.52, 173.27,
                             173.38, 175.74, 176.22, 178.09, 179.10, 179.83,
                             180.00, 180.00, 180.92, 181.82, 182.00, 183.32,
                             184.00, 186.00, 186.24, 186.47, 187.00, 187.46,
                             187.54, 188.00, 189.00, 189.01, 189.50, 189.80,
                             190.00, 190.10, 190.20, 190.25, 190.30, 190.40,
                             190.45, 190.50, 190.53 };

    // data in bar
    Vector< real > tPvap = { 1.1626, 2.0000, 2.9064, 4.104, 5.368, 7.113,
                             14.012, 15.852, 17.918, 20.145, 22.611, 26.265,
                             25.331, 29.064, 31.064, 28.278, 34.797, 31.477,
                             38.130, 34.931, 38.663, 42.716, 47.090, 51.836,
                             56.955, 62.422, 68.288, 74.567, 81.273, 92.792,
                             98.392, 88.499, 111.109, 164.850, 199.173, 233.791,
                             269.977, 308.321, 336.270, 378.341, 405.309, 402.073,
                             423.059, 512.790, 490.332, 581.240, 588.399, 655.575,
                             686.465, 724.613, 800.026, 784.532, 892.405, 900.545,
                             950.068, 990.472, 1037.25, 1078.73, 1186.60, 1163.17,
                             1284.67, 1328.60, 1372.93, 1395.49, 1471.00, 1506.10,
                             1657.32, 1675.56, 1863.26, 1893.17, 2079.01, 2126.38,
                             2333.98, 2381.05, 2624.26, 2716.44, 2759.69, 2981.22,
                             3092.33, 3285.23, 3276.21, 3286.85, 3393.10, 3500.97,
                             3495.38, 3508.87, 3668.47, 3742.28, 3987.82, 4026.12,
                             4030.83, 4115.45, 4181.85, 4190.38, 4246.55, 4381.37,
                             4392.20, 4450.30, 4492.20, 4520.37, 4534.54, 4548.76,
                             4559.60, 4563.04, 4577.38, 4584.58, 4591.80, 4596.15 };

    // convert to Pascal
    tPvap *= 1e3;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function p_vap( T_vap )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    uint tNumPoints = tTvap.length() ;

    Vector< real > tY( tNumPoints );
    Vector< real > tY0( tNumPoints );

    for( uint k=0; k<tNumPoints; ++k )
    {
        tY0( k ) = std::log( tPvap( k ) );
        tY( k ) = std::log( tCH4.p_vap( tTvap( k ) ) );
    }

    EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-3 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function dp_vap/dT( T_vap )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // test derivative
    for( uint k=0; k<tNumPoints; ++k )
    {
        tY0( k ) =  ( tCH4.p_vap( tTvap( k ) + 0.01 ) - tCH4.p_vap( tTvap( k ) - 0.01 ) )/0.02 ;
        tY( k )  = tCH4.dpvap_dT( tTvap( k ), tCH4.p_vap( tTvap( k ) ), tCH4.pi_vap( tTvap( k ) ) );
    }
    EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-6 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function T_vap(p_vap) )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // test derivative
    for( uint k=0; k<tNumPoints; ++k )
    {
        tY( k )  =  tCH4.T_vap( tCH4.p_vap( tTvap( k ) ) );
    }
    EXPECT_NEAR( r2( tY, tTvap ), 1.0, 1e-6 );

}