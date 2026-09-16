/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Elastic moduli of the pure metals: the quasi-harmonic closure of Metal::create_mech.
 *
 * K and G soften with the volumetric thermal strain, E and nu follow. The served
 * Poisson ratio must not fall with temperature beyond data accuracy, must have zero
 * slope at 0 K and stay inside ( -1, 1/2 ); the bulk modulus must soften by the
 * amount the source data show. Copper is pinned to Ledbetter 1981
 * ( 10.1002/pssa.2210660209, Table 2 ) converted to isothermal values.
 */

#include <cmath>

#include <gtest/gtest.h>

#include "cl_Material.hpp"
#include "cl_MaterialFactory.hpp"

namespace
{
    struct MetalRatio
    {
        const char * label ;
        double       kRatio ;   // isothermal K( 4 K ) / K( 293 K ) of the source data, R1 fit 2026-09-15
    };

    // chromium is served with a constant Poisson ratio ( no source ratio ), tested separately
    const MetalRatio gMetals[] = {
        { "copper",   1.0642 },   // Ledbetter 1981, 10.1002/pssa.2210660209
        { "aluminum", 1.0847 },   // Kamm & Alers 1964, 10.1063/1.1713309
        { "silver",   1.0891 },   // Neighbours & Alers 1958, 10.1103/PhysRev.111.707
        { "indium",   1.1725 },   // Kim & Ledbetter 1998, 10.1016/S0921-5093(98)00490-0
        { "lead",     1.1669 },   // Waldorf & Alers 1962, 10.1063/1.1931149
        { "tin",      1.0970 },   // Rayne & Chandrasekhar 1960, 10.1103/PhysRev.120.1658
        { "iron",     1.0490 },   // Rayne & Chandrasekhar 1961, 10.1103/PhysRev.122.1714
        { "nickel",   1.0416 } }; // Alers, Neighbours & Sato 1960, 10.1016/0022-3697(60)90125-6
}

TEST( MetalElastic, PoissonRatioIsMonotoneAndBounded )
{
    belfem::MaterialFactory tFactory ;

    for ( const MetalRatio & tEntry : gMetals )
    {
        belfem::Material * tMetal = tFactory.create_material( tEntry.label, 50.0, false );

        const belfem::real tTmax = tMetal->constant_property( belfem::MaterialProperty::T_max );

        // no decrease larger than 1e-3 anywhere on the 4 K grid; iron's isothermal nu is flat
        belfem::real tNuMax = -1.0 ;
        for ( belfem::real T = 4.0; T <= tTmax; T += 4.0 )
        {
            const belfem::real tNu = tMetal->nu( T );
            EXPECT_GT( tNu, -1.0 ) << tEntry.label << " at " << T << " K" ;
            EXPECT_LT( tNu,  0.5 ) << tEntry.label << " at " << T << " K" ;
            tNuMax = std::max( tNuMax, tNu );
            EXPECT_LE( tNuMax - tNu, 1e-3 ) << tEntry.label << " at " << T << " K" ;
        }

        // zero slope at 0 K: alpha ~ T, T^3 there
        EXPECT_LT( std::abs( tMetal->nu( 8.0 ) - tMetal->nu( 4.0 ) ), 2e-5 ) << tEntry.label ;

        // both moduli soften on warming, the bulk modulus by the amount the data show ( 2 % band )
        EXPECT_GT( tMetal->E( 4.0 ), tMetal->E( 293.0 ) ) << tEntry.label ;
        EXPECT_NEAR( tMetal->K( 4.0 ) / tMetal->K( 293.0 ), tEntry.kRatio, 0.02 * tEntry.kRatio ) << tEntry.label ;

        delete tMetal ;
    }
}

TEST( MetalElastic, CopperMatchesLedbetter1981Isothermal )
{
    belfem::MaterialFactory tFactory ;
    belfem::Material * tCopper = tFactory.create_material( "copper", 100.0, false );

    // Ledbetter 1981, 10.1002/pssa.2210660209, Table 2: G( 5 K ) = 51.72, G( 295 K ) = 47.57,
    // B( 5 K ) = 144.46, B( 295 K ) = 139.74 GPa ( adiabatic ). G is state independent; the bulk
    // modulus is served isothermal, K_T = K_S / ( 1 + alpha_V^2 T K_S / ( rho cp ) ): 135.7 GPa at 295 K,
    // unchanged at 5 K. E and nu follow: 127.8 GPa and 0.343 at 295 K.
    EXPECT_NEAR( tCopper->G( 5.0 )   * 1e-9,  51.72, 0.5 );
    EXPECT_NEAR( tCopper->G( 295.0 ) * 1e-9,  47.57, 0.5 );
    EXPECT_NEAR( tCopper->K( 5.0 )   * 1e-9, 144.46, 1.5 );
    EXPECT_NEAR( tCopper->K( 295.0 ) * 1e-9, 135.7,  1.5 );
    EXPECT_NEAR( tCopper->E( 295.0 ) * 1e-9, 127.8,  1.0 );
    EXPECT_NEAR( tCopper->nu( 295.0 ),         0.343, 0.002 );

    // the sign that started this: nu rises from 5 K to 295 K
    EXPECT_GT( tCopper->nu( 295.0 ), tCopper->nu( 5.0 ) );

    delete tCopper ;
}

TEST( MetalElastic, ChromiumServesConstantPoissonRatio )
{
    belfem::MaterialFactory tFactory ;
    belfem::Material * tChromium = tFactory.create_material( "chromium", 50.0, false );

    // deltaG = deltaK by decision ( the spin density wave anomalies are not represented ):
    // nu = 0.2371 at every temperature, E and K still soften
    EXPECT_NEAR( tChromium->nu( 4.0 ),   0.2371, 1e-3 );
    EXPECT_NEAR( tChromium->nu( 293.0 ), 0.2371, 1e-3 );
    EXPECT_LT( std::abs( tChromium->nu( 293.0 ) - tChromium->nu( 4.0 ) ), 1e-4 );
    EXPECT_GT( tChromium->E( 4.0 ), tChromium->E( 293.0 ) );

    delete tChromium ;
}

TEST( MetalElastic, FormulaAlloyInheritsBoundedPoissonRatio )
{
    belfem::MaterialFactory tFactory ;
    belfem::Material * tSolder = tFactory.create_material( "Sn60Pb40", 50.0, false );

    for ( belfem::real T : { 77.0, 293.0 } )
    {
        EXPECT_GT( tSolder->nu( T ), -1.0 ) << "at " << T << " K" ;
        EXPECT_LT( tSolder->nu( T ),  0.5 ) << "at " << T << " K" ;
        EXPECT_GT( tSolder->E( T ),   0.0 ) << "at " << T << " K" ;
    }

    delete tSolder ;
}
