/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Field dependence of the pure-metal thermal conductivity.
 *
 * Metal::set_RRR() rebuilds the lambda spline, and the spline setter resets
 * the property's dependency flags to T-only. The Kohler flags registered by
 * create_kohler() must survive that rebuild, or neither the material tool nor
 * the FEM calculator ever selects lambda( T, B, beta ).
 */

#include <cmath>

#include <gtest/gtest.h>

#include "cl_Material.hpp"
#include "cl_Material_Metal.hpp"
#include "cl_MaterialFactory.hpp"
#include "constants.hpp"
#include "fn_hust.hpp"

namespace
{
    // every pure metal that carries a Kohler magnetoresistance model
    const char * gPureMetals[] = { "copper", "aluminum", "silver", "nickel", "iron",
                                   "chromium", "indium", "lead", "tin" };
}

TEST( MetalFieldDependence, LambdaKeepsKohlerFlagsAfterSetRRR )
{
    belfem::MaterialFactory tFactory ;

    for ( const char * tLabel : gPureMetals )
    {
        // tables off: the analytic Kohler path, no database build
        belfem::Material * tMetal = tFactory.create_material( tLabel, 50.0, false );

        EXPECT_TRUE( tMetal->depends( belfem::MaterialProperty::rho,
                                      belfem::MaterialDependency::normB ) ) << tLabel ;
        EXPECT_TRUE( tMetal->depends( belfem::MaterialProperty::lambda,
                                      belfem::MaterialDependency::normB ) ) << tLabel ;
        EXPECT_TRUE( tMetal->depends( belfem::MaterialProperty::lambda,
                                      belfem::MaterialDependency::angleBxJ ) ) << tLabel ;

        delete tMetal ;
    }
}

TEST( MetalFieldDependence, LambdaFollowsWiedemannFranzInField )
{
    belfem::MaterialFactory tFactory ;

    const belfem::real T = 4.0 ;
    const belfem::real B = 20.0 ;

    for ( const char * tLabel : gPureMetals )
    {
        belfem::Material * tMetal = tFactory.create_material( tLabel, 50.0, false );

        for ( belfem::real tBeta : { 0.0, 0.5 * belfem::constant::pi } )
        {
            belfem::real tLambda0 = tMetal->lambda( T );
            belfem::real tLambdaB = tMetal->lambda( T, B, tBeta );
            belfem::real tRho0    = tMetal->rho( T );
            belfem::real tRhoB    = tMetal->rho( T, B, tBeta );

            ASSERT_TRUE( std::isfinite( tLambdaB ) ) << tLabel ;
            ASSERT_GT( tRhoB, tRho0 ) << tLabel ;

            // magnetoresistance lowers the electronic conductivity
            EXPECT_LT( tLambdaB, tLambda0 ) << tLabel << " beta = " << tBeta ;

            // lambda( T, B, beta ) = lambda( T ) * rho( T ) / rho( T, B, beta )
            EXPECT_NEAR( tLambdaB * tRhoB, tLambda0 * tRho0,
                         1e-9 * tLambda0 * tRho0 ) << tLabel << " beta = " << tBeta ;
        }

        delete tMetal ;
    }
}

// the lambda spline is built on 4 K knots with a tangent start condition; the
// tangent must be the slope of the Hust curve at the origin, L0 / rho_0, or
// the spline rings through the first knots ( copper: -79 % at 0.5 K, +3 % at
// 5 K with the reciprocal ). The 1 % band is what the knot spacing leaves for
// indium and lead between the knots; the other metals sit below 0.1 %
TEST( MetalFieldDependence, LambdaSplineMatchesHustBelowFirstKnots )
{
    belfem::MaterialFactory tFactory ;

    for ( const char * tLabel : gPureMetals )
    {
        belfem::material::Metal * tMetal = dynamic_cast< belfem::material::Metal * >(
                tFactory.create_material( tLabel, 50.0, false ) );
        ASSERT_NE( tMetal, nullptr ) << tLabel ;

        // same expression as Metal::lambda_custom, evaluated here so that the
        // spline is checked against the analytic curve it was sampled from
        const belfem::Vector< belfem::real > & p = tMetal->lambda_coefficients();
        const belfem::real tRho0 = tMetal->constant_property( belfem::MaterialProperty::rho_0 );
        const belfem::real tBeta = tRho0 / belfem::constant::L0 ;
        const belfem::real tC    = p( 7 ) / std::pow( tBeta / 0.0003, p( 0 ) );

        for ( belfem::real T : { 0.5, 1.0, 2.0, 3.0, 5.0, 6.0 } )
        {
            belfem::real tW0  = tBeta / T ;
            belfem::real tWi  = belfem::material::hust( p, T );
            belfem::real tWi0 = tC * tWi * tW0 / ( tWi + tW0 );
            belfem::real tLambdaHust = 1.0 / ( tW0 + tWi + tWi0 );

            EXPECT_NEAR( tMetal->lambda( T ), tLambdaHust, 1e-2 * tLambdaHust )
                << tLabel << " T = " << T ;
        }

        delete tMetal ;
    }
}
