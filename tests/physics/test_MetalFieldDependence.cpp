/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Pure-metal resistivity and thermal conductivity around Metal::set_RRR().
 *
 * MetalFieldDependence: set_RRR() rebuilds the lambda spline, and the spline
 * setter resets the property's dependency flags to T-only. The Kohler flags
 * registered by create_kohler() must survive that rebuild, or neither the
 * material tool nor the FEM calculator ever selects lambda( T, B, beta ).
 *
 * MetalRRR: set_RRR() solves ( rho_i( T_ref ) + rho_0 ) / rho_0 = RRR for the
 * residual resistivity rho_0 by a bracketed regula-falsi iteration; the stored rho_0 must
 * satisfy that identity for the RRR the material was created with and again
 * after the RRR is changed on the same object.
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

// the ratio the regula-falsi residual in set_RRR() drives to zero, evaluated
// on the rho_0 the material stores afterwards
TEST( MetalRRR, SecantSatisfiesTheRRRIdentity )
{
    belfem::MaterialFactory tFactory ;

    for ( const char * tLabel : gPureMetals )
    {
        SCOPED_TRACE( tLabel );

        belfem::material::Metal * tMetal = dynamic_cast< belfem::material::Metal * >(
                tFactory.create_material( tLabel, 50.0, false ) );
        ASSERT_NE( tMetal, nullptr );

        const belfem::real tTref = tMetal->constant_property( belfem::MaterialProperty::T_ref_rho_i );
        EXPECT_TRUE( std::isfinite( tTref ) );
        EXPECT_GT( tTref, 0.0 );

        // RRR from the factory argument
        belfem::real tRho0 = tMetal->constant_property( belfem::MaterialProperty::rho_0 );
        EXPECT_TRUE( std::isfinite( tRho0 ) );
        EXPECT_GT( tRho0, 0.0 );
        EXPECT_NEAR( ( tMetal->rho_i_custom( tTref ) + tRho0 ) / tRho0, 50.0, 1e-9 );

        // RRR changed on the same object
        tMetal->set_RRR( 120.0 );
        tRho0 = tMetal->constant_property( belfem::MaterialProperty::rho_0 );
        EXPECT_TRUE( std::isfinite( tRho0 ) );
        EXPECT_GT( tRho0, 0.0 );
        EXPECT_NEAR( ( tMetal->rho_i_custom( tTref ) + tRho0 ) / tRho0, 120.0, 1e-9 );

        delete tMetal ;
    }
}

// Lead and tin switch the transverse Kohler term from a low-field cubic to
// exp( q1( ln BS ) ) at the point where the cubic matches value, slope and
// curvature of the power law. The switch must carry that matching point: the
// far-point checks at x / 4 and 2 x pin which branch is taken on each side,
// against an oracle built here from the published q1. Continuity at x alone is
// no mutation gate, since a switch outside the sweep selects one smooth branch
// throughout and both branches are increasing there
TEST( MetalFieldDependence, KohlerTransverseBranchesMeetAtTheStoredSwitch )
{
    struct KohlerTrans
    {
        const char * mLabel ;
        belfem::real mExponent ;   // q1( 0 ): A = exp( mOffset ) * BS^mExponent
        belfem::real mOffset ;     // q1( 1 )
        belfem::real mLogSwitch ;  // ln of the matching point
    };

    const KohlerTrans tCases[] = { { "lead", 1.458,  -9.7936, 3.0 },
                                   { "tin",  1.8073, -12.545, 8.4 } };

    const belfem::real tBeta = 0.5 * belfem::constant::pi ;

    belfem::MaterialFactory tFactory ;

    for ( const KohlerTrans & tCase : tCases )
    {
        SCOPED_TRACE( tCase.mLabel );

        belfem::material::Metal * tMetal = dynamic_cast< belfem::material::Metal * >(
                tFactory.create_material( tCase.mLabel, 50.0, false ) );
        ASSERT_NE( tMetal, nullptr );

        const belfem::real k = tCase.mExponent ;
        const belfem::real x = std::exp( tCase.mLogSwitch );

        // power law and its first two derivatives at the matching point
        const belfem::real F  = std::exp( tCase.mOffset ) * std::pow( x, k );
        const belfem::real F1 = k * F / x ;
        const belfem::real F2 = k * ( k - 1.0 ) * F / ( x * x );

        // a t^3 + b t^2 + c t with value, slope and curvature F, F1, F2 at t = x
        const belfem::real a = ( F / x - F1 + 0.5 * F2 * x ) / ( x * x );
        const belfem::real b = 0.5 * F2 - 3.0 * a * x ;
        const belfem::real c = F / x - a * x * x - b * x ;

        auto tHigh  = [ & ]( const belfem::real t ) { return std::exp( tCase.mOffset ) * std::pow( t, k ); };
        auto tCubic = [ & ]( const belfem::real t ) { return ( ( a * t + b ) * t + c ) * t ; };

        const belfem::real tLow = 0.25 * x ;
        const belfem::real tUp  = 2.0 * x ;

        // the far points discriminate: the two branches differ there
        EXPECT_GT( std::abs( tCubic( tLow ) - tHigh( tLow ) ), 1e-2 * tHigh( tLow ) );
        EXPECT_GT( std::abs( tCubic( tUp ) - tHigh( tUp ) ), 1e-2 * tHigh( tUp ) );

        const belfem::real tKLow = tMetal->kohler( tLow, 1.0, tBeta );
        const belfem::real tKUp  = tMetal->kohler( tUp, 1.0, tBeta );
        EXPECT_NEAR( tKLow, tCubic( tLow ), 1e-10 * std::abs( tCubic( tLow ) ) );
        EXPECT_NEAR( tKUp, tHigh( tUp ), 1e-10 * tHigh( tUp ) );

        // control: the branches meet at the switch. Each side is compared with the
        // power law at its own point; the C2 match leaves an O( delta^3 ) gap,
        // while the slope alone separates the two points by O( delta )
        for ( belfem::real t : { x * ( 1.0 - 1e-4 ), x * ( 1.0 + 1e-4 ) } )
        {
            EXPECT_NEAR( tMetal->kohler( t, 1.0, tBeta ), tHigh( t ), 1e-6 * tHigh( t ) )
                << "BS = " << t ;
        }

        // longitudinal smoke control, below the first and above the last switch;
        // tin's low-field longitudinal cubic is negative below about 4 T
        for ( belfem::real B : { 6.0, 1e7 } )
        {
            const belfem::real tKLong = tMetal->kohler( B, 1.0, 0.0 );
            EXPECT_TRUE( std::isfinite( tKLong ) ) << "B = " << B ;
            EXPECT_GT( tKLong, 0.0 ) << "B = " << B ;
        }

        delete tMetal ;
    }
}
