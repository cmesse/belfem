/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * The convenience overload carries BOTH end slopes to the spline.
 *
 * Material::create_spline( property, dYdX0, dYdX1 ) turns each finite slope
 * into a Tangent end condition and forwards the value through the virtual
 * SplineLookupTable::create_spline into the Spline constructor, whose own
 * arguments default to zero. A dropped forward would therefore not fail to
 * build; it would clamp an end to zero silently. Production passes the START
 * slope ( YBCO and Metal hand L0 / rho_0 to lambda, YBCO 0 to E, nu and cp )
 * but no caller passes the END slope, so these cases are the only exercise of
 * that argument.
 *
 * The discriminator is uniqueness: a clamped cubic spline that interpolates a
 * cubic with the cubic's own end slopes reproduces it identically, because
 * the cubic satisfies every equation of the system. Any other end value
 * perturbs the slopes, and the control below shows it does.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "cl_Material_SplineLookupTable.hpp"
#include "cl_Spline.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Vector.hpp"

namespace
{
    // a cubic with no zero coefficient, positive over [ 0, 40 ]
    const belfem::real tA = 412.7;
    const belfem::real tB = 3.19;
    const belfem::real tC = -0.0417;
    const belfem::real tD = 2.63e-4;

    const belfem::real tTmax = 40.0;

    belfem::real
    cubic( const belfem::real aT )
    {
        return ( ( tD * aT + tC ) * aT + tB ) * aT + tA;
    }

    belfem::real
    dcubic( const belfem::real aT )
    {
        return ( 3.0 * tD * aT + 2.0 * tC ) * aT + tB;
    }

    belfem::real
    ddcubic( const belfem::real aT )
    {
        return 6.0 * tD * aT + 2.0 * tC;
    }

    // a lookup table whose raw cp law is that cubic. The end slope is a
    // constructor argument: the default NaN selects a Parabolic end, which is
    // the control
    class CubicLaw : public belfem::material::SplineLookupTable
    {
    public:

        explicit CubicLaw( const belfem::real aEndSlope
                           = std::numeric_limits< belfem::real >::quiet_NaN() )
            : SplineLookupTable( belfem::MaterialType::NonMetal, true )
        {
            this->set_label( "CubicLaw" );
            this->set_constant( belfem::MaterialProperty::T_max, tTmax );
            this->set_have( belfem::MaterialProperty::cp );
            this->create_spline( belfem::MaterialProperty::cp, tB, aEndSlope );
        }

    protected:

        belfem::real
        cp_custom( const belfem::real T ) const override
        {
            return cubic( T );
        }
    };

    // the knots are 4 K apart from 0; none of these sits on one
    const std::array< belfem::real, 5 > tSamples = { 1.3, 7.9, 18.6, 27.1, 38.4 };
}

TEST( SplineLookupTable, ConvenienceOverloadForwardsBothEndSlopes )
{
    const CubicLaw tMaterial( dcubic( tTmax ) );

    for ( const belfem::real tT : tSamples )
    {
        EXPECT_NEAR( tMaterial.cp( tT ), cubic( tT ), 1e-10 * std::abs( cubic( tT ) ) )
                << "T = " << tT;

        EXPECT_NEAR( tMaterial.dcpdT( tT ), dcubic( tT ), 1e-9 * std::abs( dcubic( tT ) ) )
                << "T = " << tT;

        EXPECT_NEAR( tMaterial.d2cpdT2( tT ), ddcubic( tT ), 1e-8 * std::abs( ddcubic( tT ) ) )
                << "T = " << tT;
    }

    // both ends are clamped, so the interpolant is the cubic there too
    EXPECT_NEAR( tMaterial.dcpdT( 0.0 ), tB, 1e-9 * std::abs( tB ) );
    EXPECT_NEAR( tMaterial.dcpdT( tTmax ), dcubic( tTmax ), 1e-9 * std::abs( dcubic( tTmax ) ) );
}

TEST( SplineLookupTable, ParabolicEndIsNotTheCubic )
{
    const CubicLaw tMaterial;

    // the end slope is the only thing that changed, so a forwarded default of
    // zero - or a forward that never happens - could not be told from the
    // true slope if this case were green
    EXPECT_GT( std::abs( tMaterial.dcpdT( tTmax ) - dcubic( tTmax ) ) / std::abs( dcubic( tTmax ) ),
               1e-4 );

    EXPECT_GT( std::abs( tMaterial.cp( 38.4 ) - cubic( 38.4 ) ) / std::abs( cubic( 38.4 ) ),
               1e-6 );

    // the start is still clamped, and the perturbation decays away from the
    // end it was introduced at
    EXPECT_NEAR( tMaterial.dcpdT( 0.0 ), tB, 1e-9 * std::abs( tB ) );
}

/*
 * set_spline takes ownership of the spline it is handed and points the
 * property's dispatch at it. Called again with the default nullptr it keeps
 * the stored spline and only refreshes the dispatch; called with a different
 * spline it deletes the stored one and keeps the new one.
 *
 * Not observable here: that the replaced spline is freed, and that setting
 * the SAME pointer again keeps it alive rather than deleting and re-storing
 * it. Both need a memory checker.
 */
namespace
{
    // a natural cubic spline through aFunction on [ 1, 10 ], 20 points; the
    // caller owns it
    template< typename F >
    belfem::Spline *
    make_spline( F && aFunction )
    {
        const belfem::uint tN = 20 ;
        belfem::Vector< belfem::real > tX( tN );
        belfem::Vector< belfem::real > tY( tN );
        for ( belfem::uint k = 0; k < tN; ++k )
        {
            tX( k ) = 1.0 + 9.0 * k / ( tN - 1 );
            tY( k ) = aFunction( tX( k ) );
        }

        belfem::SpMatrix tHelp ;
        belfem::spline::create_helpmatrix( tN, tX( 1 ) - tX( 0 ), tHelp );

        return new belfem::Spline( tX, tY, tHelp );
    }

    const belfem::real gProbeT = 4.3 ;
}

TEST( SplineLookupTable, SetSplineStoresAndDispatches )
{
    belfem::material::SplineLookupTable tTable( belfem::MaterialType::PureMetal );

    EXPECT_FALSE( tTable.have( belfem::MaterialProperty::cp ) );

    belfem::Spline * tSplineA = make_spline( []( const belfem::real aX ) { return aX * aX ; } );
    tTable.set_spline( belfem::MaterialProperty::cp, tSplineA );

    EXPECT_TRUE( tTable.have( belfem::MaterialProperty::cp ) );
    ASSERT_EQ( tTable.spline( belfem::MaterialProperty::cp ), tSplineA );
    EXPECT_EQ( tTable.cp( gProbeT ), tSplineA->eval( gProbeT ) );
}

TEST( SplineLookupTable, NullRefreshKeepsTheSpline )
{
    belfem::material::SplineLookupTable tTable( belfem::MaterialType::PureMetal );

    belfem::Spline * tSplineA = make_spline( []( const belfem::real aX ) { return aX * aX ; } );
    tTable.set_spline( belfem::MaterialProperty::cp, tSplineA );
    const belfem::real tBefore = tTable.cp( gProbeT );

    tTable.set_spline( belfem::MaterialProperty::cp );

    // checked before any evaluation: a refresh that dropped the spline would
    // leave a null pointer behind the dispatch
    ASSERT_EQ( tTable.spline( belfem::MaterialProperty::cp ), tSplineA );
    EXPECT_TRUE( tTable.have( belfem::MaterialProperty::cp ) );
    EXPECT_EQ( tTable.cp( gProbeT ), tBefore );
}

TEST( SplineLookupTable, ReplacementTakesOwnership )
{
    belfem::material::SplineLookupTable tTable( belfem::MaterialType::PureMetal );

    belfem::Spline * tSplineA = make_spline( []( const belfem::real aX ) { return aX * aX ; } );
    tTable.set_spline( belfem::MaterialProperty::cp, tSplineA );
    const belfem::real tOld = tTable.cp( gProbeT );

    belfem::Spline * tSplineB = make_spline( []( const belfem::real aX ) { return 3.0 * aX + 1.0 ; } );
    const belfem::real tNew = tSplineB->eval( gProbeT );

    // the two laws differ at the probe, so the value below tells them apart
    ASSERT_GT( std::abs( tNew - tOld ), 1.0 );

    // tSplineA is deleted by this call and must not be touched again
    tTable.set_spline( belfem::MaterialProperty::cp, tSplineB );

    ASSERT_EQ( tTable.spline( belfem::MaterialProperty::cp ), tSplineB );
    EXPECT_EQ( tTable.cp( gProbeT ), tNew );
}
