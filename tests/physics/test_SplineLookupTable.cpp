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
