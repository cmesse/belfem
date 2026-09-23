/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * The HEX8TS quadrature volume: sum_w() is the reference measure of
 * [-1,1]^3 ( 8 ), abs_det_J() is 0.5 |dX/dxi x dX/deta| times the layer
 * thickness, and the thickness comes from the block, not from the node
 * coordinates. The HEX8TS twin of InterfaceOrientation.Penta6TsVolumeIdentity.
 *
 * Evaluated at the element center: for a planar bilinear quad
 * det J( xi, eta ) = a0 + a1 xi + a2 eta, and its integral over [-1,1]^2
 * is 4 a0 = A, so the center value is exactly A / 4 even for the fixture's
 * non-parallelogram default quad.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Matrix.hpp"
#include "support/cl_TS_TestStack.hpp"
#include "nedelec/cl_EF_HEX8TS.hpp"

using namespace belfem;

namespace
{
    const real gEpsilon = 1.0e-12 ;

    // shoelace area of the four mid-surface facet nodes ( z = 0 plane,
    // counterclockwise from +z )
    real
    facet_area( const fem::Element * aElement )
    {
        const mesh::Facet * tFacet = aElement->facet();
        real tTwiceArea = 0.0 ;
        for ( uint k = 0; k < 4; ++k )
        {
            const mesh::Node * tA = tFacet->node( k );
            const mesh::Node * tB = tFacet->node( ( k + 1 ) % 4 );
            tTwiceArea += tA->x( 0 ) * tB->x( 1 ) - tB->x( 0 ) * tA->x( 1 );
        }
        return 0.5 * tTwiceArea ;
    }

    // abs_det_J() at the element center; E() is what refreshes the
    // per-point Jacobian
    real
    center_abs_det_j( fem::EF_HEX8TS & aEF, fem::Element * aElement )
    {
        Matrix< real > tXi( 3, 1, 0.0 );
        aEF.precompute( tXi );
        aEF.link( aElement );
        aEF.E( 0 );
        return aEF.abs_det_J();
    }
}

//------------------------------------------------------------------------------

TEST( Hex8TsVolume, VolumeIdentityAtTwoThicknesses )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    for ( const real tThickness : { 0.1, 0.4 } )
    {
        fem::test::TS_TestPrism tPrism( ElementType::HEX8TS, tThickness );

        fem::EF_HEX8TS tEF ;
        const real tDetJ = center_abs_det_j( tEF, tPrism.element() );
        const real tArea = facet_area( tPrism.element() );

        ASSERT_GT( tArea, 0.0 );
        EXPECT_NEAR( tEF.sum_w(), 8.0, gEpsilon ) << "thickness " << tThickness ;
        EXPECT_NEAR( tDetJ * tEF.sum_w(), tArea * tThickness, gEpsilon )
            << "thickness " << tThickness ;
    }
}

//------------------------------------------------------------------------------

// one thickness cannot tell a thickness factor from a constant: four times
// the thickness must give four times the Jacobian on the same facet
TEST( Hex8TsVolume, DetJScalesWithTheBlockThickness )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::test::TS_TestPrism tThin( ElementType::HEX8TS, 0.1 );
    fem::test::TS_TestPrism tThick( ElementType::HEX8TS, 0.4 );

    fem::EF_HEX8TS tEFThin ;
    fem::EF_HEX8TS tEFThick ;
    const real tDetThin  = center_abs_det_j( tEFThin, tThin.element() );
    const real tDetThick = center_abs_det_j( tEFThick, tThick.element() );

    ASSERT_GT( tDetThin, 0.0 );
    EXPECT_NEAR( tDetThick / tDetThin, 4.0, gEpsilon );
    EXPECT_EQ( tEFThin.sum_w(), tEFThick.sum_w() );
}
