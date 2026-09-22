/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the tensor-grid lookup table in src/physics/database. The
 * fixtures give each axis its own node step, so a derivative scaled by the
 * wrong axis is visible, and drive the field with a polynomial the element
 * space reproduces exactly, so value and partials have a closed-form oracle.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Database.hpp"
#include "commtools.hpp"

namespace
{
    using belfem::real;
    using belfem::index_t;

    // exactly reproduced fields, so the tolerance only covers roundoff
    const real tTolLinear = 1e-9;
    const real tTolQuad   = 1e-8;

    // -- trilinear, the full HEX8 space -------------------------------------

    real hex8_f( const real aX, const real aY, const real aZ )
    {
        return 1.0 + 2.0 * aX - 3.0 * aY + 5.0 * aZ
             + 0.5 * aX * aY - 0.25 * aY * aZ + 0.75 * aX * aZ
             + 0.125 * aX * aY * aZ;
    }

    real hex8_dfdx( const real, const real aY, const real aZ )
    {
        return 2.0 + 0.5 * aY + 0.75 * aZ + 0.125 * aY * aZ;
    }

    real hex8_dfdy( const real aX, const real, const real aZ )
    {
        return -3.0 + 0.5 * aX - 0.25 * aZ + 0.125 * aX * aZ;
    }

    real hex8_dfdz( const real aX, const real aY, const real )
    {
        return 5.0 + 0.75 * aX - 0.25 * aY + 0.125 * aX * aY;
    }

    // -- bilinear, the full QUAD4 space --------------------------------------

    real quad4_f( const real aX, const real aY )
    {
        return 1.0 + 2.0 * aX - 3.0 * aY + 0.5 * aX * aY;
    }

    real quad4_dfdx( const real, const real aY ) { return  2.0 + 0.5 * aY; }
    real quad4_dfdy( const real aX, const real ) { return -3.0 + 0.5 * aX; }

    // -- separable quadratics, inside the HEX27 and QUAD9 spaces -------------

    real fx ( const real aX ) { return 1.0 + aX + 0.5 * aX * aX; }
    real dfx( const real aX ) { return 1.0 + aX; }
    real fy ( const real aY ) { return 2.0 - aY + 0.25 * aY * aY; }
    real dfy( const real aY ) { return -1.0 + 0.5 * aY; }
    real fz ( const real aZ ) { return 1.0 + 2.0 * aZ - 0.5 * aZ * aZ; }
    real dfz( const real aZ ) { return 2.0 - aZ; }
}

// =============================================================================
// order 1
// =============================================================================

TEST( Database, Hex8DerivativesUseEachAxisStep )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // node steps 0.5, 2.0 and 0.125: a derivative scaled by axis 0 instead of
    // its own is off by a factor of 4 in y and 1/4 in z
    belfem::Mesh tMesh( 1,
        belfem::Vector< index_t >( { 4, 3, 5 } ),
        belfem::Vector< real >( { 0.5, 2.0, 0.125 } ),
        belfem::Vector< real >( { -1.0, 2.0, 0.25 } ) );

    belfem::Vector< real > & tField = tMesh.create_field( "f" );

    for ( index_t i = 0; i < 4; ++i )
        for ( index_t j = 0; j < 3; ++j )
            for ( index_t k = 0; k < 5; ++k )
            {
                belfem::mesh::Node * tNode = tMesh.node( i, j, k );
                tField( tNode->index() ) =
                    hex8_f( tNode->x(), tNode->y(), tNode->z() );
            }

    belfem::Database tDatabase( &tMesh, "f", false );

    // interior points, all coordinates non-zero, none on an element boundary,
    // and no two parametric coordinates of a point equal
    const real tX[ 3 ] = { -0.85, 0.375,     -0.175   };
    const real tY[ 3 ] = {  3.6,  4.75,       2.3     };
    const real tZ[ 3 ] = {  0.325, 0.734375,  0.53125 };

    for ( belfem::uint p = 0; p < 3; ++p )
    {
        const real x = tX[ p ];
        const real y = tY[ p ];
        const real z = tZ[ p ];

        EXPECT_NEAR( tDatabase.evaluate( x, y, z ),
                     hex8_f( x, y, z ), tTolLinear ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivx( x, y, z ),
                     hex8_dfdx( x, y, z ), tTolLinear ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivy( x, y, z ),
                     hex8_dfdy( x, y, z ), tTolLinear ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivz( x, y, z ),
                     hex8_dfdz( x, y, z ), tTolLinear ) << "point " << p;
    }
}

TEST( Database, Quad4DerivativesUseEachAxisStep )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // the 2-D overloads carry their own copy of the per-axis scaling
    belfem::Mesh tMesh( 1,
        belfem::Vector< index_t >( { 4, 6 } ),
        belfem::Vector< real >( { 0.5, 3.0 } ),
        belfem::Vector< real >( { 2.0, -1.0 } ) );

    belfem::Vector< real > & tField = tMesh.create_field( "f" );

    // node( i, j ): the three-argument overload asserts a 3-D mesh
    for ( index_t i = 0; i < 4; ++i )
        for ( index_t j = 0; j < 6; ++j )
        {
            belfem::mesh::Node * tNode = tMesh.node( i, j );
            tField( tNode->index() ) = quad4_f( tNode->x(), tNode->y() );
        }

    belfem::Database tDatabase( &tMesh, "f", false );

    const real tX[ 3 ] = { 2.15, 3.375,  2.825 };
    const real tY[ 3 ] = { 1.4,  12.125, 5.45  };

    for ( belfem::uint p = 0; p < 3; ++p )
    {
        const real x = tX[ p ];
        const real y = tY[ p ];

        EXPECT_NEAR( tDatabase.evaluate( x, y ),
                     quad4_f( x, y ), tTolLinear ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivx( x, y ),
                     quad4_dfdx( x, y ), tTolLinear ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivy( x, y ),
                     quad4_dfdy( x, y ), tTolLinear ) << "point " << p;
    }
}

// =============================================================================
// order 2
// =============================================================================

TEST( Database, Hex27ReproducesATriquadraticExactly )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // a separable triquadratic lies in the space HEX27 spans, so value and all
    // three partials are exact. Every one of the 27 nodal weights is non-zero
    // at a generic point, so a wrong entry in any of the four tables -- the
    // value table or one of the three derivative tables -- shows up here.
    belfem::Mesh tMesh( 2,
        belfem::Vector< index_t >( { 5, 7, 3 } ),
        belfem::Vector< real >( { 0.5, 0.25, 2.0 } ),
        belfem::Vector< real >( { 1.0, -2.0, 0.0 } ) );

    belfem::Vector< real > & tField = tMesh.create_field( "f" );

    for ( index_t i = 0; i < 5; ++i )
        for ( index_t j = 0; j < 7; ++j )
            for ( index_t k = 0; k < 3; ++k )
            {
                belfem::mesh::Node * tNode = tMesh.node( i, j, k );
                tField( tNode->index() ) =
                    fx( tNode->x() ) * fy( tNode->y() ) * fz( tNode->z() );
            }

    belfem::Database tDatabase( &tMesh, "f", false );

    const real tX[ 3 ] = { 1.3,  2.75,    1.65   };
    const real tY[ 3 ] = { -1.6, -0.8125, -1.425 };
    const real tZ[ 3 ] = { 2.4,  3.5,     1.0    };

    for ( belfem::uint p = 0; p < 3; ++p )
    {
        const real x = tX[ p ];
        const real y = tY[ p ];
        const real z = tZ[ p ];

        EXPECT_NEAR( tDatabase.evaluate( x, y, z ),
                     fx( x ) * fy( y ) * fz( z ), tTolQuad ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivx( x, y, z ),
                     dfx( x ) * fy( y ) * fz( z ), tTolQuad ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivy( x, y, z ),
                     fx( x ) * dfy( y ) * fz( z ), tTolQuad ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivz( x, y, z ),
                     fx( x ) * fy( y ) * dfz( z ), tTolQuad ) << "point " << p;
    }
}

TEST( Database, Quad9ReproducesABiquadraticExactly )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    belfem::Mesh tMesh( 2,
        belfem::Vector< index_t >( { 5, 5 } ),
        belfem::Vector< real >( { 0.5, 2.0 } ),
        belfem::Vector< real >( { 1.0, -2.0 } ) );

    belfem::Vector< real > & tField = tMesh.create_field( "f" );

    for ( index_t i = 0; i < 5; ++i )
        for ( index_t j = 0; j < 5; ++j )
        {
            belfem::mesh::Node * tNode = tMesh.node( i, j );
            tField( tNode->index() ) = fx( tNode->x() ) * fy( tNode->y() );
        }

    belfem::Database tDatabase( &tMesh, "f", false );

    const real tX[ 3 ] = { 1.3, 2.75, 1.65 };
    const real tY[ 3 ] = { 1.2, 3.5,  2.6  };

    for ( belfem::uint p = 0; p < 3; ++p )
    {
        const real x = tX[ p ];
        const real y = tY[ p ];

        EXPECT_NEAR( tDatabase.evaluate( x, y ),
                     fx( x ) * fy( y ), tTolQuad ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivx( x, y ),
                     dfx( x ) * fy( y ), tTolQuad ) << "point " << p;
        EXPECT_NEAR( tDatabase.evaluate_derivy( x, y ),
                     fx( x ) * dfy( y ), tTolQuad ) << "point " << p;
    }
}
