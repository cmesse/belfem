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
 * The automatic integration order. IntegrationData::populate( 0 ) routes to
 * auto_integration_order( type ), which maps the element's interpolation
 * order to a rule degree: constant 1, linear 4, serendipity and quadratic
 * 7, cubic 10, quartic 13, quintic 16. The rule of degree d integrates
 * every monomial of total degree <= d exactly ( IntegrationExactness locks
 * the tables ); here the chain from element type to rule is checked.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"
#include "meshtools.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "en_IntegrationScheme.hpp"
#include "cl_IF_IntegrationData.hpp"

using namespace belfem;

namespace
{
    // every type with an interpolation order that the fem suite builds
    const ElementType gTypes[] = {
            ElementType::LINE2,   ElementType::LINE3,   ElementType::LINE4,
            ElementType::TRI3,    ElementType::TRI6,    ElementType::TRI10,
            ElementType::QUAD4,   ElementType::QUAD8,   ElementType::QUAD9,   ElementType::QUAD16,
            ElementType::QUAD4TS, ElementType::QUAD9TS,
            ElementType::TET4,    ElementType::TET10,   ElementType::TET20,
            ElementType::PENTA6,  ElementType::PENTA15, ElementType::PENTA18,
            ElementType::PENTA6TS, ElementType::PENTA18TS,
            ElementType::HEX8,    ElementType::HEX20,   ElementType::HEX27,   ElementType::HEX64,
            ElementType::HEX8TS,  ElementType::HEX8TB };

    uint
    expected_degree( const InterpolationOrder aOrder )
    {
        switch ( aOrder )
        {
            case InterpolationOrder::CONSTANT    : return 1 ;
            case InterpolationOrder::LINEAR      : return 4 ;
            case InterpolationOrder::SERENDIPITY :
            case InterpolationOrder::QUADRATIC   : return 7 ;
            case InterpolationOrder::CUBIC       : return 10 ;
            case InterpolationOrder::QUARTIC     : return 13 ;
            case InterpolationOrder::QUINTIC     : return 16 ;
            default                              : return 0 ;
        }
    }

    real
    factorial( const uint n )
    {
        real f = 1.0 ;
        for ( uint k = 2; k <= n; ++k ) f *= k ;
        return f ;
    }

    // integral of t^a over [-1,1]
    real
    line_moment( const uint a )
    {
        return ( a % 2 == 0 ) ? 2.0 / ( a + 1 ) : 0.0 ;
    }

    // integral of x^a y^b over the unit triangle
    real
    tri_moment( const uint a, const uint b )
    {
        return factorial( a ) * factorial( b ) / factorial( a + b + 2 );
    }

    // integral of x^a y^b z^c over the unit tetrahedron
    real
    tet_moment( const uint a, const uint b, const uint c )
    {
        return factorial( a ) * factorial( b ) * factorial( c ) / factorial( a + b + c + 3 );
    }

    // reference domains as in IntegrationExactness: LINE [-1,1], QUAD and
    // HEX [-1,1]^d, TRI and TET the unit simplex, PENTA triangle x [-1,1]
    real
    exact_moment( const GeometryType aGeometry, const uint a, const uint b, const uint c )
    {
        switch ( aGeometry )
        {
            case GeometryType::LINE  : return line_moment( a );
            case GeometryType::QUAD  : return line_moment( a ) * line_moment( b );
            case GeometryType::HEX   : return line_moment( a ) * line_moment( b ) * line_moment( c );
            case GeometryType::TRI   : return tri_moment( a, b );
            case GeometryType::TET   : return tet_moment( a, b, c );
            case GeometryType::PENTA : return tri_moment( a, b ) * line_moment( c );
            default                  : return BELFEM_QUIET_NAN ;
        }
    }

    uint
    dimension_of( const GeometryType aGeometry )
    {
        switch ( aGeometry )
        {
            case GeometryType::LINE : return 1 ;
            case GeometryType::TRI  :
            case GeometryType::QUAD : return 2 ;
            default                 : return 3 ;
        }
    }
}

//------------------------------------------------------------------------------

TEST( IntegrationAuto, AutoOrderFollowsTheInterpolationOrderTable )
{
    for ( const ElementType tType : gTypes )
    {
        const InterpolationOrder tOrder = mesh::interpolation_order( tType );
        ASSERT_NE( tOrder, InterpolationOrder::UNDEFINED ) << "type " << ( int ) tType ;

        uint tDegree = 0 ;
        ASSERT_NO_THROW( tDegree = auto_integration_order( tType ) ) << "type " << ( int ) tType ;
        EXPECT_EQ( tDegree, expected_degree( tOrder ) ) << "type " << ( int ) tType ;

        // the rule must at least integrate a mass matrix of the element's own
        // basis exactly
        EXPECT_GE( tDegree, 2 * mesh::interpolation_order_numeric( tType ) )
            << "type " << ( int ) tType ;
    }
}

TEST( IntegrationAuto, PopulateOrderZeroIsExactToTwiceTheInterpolationDegree )
{
    for ( const ElementType tType : gTypes )
    {
        const GeometryType tGeometry = mesh::geometry_type( tType );
        const uint tDegree = auto_integration_order( tType );

        fem::IntegrationData tData( tType );
        tData.populate( 0 );

        // the same rule as asked for explicitly
        Vector< real > tW ;
        Matrix< real > tP ;
        intpoints( IntegrationScheme::GAUSS, tGeometry, tDegree, tW, tP );

        const Vector< real > & tDataW = tData.weights() ;
        const Matrix< real > & tDataP = tData.points() ;

        ASSERT_EQ( tDataW.length(), tW.length() ) << "type " << ( int ) tType ;
        ASSERT_GT( tDataW.length(), 0u ) << "type " << ( int ) tType ;

        for ( uint k = 0; k < tW.length(); ++k )
        {
            EXPECT_EQ( tDataW( k ), tW( k ) ) << "type " << ( int ) tType << " weight " << k ;
            for ( uint i = 0; i < tP.n_rows(); ++i )
            {
                EXPECT_EQ( tDataP( i, k ), tP( i, k ) ) << "type " << ( int ) tType << " point " << k ;
            }
        }

        // every monomial of TOTAL degree up to twice the interpolation
        // degree, the reference measure being degree 0. That is the whole
        // mass matrix on simplices; on tensor-product elements the product
        // of two basis functions reaches a higher total degree ( QUAD9:
        // xi^4 eta^4 ), which this does not claim
        const uint tDim = dimension_of( tGeometry );
        const uint tMax = 2 * mesh::interpolation_order_numeric( tType );

        for ( uint d = 0; d <= tMax; ++d )
        {
            for ( uint a = 0; a <= d; ++a )
            {
                for ( uint b = 0; b <= d - a; ++b )
                {
                    const uint c = d - a - b ;
                    if ( tDim < 3 && c > 0 ) continue ;
                    if ( tDim < 2 && b > 0 ) continue ;

                    real tQ = 0.0 ;
                    for ( uint k = 0; k < tDataW.length(); ++k )
                    {
                        real tTerm = tDataW( k ) * std::pow( tDataP( 0, k ), a );
                        if ( tDim > 1 ) tTerm *= std::pow( tDataP( 1, k ), b );
                        if ( tDim > 2 ) tTerm *= std::pow( tDataP( 2, k ), c );
                        tQ += tTerm ;
                    }

                    const real tExact = exact_moment( tGeometry, a, b, c );
                    EXPECT_NEAR( tQ, tExact, 1e-12 * std::max( 1.0, std::abs( tExact ) ) )
                        << "type " << ( int ) tType
                        << " monomial x^" << a << " y^" << b << " z^" << c ;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

// Contract, red on the tree today: every element type the fem suite builds
// must have an automatic integration order. mesh::interpolation_order()
// has no branch for the pyramids and no QUARTIC branch at all ( TET35 ),
// so auto_integration_order() refuses them with BELFEM_ERROR. The fix is on
// the mesh side; this case turns green, and starts covering the QUARTIC
// arm of the table, the day it lands.
TEST( IntegrationAuto, DISABLED_AutoOrderCoversPyramidsAndQuarticTets )
{
    const ElementType tRed[] = {
            ElementType::PYRA5, ElementType::PYRA13, ElementType::PYRA14,
            ElementType::TET35 };

    const uint tExpected[] = { 4, 7, 7, 13 };

    for ( uint k = 0; k < 4; ++k )
    {
        uint tDegree = 0 ;
        EXPECT_NO_THROW( tDegree = auto_integration_order( tRed[ k ] ) )
            << "type " << ( int ) tRed[ k ] ;
        EXPECT_EQ( tDegree, tExpected[ k ] ) << "type " << ( int ) tRed[ k ] ;
    }
}
