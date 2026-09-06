/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_gauss_tet10.hpp"
#include "en_IntegrationScheme.hpp"

using namespace belfem;

namespace
{
    // The contract of intpoints(): a rule requested for order p integrates every
    // monomial of total degree <= p exactly on the reference domain the point
    // rows are written in. This test locks that contract for every geometry
    // and every order the dispatch provides, so a table whose comment or
    // dispatch case overstates its degree fails here rather than in a run.
    //
    // Reference domains, matching the tables and the Lagrange functions:
    //   LINE          xi in [-1,1]
    //   QUAD, HEX     [-1,1]^d
    //   TRI           rows are barycentric ( 3 rows ); rows 0,1 are x,y of the
    //                 unit triangle x,y >= 0, x + y <= 1
    //   TET           rows are barycentric ( 4 rows ); rows 0..2 are x,y,z of
    //                 the unit tetrahedron
    //   PENTA         ( r, s ) on the unit triangle, t in [-1,1]
    //   PYRA          base [-1,1]^2 at zeta = 0, apex ( 0, 0, 1 )

    real
    factorial( const uint n )
    {
        real f = 1.0 ;
        for( uint k = 2; k <= n; ++k ) f *= k ;
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

    // integral of x^a y^b z^c over the pyramid with base [-1,1]^2 at z = 0 and
    // apex ( 0, 0, 1 ): the cross section at height z is [-(1-z),(1-z)]^2
    real
    pyra_moment( const uint a, const uint b, const uint c )
    {
        return line_moment( a ) * line_moment( b )
             * factorial( c ) * factorial( a + b + 2 ) / factorial( a + b + c + 3 );
    }

    real
    exact_moment( const GeometryType aGeometry, const uint a, const uint b, const uint c )
    {
        switch( aGeometry )
        {
            case( GeometryType::LINE )  : return line_moment( a );
            case( GeometryType::QUAD )  : return line_moment( a ) * line_moment( b );
            case( GeometryType::HEX )   : return line_moment( a ) * line_moment( b ) * line_moment( c );
            case( GeometryType::TRI )   : return tri_moment( a, b );
            case( GeometryType::TET )   : return tet_moment( a, b, c );
            case( GeometryType::PENTA ) : return tri_moment( a, b ) * line_moment( c );
            case( GeometryType::PYRA )  : return pyra_moment( a, b, c );
            default : return BELFEM_QUIET_NAN ;
        }
    }

    uint
    dimension_of( const GeometryType aGeometry )
    {
        switch( aGeometry )
        {
            case( GeometryType::LINE ) : return 1 ;
            case( GeometryType::QUAD ) :
            case( GeometryType::TRI )  : return 2 ;
            default : return 3 ;
        }
    }

    // every monomial of total degree <= aOrder must be integrated to
    // machine precision. The tolerance is 1e-12 absolute below |exact| = 1
    // and relative above it: the largest exact moment is 8 ( HEX, degree 0 ),
    // the roundoff of a 1331-point sum of cancelling O(1) terms is ~1e-13,
    // and an under-integrating table misses by 1e-8 or more, so the margin
    // is four decades on either side
    void
    check_moments( const GeometryType aGeometry, const uint aOrder,
                   const Vector< real > & tW, const Matrix< real > & tP )
    {
        const uint tDim = dimension_of( aGeometry );
        const uint tN   = tW.length();
        ASSERT_EQ( tP.n_cols(), tN );

        for( uint d = 0; d <= aOrder; ++d )
        {
            for( uint a = 0; a <= d; ++a )
            {
                for( uint b = 0; b <= d - a; ++b )
                {
                    const uint c = d - a - b ;
                    if( tDim < 3 && c > 0 ) continue ;
                    if( tDim < 2 && b > 0 ) continue ;

                    real tQ = 0.0 ;
                    for( uint k = 0; k < tN; ++k )
                    {
                        real tTerm = tW( k ) * std::pow( tP( 0, k ), a );
                        if( tDim > 1 ) tTerm *= std::pow( tP( 1, k ), b );
                        if( tDim > 2 ) tTerm *= std::pow( tP( 2, k ), c );
                        tQ += tTerm ;
                    }
                    const real tExact = exact_moment( aGeometry, a, b, c );
                    EXPECT_NEAR( tQ, tExact, 1e-12 * std::max( 1.0, std::abs( tExact ) ) )
                        << "geometry " << (int) aGeometry << " order " << aOrder
                        << " monomial x^" << a << " y^" << b << " z^" << c ;
                }
            }
        }
    }

    // the rule the dispatch returns for a requested order
    void
    check_exactness( const GeometryType aGeometry, const uint aOrder )
    {
        Vector< real > tW ;
        Matrix< real > tP ;
        intpoints( IntegrationScheme::GAUSS, aGeometry, aOrder, tW, tP );
        check_moments( aGeometry, aOrder, tW, tP );
    }

    void
    check_all_orders( const GeometryType aGeometry, const uint aMaxOrder )
    {
        for( uint p = 0; p <= aMaxOrder; ++p )
        {
            check_exactness( aGeometry, p );
        }
    }
}

// the maximum order per geometry is the highest case the dispatch in
// fn_intpoints.cpp provides; raising a table there means raising it here
TEST( IntegrationExactness, LINE  ) { check_all_orders( GeometryType::LINE,  20 ); }
TEST( IntegrationExactness, TRI   ) { check_all_orders( GeometryType::TRI,   20 ); }
TEST( IntegrationExactness, QUAD  ) { check_all_orders( GeometryType::QUAD,  20 ); }
TEST( IntegrationExactness, TET   ) { check_all_orders( GeometryType::TET,   14 ); }
TEST( IntegrationExactness, HEX   ) { check_all_orders( GeometryType::HEX,   20 ); }
TEST( IntegrationExactness, PENTA ) { check_all_orders( GeometryType::PENTA,  9 ); }
TEST( IntegrationExactness, PYRA  ) { check_all_orders( GeometryType::PYRA,   9 ); }

// the Shunn and Ham 10-point tetrahedron table is not dispatched, so the
// loop above never sees it; lock it directly at its degree ( 3 )
TEST( IntegrationExactness, TET10_table )
{
    Vector< real > tW ;
    Matrix< real > tP ;
    integration::gauss_tet10( tW, tP );
    check_moments( GeometryType::TET, 3, tW, tP );
}
