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
 * EF_HEX8 edge-direction signs under a flipped mesh edge. The circulation
 * of dof k along mesh edge e must be s_k on its own edge and zero on every
 * other; with every s_k = +1, as in the unflipped batteries of
 * test_EdgeFunctions.cpp, an exchange of two signs is invisible. Edge 7
 * ( top loop ) and edge 11 ( vertical ) are the pair the row names, whose
 * sign lookups share the same four-block layout in cl_EF_HEX8.cpp, so each
 * is flipped alone here: once through the circulation identity, and once against the
 * unflipped element, where the E, curl and gradient columns of the flipped
 * edge must change sign and no other column may.
 *
 * The corner and edge tables below are the ones test_EdgeFunctions.cpp
 * reads from cl_Element_HEX8.hpp.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "support/cl_EF_TestVolume.hpp"
#include "cl_EdgeFunctionFactory.hpp"

using namespace belfem;

namespace
{
    const real gEpsilon = 1.0e-11 ;

    // bottom loop, top loop, verticals
    const uint gHexEdges[ 12 ][ 2 ] = {
            { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
            { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
            { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };

    // reference corners of [-1,1]^3, EXODUS order
    const real gHexXi[ 8 ][ 3 ] = {
            { -1, -1, -1 }, { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 },
            { -1, -1,  1 }, { 1, -1,  1 }, { 1, 1,  1 }, { -1, 1,  1 } };

    // an axis-aligned box with three different edge lengths
    Matrix< real >
    box_corners()
    {
        const real tScale[ 3 ] = { 1.2, 0.9, 1.1 };
        Matrix< real > tX( 3, 8 );
        for ( uint n = 0; n < 8; ++n )
        {
            for ( uint d = 0; d < 3; ++d )
            {
                tX( d, n ) = tScale[ d ] * gHexXi[ n ][ d ];
            }
        }
        return tX ;
    }

    // circulation of every dof along every edge ( 2-point Gauss, straight
    // edges ), expected s_k delta_jk, with mesh edge aFlipEdge reversed
    void
    flipped_circulation_battery( const int aFlipEdge )
    {
        const Matrix< real > tX = box_corners();
        fem::test::EF_TestVolume tFixture( ElementType::HEX8, tX, aFlipEdge );

        const real tG = 1.0 / std::sqrt( 3.0 );
        Matrix< real > tXi( 3, 24 );
        for ( uint e = 0; e < 12; ++e )
        {
            const uint a = gHexEdges[ e ][ 0 ];
            const uint b = gHexEdges[ e ][ 1 ];
            for ( uint g = 0; g < 2; ++g )
            {
                const real t = ( g == 0 ) ? -tG : tG ;
                for ( uint d = 0; d < 3; ++d )
                {
                    tXi( d, 2 * e + g ) =
                          0.5 * ( gHexXi[ a ][ d ] + gHexXi[ b ][ d ] )
                        + 0.5 * t * ( gHexXi[ b ][ d ] - gHexXi[ a ][ d ] );
                }
            }
        }

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tEF = tFactory.create_edge_function( ElementType::HEX8 );
        tEF->precompute( tXi );
        tEF->link( tFixture.element() );

        real tS[ 12 ];
        tFixture.element()->edge_directions( tS );

        // the fixture really flipped the edge, and only that one; otherwise
        // the battery below degenerates into the unflipped one
        for ( uint k = 0; k < 12; ++k )
        {
            ASSERT_EQ( tS[ k ], ( int ) k == aFlipEdge ? -1.0 : 1.0 ) << "edge " << k ;
        }

        for ( uint e = 0; e < 12; ++e )
        {
            const uint a = gHexEdges[ e ][ 0 ];
            const uint b = gHexEdges[ e ][ 1 ];
            real tDx[ 3 ];
            for ( uint d = 0; d < 3; ++d )
            {
                tDx[ d ] = 0.5 * ( tX( d, b ) - tX( d, a ) );
            }
            for ( uint k = 0; k < 12; ++k )
            {
                real tCirc = 0.0 ;
                for ( uint g = 0; g < 2; ++g )
                {
                    const Matrix< real > & tE = tEF->E( 2 * e + g );
                    tCirc += tE( 0, k ) * tDx[ 0 ]
                           + tE( 1, k ) * tDx[ 1 ]
                           + tE( 2, k ) * tDx[ 2 ];
                }
                // along the CANONICAL edge direction: s_k on the own edge
                EXPECT_NEAR( tCirc, e == k ? tS[ k ] : 0.0, gEpsilon )
                        << "flip " << aFlipEdge << " edge " << e << " dof " << k ;
            }
        }
        delete tEF ;
    }

    // E, curl and gradient of a flipped element against the unflipped one
    // at the same interior points: column k changes sign exactly when k is
    // the flipped edge, and is unchanged otherwise. The circulation above
    // reads E alone; this also covers the derivative tables C() and G()
    // are built from
    void
    flipped_column_signs( const int aFlipEdge )
    {
        const Matrix< real > tX = box_corners();
        fem::test::EF_TestVolume tPlain( ElementType::HEX8, tX );
        fem::test::EF_TestVolume tFlipped( ElementType::HEX8, tX, aFlipEdge );

        // three pinned, asymmetric interior points
        Matrix< real > tXi( 3, 3 );
        const real tPoints[ 3 ][ 3 ] = {
                { 0.3, -0.2, 0.5 }, { -0.6, 0.4, -0.1 }, { 0.1, 0.7, -0.4 } };
        for ( uint p = 0; p < 3; ++p )
        {
            for ( uint d = 0; d < 3; ++d )
            {
                tXi( d, p ) = tPoints[ p ][ d ];
            }
        }

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tRef  = tFactory.create_edge_function( ElementType::HEX8 );
        fem::EdgeFunction * tFlip = tFactory.create_edge_function( ElementType::HEX8 );
        tRef->precompute( tXi );
        tFlip->precompute( tXi );
        tRef->link( tPlain.element() );
        tFlip->link( tFlipped.element() );

        for ( uint p = 0; p < 3; ++p )
        {
            const Matrix< real > tERef  = tRef->E( p );
            const Matrix< real > tEFlip = tFlip->E( p );
            const Matrix< real > tCRef  = tRef->C( p );
            const Matrix< real > tCFlip = tFlip->C( p );
            const Matrix< real > tGRef  = tRef->G( p );
            const Matrix< real > tGFlip = tFlip->G( p );

            for ( uint k = 0; k < 12; ++k )
            {
                const real tSign = ( int ) k == aFlipEdge ? -1.0 : 1.0 ;
                for ( uint i = 0; i < 3; ++i )
                {
                    EXPECT_NEAR( tEFlip( i, k ), tSign * tERef( i, k ), gEpsilon )
                        << "E flip " << aFlipEdge << " point " << p << " dof " << k ;
                    EXPECT_NEAR( tCFlip( i, k ), tSign * tCRef( i, k ), gEpsilon )
                        << "C flip " << aFlipEdge << " point " << p << " dof " << k ;
                }
                for ( uint i = 0; i < tGRef.n_rows(); ++i )
                {
                    EXPECT_NEAR( tGFlip( i, k ), tSign * tGRef( i, k ), gEpsilon )
                        << "G flip " << aFlipEdge << " point " << p << " dof " << k ;
                }
            }
        }
        delete tFlip ;
        delete tRef ;
    }
}

//------------------------------------------------------------------------------

TEST( EdgeFunctions, Hex8CirculationFlipEdge7 )
{
    flipped_circulation_battery( 7 );
}

//------------------------------------------------------------------------------

TEST( EdgeFunctions, Hex8CirculationFlipEdge11 )
{
    flipped_circulation_battery( 11 );
}

//------------------------------------------------------------------------------

TEST( EdgeFunctions, Hex8FlipEdge7NegatesItsOwnColumnOnly )
{
    flipped_column_signs( 7 );
}

//------------------------------------------------------------------------------

TEST( EdgeFunctions, Hex8FlipEdge11NegatesItsOwnColumnOnly )
{
    flipped_column_signs( 11 );
}
