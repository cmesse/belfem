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
 * The node-to-control-point weights of the cubic tensor meshes ( QUAD16,
 * HEX64 ), built end to end through Mesh. Each node is a weighted sum of the
 * uniform cubic B-splines of its elements, the Cox-de Boor recursion
 * ( Zienkiewicz et al. 2013, Eq. 17.4 ) evaluated for uniform knots. The
 * control values that make the spline equal x,
 * x^2 and x^3 are the polar forms of those monomials over the knots g-h, g,
 * g+h: g, g^2 - h^2/3 and g^3 - g h^2, where g = x0 + ( a - 1 ) h is the
 * Greville abscissa of control point a ( a = 0 ... ne+2 per direction ) and
 * h the element size. With them the weights must reproduce the monomials
 * exactly at every node. The test checks the weights only; it computes g
 * from the lattice index and does not read the stored ControlPoint
 * coordinates.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_ControlPoint.hpp"

using namespace belfem;

namespace
{
    const real tTol = 1e-11;

    // B-spline control values that reproduce x^p for a uniform cubic
    real
    control_value( const int aPower, const real aG, const real aH )
    {
        switch ( aPower )
        {
            case 1 : return aG ;
            case 2 : return aG * aG - aH * aH / 3.0 ;
            default: return aG * aG * aG - aG * aH * aH ;
        }
    }

    // checks partition of unity, control-point coverage, and exact
    // reproduction of x^p, p = 1..3, along every direction
    void
    check_cubic_mesh( Mesh & aMesh,
                      const Vector< index_t > & aNumNodes,
                      const Vector< real > & aStep,
                      const Vector< real > & aOrigin )
    {
        const uint tDim = aNumNodes.length();

        Cell< index_t > tNumCP( tDim, 0 );
        for ( uint d = 0; d < tDim; ++d )
        {
            tNumCP( d ) = ( aNumNodes( d ) - 1 ) / 3 + 3 ;
        }

        index_t tExpected = 1 ;
        for ( uint d = 0; d < tDim; ++d )
        {
            tExpected *= tNumCP( d );
        }

        const index_t tNumControlPoints = aMesh.number_of_control_points();
        ASSERT_EQ( tNumControlPoints, tExpected );

        Cell< bool > tUsed( tNumControlPoints, false );

        for ( mesh::Node * tNode : aMesh.nodes() )
        {
            const uint tN = tNode->number_of_sources();
            ASSERT_GT( tN, 0u ) << "node " << tNode->id() ;

            real tSum = 0.0 ;
            for ( uint s = 0; s < tN; ++s )
            {
                tSum += tNode->weight( s );
            }
            EXPECT_NEAR( tSum, 1.0, tTol ) << "node " << tNode->id() ;

            for ( uint d = 0; d < tDim; ++d )
            {
                for ( int p = 1; p <= 3; ++p )
                {
                    real tValue = 0.0 ;
                    for ( uint s = 0; s < tN; ++s )
                    {
                        // control point ids, last index fastest: id - 1 = ( i nj + j ) nk + k
                        index_t tIndex = tNode->source( s )->id() - 1 ;
                        tUsed( tIndex ) = true ;

                        index_t tLattice = tIndex ;
                        for ( uint e = tDim - 1; e > d; --e )
                        {
                            tLattice /= tNumCP( e );
                        }
                        tLattice %= tNumCP( d );

                        const real tH = 3.0 * aStep( d ) ;
                        const real tG = aOrigin( d ) + ( ( real ) tLattice - 1.0 ) * tH ;
                        tValue += tNode->weight( s ) * control_value( p, tG, tH );
                    }

                    const real tX = tNode->x( d );
                    EXPECT_NEAR( tValue, std::pow( tX, p ), tTol * std::max( 1.0, std::pow( std::abs( tX ), p ) ) )
                        << "node " << tNode->id() << ", direction " << d << ", power " << p ;
                }
            }
        }

        for ( index_t c = 0; c < tNumControlPoints; ++c )
        {
            EXPECT_TRUE( tUsed( c ) ) << "control point " << c + 1 << " has no node" ;
        }
    }
}

//------------------------------------------------------------------------------

TEST( TensorMeshBSpline, Quad16ReproducesCubics )
{
    // 3 x 3 elements, 10 x 10 nodes
    const Vector< index_t > tNumNodes( { 10, 10 } );
    const Vector< real > tStep( { 0.5, 0.5 } );
    const Vector< real > tOrigin( { 0.0, 0.0 } );
    Mesh tMesh( 3, tNumNodes, tStep, tOrigin );
    check_cubic_mesh( tMesh, tNumNodes, tStep, tOrigin );
}

//------------------------------------------------------------------------------

TEST( TensorMeshBSpline, Quad16TranslatedAnisotropic )
{
    // 3 x 2 elements, unequal steps, origin off zero
    const Vector< index_t > tNumNodes( { 10, 7 } );
    const Vector< real > tStep( { 0.25, 0.7 } );
    const Vector< real > tOrigin( { -1.0, 2.5 } );
    Mesh tMesh( 3, tNumNodes, tStep, tOrigin );
    check_cubic_mesh( tMesh, tNumNodes, tStep, tOrigin );
}

//------------------------------------------------------------------------------

TEST( TensorMeshBSpline, Hex64ReproducesCubics )
{
    // 2 x 2 x 2 elements, 7 x 7 x 7 nodes
    const Vector< index_t > tNumNodes( { 7, 7, 7 } );
    const Vector< real > tStep( { 0.5, 0.4, 0.3 } );
    const Vector< real > tOrigin( { 0.0, 1.0, -2.0 } );
    Mesh tMesh( 3, tNumNodes, tStep, tOrigin );
    check_cubic_mesh( tMesh, tNumNodes, tStep, tOrigin );
}
