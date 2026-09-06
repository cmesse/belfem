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
 * Interface / orientation regression battery (sign, orientation,
 * master-slave, indexing bug class). See todo/falsification_tooling.md.
 *
 * Physics-load-bearing expected values are derived in comments and flagged
 * EXPECTED: pending Christian sign-off. Reviewer agreement does not settle
 * physics; the flag is replaced by a provenance comment on sign-off.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "support/cl_TS_TestStack.hpp"
#include "mt_maxwell_h.hpp"
#include "cl_TimestepMatrices.hpp"
#include "nedelec/cl_EF_QUAD4TS.hpp"
#include "nedelec/cl_EF_PENTA6TS.hpp"
#include "nedelec/cl_EF_HEX8TS.hpp"
#include "nedelec/cl_EF_HEX8TB.hpp"
#include "meshtools.hpp"

using namespace belfem;

namespace
{
    const real tEpsilon = 1.0e-12 ;

    /**
     * Evaluation points on the QUAD4TS reference element:
     * columns 0,1 : 2-point Gauss along the bottom edge ( eta = -1 )
     * columns 2,3 : 2-point Gauss along the top edge    ( eta = +1 )
     * column  4   : element center
     */
    Matrix< real >
    quad4ts_eval_points()
    {
        const real tG = 1.0 / std::sqrt( 3.0 );
        Matrix< real > aXi( 2, 5 );
        aXi( 0, 0 ) = -tG ; aXi( 1, 0 ) = -1.0 ;
        aXi( 0, 1 ) =  tG ; aXi( 1, 1 ) = -1.0 ;
        aXi( 0, 2 ) = -tG ; aXi( 1, 2 ) =  1.0 ;
        aXi( 0, 3 ) =  tG ; aXi( 1, 3 ) =  1.0 ;
        aXi( 0, 4 ) =  0.0 ; aXi( 1, 4 ) =  0.0 ;
        return aXi ;
    }

    /**
     * circulation of E-column aDof along the edge at rows aFirstPoint,
     * aFirstPoint+1 of quad4ts_eval_points(): 2-point Gauss with the
     * physical tangent, sum_g w_g ( E_k . t ) * ( L / 2 ), w_g = 1.
     */
    real
    edge_circulation(
            fem::EF_QUAD4TS  & aEF,
            const uint         aDof,
            const uint         aFirstPoint,
            const real         aTx,
            const real         aTy,
            const real         aLength )
    {
        real tCirc = 0.0 ;
        for ( uint g = 0; g < 2; ++g )
        {
            const Matrix< real > & tE = aEF.E( aFirstPoint + g );
            tCirc += ( tE( 0, aDof ) * aTx + tE( 1, aDof ) * aTy )
                     * 0.5 * aLength ;
        }
        return tCirc ;
    }
}

//------------------------------------------------------------------------------
// QUAD4TS: unit circulation per edge dof
//------------------------------------------------------------------------------

// Derivation: E_k = s_k * f_k( eta ) * nabla_xi with f_0 = (1-eta)/2,
// f_1 = (1+eta)/2 and nabla_xi = t_hat / L. Along an edge the tangential
// trace is s_k * f_k / L, so the edge circulation is s_k * f_k( eta_edge ):
// bottom edge -> ( s_0, 0 ), top edge -> ( 0, s_1 ). Both edges circulate
// along +tape (the 2026-08-04 convention that keeps shared inter-layer
// edges continuous, commit 4a42d982).
// EXPECTED: pending Christian sign-off
TEST( InterfaceOrientation, Quad4TsUnitCirculation )
{
    for ( real tAngle : { 0.0, 0.7 } )
    {
        fem::test::TS_TestStack2D tStack( 1, 2.0, 0.1, tAngle );

        fem::EF_QUAD4TS tEF ;
        Matrix< real > tXi = quad4ts_eval_points();
        tEF.precompute( tXi );
        tEF.link( tStack.element( 0 ) );

        const real tTx = std::cos( tAngle );
        const real tTy = std::sin( tAngle );
        const real tL  = tStack.length();

        // bottom edge: dof 0 carries unit circulation, dof 1 nothing
        EXPECT_NEAR( edge_circulation( tEF, 0, 0, tTx, tTy, tL ),  1.0, tEpsilon );
        EXPECT_NEAR( edge_circulation( tEF, 1, 0, tTx, tTy, tL ),  0.0, tEpsilon );

        // top edge: dof 1 carries unit circulation ALONG +TAPE, dof 0 nothing.
        // The historical defect ( -mS[1], pre-4a42d982 ) returns -1 here.
        EXPECT_NEAR( edge_circulation( tEF, 0, 2, tTx, tTy, tL ),  0.0, tEpsilon );
        EXPECT_NEAR( edge_circulation( tEF, 1, 2, tTx, tTy, tL ),  1.0, tEpsilon );
    }
}

//------------------------------------------------------------------------------
// QUAD4TS: edge-direction sign follows the mesh edge orientation
//------------------------------------------------------------------------------

// Structural: flipping the stored node order of the shared mesh edge must
// flip the sign of that dof's circulation and nothing else.
TEST( InterfaceOrientation, Quad4TsEdgeDirectionSign )
{
    // flip row 1 ( = top edge of the single layer )
    fem::test::TS_TestStack2D tStack( 1, 2.0, 0.1, 0.0, 1 );

    fem::EF_QUAD4TS tEF ;
    Matrix< real > tXi = quad4ts_eval_points();
    tEF.precompute( tXi );
    tEF.link( tStack.element( 0 ) );

    EXPECT_NEAR( edge_circulation( tEF, 0, 0, 1.0, 0.0, 2.0 ),  1.0, tEpsilon );
    EXPECT_NEAR( edge_circulation( tEF, 1, 2, 1.0, 0.0, 2.0 ), -1.0, tEpsilon );
}

//------------------------------------------------------------------------------
// QUAD4TS: E-column face activity
//------------------------------------------------------------------------------

// Structural: the bottom dof column vanishes on the top face and vice versa
// ( f_0(+1) = f_1(-1) = 0 ).
TEST( InterfaceOrientation, Quad4TsFaceActivity )
{
    fem::test::TS_TestStack2D tStack( 1 );

    fem::EF_QUAD4TS tEF ;
    Matrix< real > tXi = quad4ts_eval_points();
    tEF.precompute( tXi );
    tEF.link( tStack.element( 0 ) );

    for ( uint g = 0; g < 2; ++g )
    {
        // top edge points: bottom column zero
        const Matrix< real > & tTop = tEF.E( 2 + g );
        EXPECT_NEAR( tTop( 0, 0 ), 0.0, tEpsilon );
        EXPECT_NEAR( tTop( 1, 0 ), 0.0, tEpsilon );
    }
    for ( uint g = 0; g < 2; ++g )
    {
        // bottom edge points: top column zero
        const Matrix< real > & tBot = tEF.E( g );
        EXPECT_NEAR( tBot( 0, 1 ), 0.0, tEpsilon );
        EXPECT_NEAR( tBot( 1, 1 ), 0.0, tEpsilon );
    }
}

//------------------------------------------------------------------------------
// QUAD4TS: Stokes consistency of E and C
//------------------------------------------------------------------------------

// Structural, coordinate-free: for every dof k the boundary circulation
// must equal the curl integral, oint_dOmega E_k . dl = C_k * Area. E has no
// through-thickness component, so the boundary reduces to bottom ( +tape )
// minus top ( +tape ): s_k * ( f_k(-1) - f_k(+1) ). A one-sided sign flip
// in either E or C breaks this identity. Catches the greg2 class even when
// only one operator is wrong.
TEST( InterfaceOrientation, Quad4TsStokesConsistency )
{
    for ( real tAngle : { 0.0, 0.7 } )
    {
        fem::test::TS_TestStack2D tStack( 1, 2.0, 0.1, tAngle );

        fem::EF_QUAD4TS tEF ;
        Matrix< real > tXi = quad4ts_eval_points();
        tEF.precompute( tXi );
        tEF.link( tStack.element( 0 ) );

        const real tTx   = std::cos( tAngle );
        const real tTy   = std::sin( tAngle );
        const real tL    = tStack.length();
        const real tArea = tStack.length() * tStack.thickness();

        // quadrature area must agree with the geometry
        EXPECT_NEAR( tEF.abs_det_J() * tEF.sum_w(), tArea, tEpsilon );

        for ( uint k = 0; k < 2; ++k )
        {
            const real tBoundary =
                    edge_circulation( tEF, k, 0, tTx, tTy, tL )
                  - edge_circulation( tEF, k, 2, tTx, tTy, tL );

            // C is constant on the element; evaluate at the center point
            const real tCurlArea = tEF.C( 4 )( 0, k ) * tArea ;

            EXPECT_NEAR( tBoundary, tCurlArea, tEpsilon );
        }
    }
}

//------------------------------------------------------------------------------
// QUAD4TS: inter-layer tangential continuity ( the greg2 mechanism )
//------------------------------------------------------------------------------

// Structural: layers l and l+1 share the interface edge ( same mesh Edge
// object ). With a unit dof on that shared edge and zero elsewhere, the
// tangential field evaluated from ABOVE ( layer l+1 at eta=-1 ) must equal
// the field evaluated from BELOW ( layer l at eta=+1 ). The pre-4a42d982
// -mS[1] convention flips the from-below trace, producing the alternating
// J/Jc Gregory observed.
TEST( InterfaceOrientation, Quad4TsInterlayerContinuity )
{
    for ( real tAngle : { 0.0, 0.7 } )
    {
        fem::test::TS_TestStack2D tStack( 2, 2.0, 0.1, tAngle );

        Matrix< real > tXi = quad4ts_eval_points();

        fem::EF_QUAD4TS tLower ;
        tLower.precompute( tXi );
        tLower.link( tStack.element( 0 ) );

        fem::EF_QUAD4TS tUpper ;
        tUpper.precompute( tXi );
        tUpper.link( tStack.element( 1 ) );

        // unit dof on the shared edge: dof 1 of the lower layer,
        // dof 0 of the upper layer
        for ( uint g = 0; g < 2; ++g )
        {
            const Matrix< real > & tFromBelow = tLower.E( 2 + g ); // eta=+1
            const real tBx = tFromBelow( 0, 1 );
            const real tBy = tFromBelow( 1, 1 );

            const Matrix< real > & tFromAbove = tUpper.E( g );     // eta=-1
            const real tAx = tFromAbove( 0, 0 );
            const real tAy = tFromAbove( 1, 0 );

            EXPECT_NEAR( tBx, tAx, tEpsilon );
            EXPECT_NEAR( tBy, tAy, tEpsilon );
        }
    }
}

//------------------------------------------------------------------------------
// QUAD4TS: Ampere telescope over an N-layer stack
//------------------------------------------------------------------------------

// Derivation: the sheet current of layer l is I_l = int C.q dA. C is
// constant, the area is t*L, and with all edge directions positive
// C = ( +1, -1 ) / ( t*L ), so I_l = h_l - h_{l+1} where h_r is the
// tangential trace dof on interface row r ( circulation convention along
// +tape, curl orientation ( nabla_eta x nabla_xi )_z ). The layer currents
// must telescope to the outer-trace difference, sum_l I_l = h_0 - h_N
// ( Ampere; dl20260804_2d_thinshell_layer_alternation.md:107 ), and a
// uniform tangential field h_r = const must carry zero current in EVERY
// layer — the historical defect made the per-layer currents alternate.
// This is the minutes-scale test that would have caught greg2.
// EXPECTED: pending Christian sign-off
TEST( InterfaceOrientation, Quad4TsAmpereTelescope )
{
    const uint tN = 4 ;                    // greg2 configuration: 4 layers

    for ( real tAngle : { 0.0, 0.7 } )
    {
        fem::test::TS_TestStack2D tStack( tN, 2.0, 0.1, tAngle );
        const real tArea = tStack.length() * tStack.thickness();

        Matrix< real > tXi = quad4ts_eval_points();

        // interface traces: uniform and linear-through-thickness
        Vector< real > tUniform( tN + 1, 3.0 );
        Vector< real > tLinear( tN + 1 );
        for ( uint r = 0; r <= tN; ++r )
        {
            tLinear( r ) = 1.0 + 0.5 * r ;
        }

        real tSumUniform = 0.0 ;
        real tSumLinear  = 0.0 ;

        for ( uint l = 0; l < tN; ++l )
        {
            fem::EF_QUAD4TS tEF ;
            tEF.precompute( tXi );
            tEF.link( tStack.element( l ) );

            const Matrix< real > & tC = tEF.C( 4 );

            // I_l = ( C . q ) * Area ( C constant over the element )
            const real tIUniform =
                    ( tC( 0, 0 ) * tUniform( l )
                    + tC( 0, 1 ) * tUniform( l + 1 ) ) * tArea ;
            const real tILinear =
                    ( tC( 0, 0 ) * tLinear( l )
                    + tC( 0, 1 ) * tLinear( l + 1 ) ) * tArea ;

            // uniform tangential field: zero current in every layer
            // ( alternation detector — the defect produced +/- 2h here )
            EXPECT_NEAR( tIUniform, 0.0, tEpsilon );

            // per-layer Ampere: I_l = h_l - h_{l+1}
            // EXPECTED: pending Christian sign-off
            EXPECT_NEAR( tILinear, tLinear( l ) - tLinear( l + 1 ), tEpsilon );

            tSumUniform += tIUniform ;
            tSumLinear  += tILinear ;
        }

        // telescope: total current = outer-trace difference
        // EXPECTED: pending Christian sign-off
        EXPECT_NEAR( tSumUniform, 0.0, tEpsilon );
        EXPECT_NEAR( tSumLinear, tLinear( 0 ) - tLinear( tN ), tEpsilon );
    }
}

//------------------------------------------------------------------------------
// 3-D thin shells: circulation matrix on the reference prism
//------------------------------------------------------------------------------

namespace
{
    // parametric corners of the thin-shell prisms ( xi, eta, zeta );
    // node ordering follows the Lagrange convention ( shape k <-> node k ),
    // bottom surface at zeta=-1, top at zeta=+1
    void
    prism_corner( const ElementType aType, const uint aNode, real * aXi )
    {
        if ( aType == ElementType::PENTA6TS )
        {
            const real tXi[ 6 ]  = { 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
            const real tEta[ 6 ] = { 0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
            aXi[ 0 ] = tXi[ aNode ];
            aXi[ 1 ] = tEta[ aNode ];
            aXi[ 2 ] = aNode < 3 ? -1.0 : 1.0 ;
        }
        else // HEX8TS
        {
            const real tXi[ 8 ]  = { -1.0,  1.0, 1.0, -1.0, -1.0,  1.0, 1.0, -1.0 };
            const real tEta[ 8 ] = { -1.0, -1.0, 1.0,  1.0, -1.0, -1.0, 1.0,  1.0 };
            aXi[ 0 ] = tXi[ aNode ];
            aXi[ 1 ] = tEta[ aNode ];
            aXi[ 2 ] = aNode < 4 ? -1.0 : 1.0 ;
        }
    }

    // index of a node within its element
    uint
    node_index( mesh::Element * aElement, const mesh::Node * aNode )
    {
        for ( uint k = 0; k < aElement->number_of_nodes(); ++k )
        {
            if ( aElement->node( k ) == aNode ) return k ;
        }
        ADD_FAILURE() << "node not on element" ;
        return 0 ;
    }

    /**
     * circulation matrix tCirc( j, k ) = int_edge_j E_k . dl for a 3-D
     * thin-shell edge function, 2-point Gauss along the canonical
     * ( get_nodes_of_edge ) direction of each straight edge.
     */
    template < typename EF >
    void
    prism_circulation_matrix(
            const ElementType aType,
            fem::Element * aFemElement,
            Matrix< real > & aCirc )
    {
        mesh::Element * tElement = aFemElement->element();
        const uint tNumEdges = tElement->number_of_edges();

        const real tG = 1.0 / std::sqrt( 3.0 );
        Cell< mesh::Node * > tEdgeNodes ;

        aCirc.set_size( tNumEdges, tNumEdges, 0.0 );

        for ( uint j = 0; j < tNumEdges; ++j )
        {
            tElement->get_nodes_of_edge( j, tEdgeNodes );
            const uint tA = node_index( tElement, tEdgeNodes( 0 ) );
            const uint tB = node_index( tElement, tEdgeNodes( 1 ) );

            real tXiA[ 3 ]; real tXiB[ 3 ];
            prism_corner( aType, tA, tXiA );
            prism_corner( aType, tB, tXiB );

            // two Gauss points along the edge, parametric
            Matrix< real > tXi( 3, 2 );
            for ( uint i = 0; i < 3; ++i )
            {
                tXi( i, 0 ) = 0.5 * ( tXiA[ i ] + tXiB[ i ] )
                            - 0.5 * tG * ( tXiB[ i ] - tXiA[ i ] );
                tXi( i, 1 ) = 0.5 * ( tXiA[ i ] + tXiB[ i ] )
                            + 0.5 * tG * ( tXiB[ i ] - tXiA[ i ] );
            }

            // physical half-edge vector ( straight edges )
            const real tDx = 0.5 * ( tEdgeNodes( 1 )->x() - tEdgeNodes( 0 )->x() );
            const real tDy = 0.5 * ( tEdgeNodes( 1 )->y() - tEdgeNodes( 0 )->y() );
            const real tDz = 0.5 * ( tEdgeNodes( 1 )->z() - tEdgeNodes( 0 )->z() );

            EF tEF ;
            tEF.precompute( tXi );
            tEF.link( aFemElement );

            for ( uint g = 0; g < 2; ++g )
            {
                const Matrix< real > & tE = tEF.E( g );
                for ( uint k = 0; k < tNumEdges; ++k )
                {
                    aCirc( j, k ) += tE( 0, k ) * tDx
                                   + tE( 1, k ) * tDy
                                   + tE( 2, k ) * tDz ;
                }
            }
        }
    }
}

// Derivation: same identity as 2-D per shell surface — each edge dof k
// carries unit circulation along its own canonical edge and none along the
// others, int_edge_j E_k . dl = s_k delta_jk ( with all fixture edges
// canonical, s = +1 ). Dof column k is expected to correspond to element
// edge k.
// EXPECTED: pending Christian sign-off
TEST( InterfaceOrientation, Penta6TsUnitCirculation )
{
    fem::test::TS_TestPrism tPrism( ElementType::PENTA6TS );

    Matrix< real > tCirc ;
    prism_circulation_matrix< fem::EF_PENTA6TS >(
            ElementType::PENTA6TS, tPrism.element(), tCirc );

    for ( uint j = 0; j < 6; ++j )
    {
        for ( uint k = 0; k < 6; ++k )
        {
            EXPECT_NEAR( tCirc( j, k ), j == k ? 1.0 : 0.0, tEpsilon )
                << "edge " << j << " dof " << k ;
        }
    }
}

// EXPECTED: pending Christian sign-off ( same identity, HEX8TS, 8 dofs )
TEST( InterfaceOrientation, Hex8TsUnitCirculation )
{
    fem::test::TS_TestPrism tPrism( ElementType::HEX8TS );

    Matrix< real > tCirc ;
    prism_circulation_matrix< fem::EF_HEX8TS >(
            ElementType::HEX8TS, tPrism.element(), tCirc );

    for ( uint j = 0; j < 8; ++j )
    {
        for ( uint k = 0; k < 8; ++k )
        {
            EXPECT_NEAR( tCirc( j, k ), j == k ? 1.0 : 0.0, tEpsilon )
                << "edge " << j << " dof " << k ;
        }
    }
}

// quadrature volume must agree with the geometry: sum_w() is the reference
// wedge measure ( 1 ), abs_det_J() the physical volume
TEST( InterfaceOrientation, Penta6TsVolumeIdentity )
{
    const real tThickness = 0.1 ;
    fem::test::TS_TestPrism tPrism( ElementType::PENTA6TS, tThickness );

    Matrix< real > tXi( 3, 1, 0.0 );
    tXi( 0, 0 ) = 1.0 / 3.0 ;
    tXi( 1, 0 ) = 1.0 / 3.0 ;

    fem::EF_PENTA6TS tEF ;
    tEF.precompute( tXi );
    tEF.link( tPrism.element() );

    // triangle area from the facet's three nodes
    const mesh::Facet * tFacet = tPrism.element()->facet();
    real tA[ 3 ]; real tB[ 3 ];
    for ( uint i = 0; i < 3; ++i )
    {
        tA[ i ] = tFacet->node( 1 )->x( i ) - tFacet->node( 0 )->x( i );
        tB[ i ] = tFacet->node( 2 )->x( i ) - tFacet->node( 0 )->x( i );
    }
    const real tNx = tA[ 1 ] * tB[ 2 ] - tA[ 2 ] * tB[ 1 ];
    const real tNy = tA[ 2 ] * tB[ 0 ] - tA[ 0 ] * tB[ 2 ];
    const real tNz = tA[ 0 ] * tB[ 1 ] - tA[ 1 ] * tB[ 0 ];
    const real tArea = 0.5 * std::sqrt( tNx * tNx + tNy * tNy + tNz * tNz );

    EXPECT_NEAR( tEF.sum_w(), 1.0, tEpsilon );
    EXPECT_NEAR( tEF.abs_det_J() * tEF.sum_w(), tArea * tThickness, tEpsilon );
}

// Structural: one flipped mesh edge flips exactly its own diagonal entry.
TEST( InterfaceOrientation, Penta6TsEdgeFlipSign )
{
    fem::test::TS_TestPrism tPrism( ElementType::PENTA6TS, 0.1, 4 );

    Matrix< real > tCirc ;
    prism_circulation_matrix< fem::EF_PENTA6TS >(
            ElementType::PENTA6TS, tPrism.element(), tCirc );

    for ( uint j = 0; j < 6; ++j )
    {
        const real tExpect = j == 4 ? -1.0 : 1.0 ;
        EXPECT_NEAR( tCirc( j, j ), tExpect, tEpsilon ) << "edge " << j ;
    }
}

//------------------------------------------------------------------------------
// HEX8TB ( side-connector wall, current exact-cuboid implementation )
//------------------------------------------------------------------------------

// Derivation: the wall's four dofs sit on the four longitudinal edges
// ( canonical order along +xi: {0,1}, {3,2}, {4,5}, {7,6} ). The h_t
// stream-function design assigns each edge dof a unit line integral along
// its own edge. Scope per O4: the CURRENT exact-cuboid EF only ( j_t == 0 );
// wall-term physics stays with hex8tb_phase2_fem_wiring.md.
// RE-ENABLED 2026-08-10: the blocker is gone — the ElementType::HEX8TB case
// landed in cl_IF_InterpolationFunctionFactory.cpp with hex8tb R4, so the
// kernel chain no longer dies with "no lagrange function available"
// Dimensions check out: the fixture's element is the
// HEX8TB wall, number_of_edges() == 4 ( meshtools.cpp ), which matches
// EF_HEX8TB's four columns and the 4x4 block asserted below; the HEX8TS
// argument selects the 8-node corner table only.
// EXPECTED: pending Christian sign-off ( same caveat as the enabled
// Penta6Ts/Hex8Ts siblings — first run of this one is the sign-off )
TEST( InterfaceOrientation, Hex8TbUnitCirculation )
{
    fem::test::TS_TestWall tWall ;

    // the wall's corner layout equals the standard hex table used for
    // HEX8TS ( nodes 0-3 bottom, 4-7 top )
    Matrix< real > tCirc ;
    prism_circulation_matrix< fem::EF_HEX8TB >(
            ElementType::HEX8TS, tWall.element(), tCirc );

    for ( uint j = 0; j < 4; ++j )
    {
        for ( uint k = 0; k < 4; ++k )
        {
            EXPECT_NEAR( tCirc( j, k ), j == k ? 1.0 : 0.0, tEpsilon )
                << "edge " << j << " dof " << k ;
        }
    }
}

//------------------------------------------------------------------------------
// get_top_nodes / get_bottom_nodes positional alignment ( 43474c9f class )
//------------------------------------------------------------------------------

// Structural, geometric: get_top_nodes( k ) must be the node vertically
// above get_bottom_nodes( k ) — the thin-shell node-tie contract. The
// 2026-08-04 top-tie defect ( fixed in 43474c9f ) returned the top nodes in
// facet order instead of master order, i.e. positionally swapped.
TEST( InterfaceOrientation, TopBottomNodeAlignment )
{
    // 2-D stack ( QUAD4TS, rotated so the through-thickness direction is
    // not axis-aligned )
    {
        fem::test::TS_TestStack2D tStack( 1, 2.0, 0.1, 0.7 );
        mesh::Element * tElement = tStack.element( 0 )->element();

        Cell< mesh::Node * > tBottom ;
        Cell< mesh::Node * > tTop ;
        mesh::get_bottom_nodes( tElement, tBottom );
        mesh::get_top_nodes( tElement, tTop );

        ASSERT_EQ( tBottom.size(), tTop.size() );

        const real tNx = -std::sin( 0.7 ) * 0.1 ;
        const real tNy =  std::cos( 0.7 ) * 0.1 ;
        for ( uint k = 0; k < tBottom.size(); ++k )
        {
            EXPECT_NEAR( tTop( k )->x() - tBottom( k )->x(), tNx, tEpsilon )
                << "position " << k ;
            EXPECT_NEAR( tTop( k )->y() - tBottom( k )->y(), tNy, tEpsilon )
                << "position " << k ;
        }
    }

    // 3-D prisms ( extruded in +z by the thickness )
    for ( ElementType tType : { ElementType::PENTA6TS, ElementType::HEX8TS } )
    {
        fem::test::TS_TestPrism tPrism( tType, 0.1 );
        mesh::Element * tElement = tPrism.element()->element();

        Cell< mesh::Node * > tBottom ;
        Cell< mesh::Node * > tTop ;
        mesh::get_bottom_nodes( tElement, tBottom );
        mesh::get_top_nodes( tElement, tTop );

        ASSERT_EQ( tBottom.size(), tTop.size() );

        for ( uint k = 0; k < tBottom.size(); ++k )
        {
            EXPECT_NEAR( tTop( k )->x() - tBottom( k )->x(), 0.0, tEpsilon );
            EXPECT_NEAR( tTop( k )->y() - tBottom( k )->y(), 0.0, tEpsilon );
            EXPECT_NEAR( tTop( k )->z() - tBottom( k )->z(), 0.1, tEpsilon );
        }
    }
}

//------------------------------------------------------------------------------
// Ghost facet element contract ( scaffold )
//------------------------------------------------------------------------------

// The h_ghost interior-penalty kernel ( mt_maxwell_h.cpp, h_ghost();
// src/fem/maxwell/doc/ghost_penalty_stabilization.md ) couples master and
// slave layer edge dofs through four blocks K++, K+-, K-+, K-- ( 2 x 6 = 12
// dofs for PENTA6TS layers ). TS_TestGhostStack stands up the full
// production Kernel/DofManager/IWG_Maxwell chain over a hand-built
// two-layer ghost mesh, so the contract runs against the real Calculator:
//
//   1. 12 dofs, and the kernel's silent 2m == n assumption ( n nedelec
//      dofs on the layer element, m on the facet )
//   2. interface locality: dofs away from the interface ( master bottom,
//      slave top ) couple only through the rho_harm-scaled gradient
//      terms — their columns scale linearly in rho while the
//      k_reg-dominated penalty columns stay put
//   3. annihilation: a constant tangential field continuous across the
//      interface lies in the kernel of K; a master-only half field does
//      not ( negative control )
//   4. symmetry: the adjoint-consistency terms make K algebraically
//      symmetric ( symmetric weighted interior penalty; ruled intentional
//      2026-08-28, and ghost_penalty_stabilization.md now says so )
namespace
{
    // dof layout in K: master element edges 0-5, then slave edges 0-5;
    // local edges 0-2 are the bottom row, 3-5 the top row, so the
    // interface carries master {3,4,5} and slave {6,7,8}
    const uint tGhostFarDofs[ 6 ]  = { 0, 1, 2, 9, 10, 11 };
    const uint tGhostNearDofs[ 6 ] = { 3, 4, 5, 6, 7, 8 };

    // L1 norm of column j
    real
    ghost_column_norm( const Matrix< real > & aK, const uint aJ )
    {
        real tNorm = 0.0 ;
        for ( uint i = 0; i < 12; ++i )
        {
            tNorm += std::abs( aK( i, aJ ) );
        }
        return tNorm ;
    }

    // lowest-order edge dof of the constant field aC: line integral along
    // the edge in the edge object's own node order
    real
    ghost_edge_dof( mesh::Edge * aEdge, const real aC[ 3 ] )
    {
        return aC[ 0 ] * ( aEdge->node( 1 )->x() - aEdge->node( 0 )->x() )
             + aC[ 1 ] * ( aEdge->node( 1 )->y() - aEdge->node( 0 )->y() )
             + aC[ 2 ] * ( aEdge->node( 1 )->z() - aEdge->node( 0 )->z() );
    }

    real
    ghost_max_Kq( const Matrix< real > & aK, const Vector< real > & aQ )
    {
        real aMax = 0.0 ;
        for ( uint i = 0; i < 12; ++i )
        {
            real tSum = 0.0 ;
            for ( uint j = 0; j < 12; ++j )
            {
                tSum += aK( i, j ) * aQ( j );
            }
            aMax = std::max( aMax, std::abs( tSum ) );
        }
        return aMax ;
    }
}

TEST( InterfaceOrientation, GhostElementContract )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // canonical edge orientation, plus one legal flipped column ( local
    // edge 1 reversed in every row of both prisms — flipping fewer rows
    // is an illegal mesh state, see the fixture comment )
    for ( int tFlip : { -1, 1 } )
    {
        fem::test::TS_TestGhostStack tStack( 0.1, tFlip );
        fem::Calculator * tCalc = tStack.calculator();
        tCalc->link( tStack.ghost_element() );

        // contract: link_to_group ran in the fixture on a kernel
        // that has no Controller — it must leave the thermal flag down
        // rather than read an uninitialized pointer
        ASSERT_FALSE( tStack.iwg()->have_thermal() )
                << "controller-less link_to_group must not report a thermal kernel" ;

        // ---- 1. the 12-dof contract, and the 2m == n assumption ----
        const uint n = mesh::number_of_nedelec_dofs( ElementType::PENTA6TS );
        const uint m = mesh::number_of_nedelec_dofs( ElementType::TRI3 );
        ASSERT_EQ( n, 6u );
        ASSERT_EQ( 2 * m, n ) << "h_ghost's Dm/Ds assembly assumes 2m == n" ;

        fem::TimestepMatrices tMat ;
        tMat.set_flag( fem::MatrixFlag::K );
        tMat.initialize( 2 * n );
        fem::maxwell::h_ghost( tCalc, & tMat );

        // copy: compared below against the 10x-rho assembly
        Matrix< real > tK = tMat.K();

        ASSERT_EQ( tK.n_rows(), 12u );
        ASSERT_EQ( tK.n_cols(), 12u );

        real tMaxK = 0.0 ;
        for ( uint i = 0; i < 12; ++i )
        {
            for ( uint j = 0; j < 12; ++j )
            {
                tMaxK = std::max( tMaxK, std::abs( tK( i, j ) ) );
            }
        }
        ASSERT_GT( tMaxK, 0.0 );

        // ---- 4. symmetry ----
        // the adjoint-consistency terms make K exactly symmetric
        // ( symmetric weighted interior penalty; ruled intentional
        // 2026-08-28, measured max|K - K^T| = 0.0 on this fixture ).
        // relative tolerance instead of exact zero so a backend swap
        // cannot flake the gate
        real tMaxAsym = 0.0 ;
        for ( uint i = 0; i < 12; ++i )
        {
            for ( uint j = i + 1; j < 12; ++j )
            {
                tMaxAsym = std::max( tMaxAsym,
                        std::abs( tK( i, j ) - tK( j, i ) ) );
            }
        }
        EXPECT_LT( tMaxAsym, 1.0e-12 * tMaxK )
            << "h_ghost must assemble a symmetric operator" ;

        // ---- 2. interface locality + rho-linearity ----
        // at rho = 1e-8 the far columns are ~5 orders below the near ones
        // ( fixture default; measured worst ratio 4.3e-5 )
        real tMaxFar  = 0.0 ;
        real tMinNear = BELFEM_REAL_MAX ;
        for ( uint k = 0; k < 6; ++k )
        {
            tMaxFar  = std::max( tMaxFar,
                    ghost_column_norm( tK, tGhostFarDofs[ k ] ) );
            tMinNear = std::min( tMinNear,
                    ghost_column_norm( tK, tGhostNearDofs[ k ] ) );
        }
        EXPECT_LT( tMaxFar, 1.0e-4 * tMinNear )
            << "far-from-interface dofs must couple only at rho scale" ;

        // 10x the resistivity: far columns are linear in rho_harm and
        // must scale by 10; near columns are k_reg-dominated and must
        // stay put ( h_ghost re-reads element_rho, no re-link needed )
        tStack.set_rho( 1.0e-7 );
        fem::TimestepMatrices tMat10 ;
        tMat10.set_flag( fem::MatrixFlag::K );
        tMat10.initialize( 2 * n );
        fem::maxwell::h_ghost( tCalc, & tMat10 );
        const Matrix< real > & tK10 = tMat10.K();

        for ( uint k = 0; k < 6; ++k )
        {
            const real tFar1  = ghost_column_norm( tK,   tGhostFarDofs[ k ] );
            const real tFar10 = ghost_column_norm( tK10, tGhostFarDofs[ k ] );
            ASSERT_GT( tFar1, 0.0 )
                << "far column " << tGhostFarDofs[ k ] << " is empty" ;
            EXPECT_NEAR( tFar10 / tFar1, 10.0, 0.1 )
                << "far column " << tGhostFarDofs[ k ] ;

            const real tNear1  = ghost_column_norm( tK,   tGhostNearDofs[ k ] );
            const real tNear10 = ghost_column_norm( tK10, tGhostNearDofs[ k ] );
            ASSERT_GT( tNear1, 0.0 )
                << "near column " << tGhostNearDofs[ k ] << " is empty" ;
            EXPECT_NEAR( tNear10 / tNear1, 1.0, 0.01 )
                << "near column " << tGhostNearDofs[ k ] ;
        }

        // ---- 3. annihilation + negative control ----
        // constant tangential field, continuous across the interface:
        // K annihilates it; the master-only half field is not continuous
        // and must NOT be annihilated
        const real tC[ 2 ][ 3 ] = { { 1.0, 0.0, 0.0 },
                                    { 0.0, 1.0, 0.0 } };
        mesh::Element * tPrisms[ 2 ] = { tStack.master(), tStack.slave() };

        for ( uint c = 0; c < 2; ++c )
        {
            Vector< real > tQ( 12, 0.0 );
            for ( uint p = 0; p < 2; ++p )
            {
                for ( uint e = 0; e < 6; ++e )
                {
                    tQ( 6 * p + e ) =
                            ghost_edge_dof( tPrisms[ p ]->edge( e ), tC[ c ] );
                }
            }
            EXPECT_LT( ghost_max_Kq( tK, tQ ), 1.0e-12 * tMaxK )
                << "continuous tangential field not annihilated, c = " << c ;

            // negative control: master half only
            for ( uint j = 6; j < 12; ++j )
            {
                tQ( j ) = 0.0 ;
            }
            EXPECT_GT( ghost_max_Kq( tK, tQ ), 1.0e-3 * tMaxK )
                << "one-sided field must not be annihilated, c = " << c ;
        }
    }
}

//------------------------------------------------------------------------------
// G-operator battery for the thin-shell pair ( QUAD4TS, PENTA6TS )
//------------------------------------------------------------------------------

namespace
{
    // FD constants for the gradient battery ( split tolerances per the
    // TET4/TRI3 round: exact identities on tEpsilon, central differences
    // on their own coarser epsilon )
    const real tDxiG       = 1.0e-4 ;
    const real tEpsilonFDG = 1.0e-8 ;

    /**
     * G-operator battery for one QUAD4TS layer:
     *   1. shape ( ASSERT: a return-mG typo yields the dead 0x0 member )
     *   2. curl tie vs compiled C() at all five eval points
     *   3. eta-derivative of E ( columns 0<->2, 1<->3 share xi, eta = -+1,
     *      E is linear in eta -> exact ), chained with the FIXTURE-side
     *      nabla eta = ( 2/t )( -sin a, cos a ) built from the test's own
     *      angle argument, never from link()'s formula
     *   4. trace == 0 ( the fixture stacks orthogonally for every angle )
     *   5. kernel: q = ( s0 c, s1 c ) -> G*q = 0 ( constant field )
     *   6. negative control: q = ( s0, 0 ) -> |G*q| > tEpsilon
     */
    void
    quad4ts_gradient_battery(
            const real aAngle,
            const int  aFlipRow )
    {
        const real tT = 0.1 ;
        fem::test::TS_TestStack2D tStack( 1, 2.0, tT, aAngle, aFlipRow );

        fem::EF_QUAD4TS tEF ;
        Matrix< real > tXi = quad4ts_eval_points();
        tEF.precompute( tXi );
        tEF.link( tStack.element( 0 ) );

        // 1. shape
        Matrix< real > tG( tEF.G( 0 ) );
        ASSERT_EQ( tG.n_rows(), ( uint ) 4 ) << "G row count" ;
        ASSERT_EQ( tG.n_cols(), ( uint ) 2 ) << "G column count" ;

        // 2. curl tie at every eval point ( G is constant, C recomputes )
        for ( uint g = 0; g < 5; ++g )
        {
            const Matrix< real > & tC = tEF.C( g );
            for ( uint k = 0; k < 2; ++k )
            {
                EXPECT_NEAR( tC( 0, k ), tG( 2, k ) - tG( 1, k ), tEpsilon )
                        << "curl tie, point " << g << " dof " << k ;
            }
        }

        // 3. eta-FD vs E: copy the bottom side first ( E aliases mE )
        const real tNablaEtaFix[ 2 ] =
                { -2.0 / tT * std::sin( aAngle ),
                   2.0 / tT * std::cos( aAngle ) };
        for ( uint p = 0; p < 2; ++p )   // the two xi Gauss stations
        {
            Matrix< real > tBot( tEF.E( p ) );
            const Matrix< real > & tTop = tEF.E( 2 + p );
            for ( uint k = 0; k < 2; ++k )
            {
                for ( uint c = 0; c < 2; ++c )        // field component
                {
                    const real tDEdEta =
                            0.5 * ( tTop( c, k ) - tBot( c, k ) );
                    for ( uint m = 0; m < 2; ++m )    // derivative direction
                    {
                        EXPECT_NEAR( tG( m + 2 * c, k ),
                                tDEdEta * tNablaEtaFix[ m ], tEpsilon )
                                << "G vs dE/deta, station " << p
                                << " dof " << k << " entry ( " << m
                                << ", " << c << " )" ;
                    }
                }
            }
        }

        // 4. trace ( orthogonal stacking )
        for ( uint k = 0; k < 2; ++k )
        {
            EXPECT_NEAR( tG( 0, k ) + tG( 3, k ), 0.0, tEpsilon )
                    << "trace of dof " << k ;
        }

        // 5. kernel + 6. negative control
        real tS[ 12 ];
        tStack.element( 0 )->edge_directions( tS );
        const real tC0 = 0.7 ;
        for ( uint r = 0; r < 4; ++r )
        {
            EXPECT_NEAR( tG( r, 0 ) * tS[ 0 ] * tC0
                       + tG( r, 1 ) * tS[ 1 ] * tC0, 0.0, tEpsilon )
                    << "kernel mode, row " << r ;
        }
        real tMaxNeg = 0.0 ;
        for ( uint r = 0; r < 4; ++r )
        {
            const real tV = std::abs( tG( r, 0 ) * tS[ 0 ] );
            if ( tV > tMaxNeg ) { tMaxNeg = tV ; }
        }
        EXPECT_GT( tMaxNeg, tEpsilon ) << "negative control ( s0, 0 )" ;
    }

//------------------------------------------------------------------------------

    /**
     * G-operator battery for the PENTA6TS reference prism:
     *   1. shape ASSERT ( 9 x 6 )
     *   2. FD vs E on the packed cluster ( col 0 base ( 0.3, 0.25, 0.4 ),
     *      cols 1..6 minus/plus per axis ) through the test-side
     *      NON-transposed Jacobian J( m, d ) = dx_m / dxi_d with columns
     *      ( n0 - n2, n1 - n2, ( t/2 ) n ) — the tau = 0.4 base makes an
     *      F0<->F1 swap visible ( F0 = 0.3, F1 = 0.7 )
     *   3. curl tie vs compiled C() at the base point ( round-off )
     *   4. trace == 0 at the base point
     *   5. kernel: phi = x + 2y + 3z on the bottom corners, local pairs
     *      { 01, 12, 20 }, q_m = s_m c_m copied to the top slots ->
     *      G*q = 0 AND C*q = 0 ( sum c_m = 0 around the closed loop )
     *   6. negative control: q_m = s_m ( all six ) is the Ampere mode ->
     *      max |G*q| > tEpsilon
     */
    void
    penta6ts_gradient_battery( const int aFlipEdge )
    {
        const real tT = 0.1 ;
        fem::test::TS_TestPrism tPrism(
                ElementType::PENTA6TS, tT, aFlipEdge );
        fem::Element * tFem = tPrism.element();
        mesh::Element * tElement = tFem->element();

        // packed FD cluster
        const real tBasePt[ 3 ] = { 0.3, 0.25, 0.4 };
        Matrix< real > tXi( 3, 7 );
        for ( uint i = 0; i < 3; ++i ) { tXi( i, 0 ) = tBasePt[ i ]; }
        uint tCount = 1 ;
        for ( uint p = 0; p < 3; ++p )
        {
            for ( uint s = 0; s < 2; ++s )
            {
                for ( uint i = 0; i < 3; ++i )
                {
                    tXi( i, tCount ) = tBasePt[ i ];
                }
                tXi( p, tCount ) += ( s == 0 ? -tDxiG : tDxiG );
                ++tCount ;
            }
        }

        fem::EF_PENTA6TS tEF ;
        tEF.precompute( tXi );
        tEF.link( tFem );

        // 1. shape
        Matrix< real > tG( tEF.G( 0 ) );
        ASSERT_EQ( tG.n_rows(), ( uint ) 9 ) << "G row count" ;
        ASSERT_EQ( tG.n_cols(), ( uint ) 6 ) << "G column count" ;

        // test-side NON-transposed Jacobian J( m, d ) = dx_m / dxi_d
        Matrix< real > tJ( 3, 3 );
        for ( uint m = 0; m < 3; ++m )
        {
            const real tD0 = ( m == 0 ) ?
                    tElement->node( 0 )->x() - tElement->node( 2 )->x() :
                    ( m == 1 ) ?
                    tElement->node( 0 )->y() - tElement->node( 2 )->y() :
                    tElement->node( 0 )->z() - tElement->node( 2 )->z() ;
            const real tD1 = ( m == 0 ) ?
                    tElement->node( 1 )->x() - tElement->node( 2 )->x() :
                    ( m == 1 ) ?
                    tElement->node( 1 )->y() - tElement->node( 2 )->y() :
                    tElement->node( 1 )->z() - tElement->node( 2 )->z() ;
            tJ( m, 0 ) = tD0 ;
            tJ( m, 1 ) = tD1 ;
        }
        // ( t/2 ) n, n = normalize( col0 x col1 )
        real tN[ 3 ] = {
            tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 2, 0 ) * tJ( 1, 1 ),
            tJ( 2, 0 ) * tJ( 0, 1 ) - tJ( 0, 0 ) * tJ( 2, 1 ),
            tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 1, 0 ) * tJ( 0, 1 ) };
        const real tNn = std::sqrt(
                tN[ 0 ] * tN[ 0 ] + tN[ 1 ] * tN[ 1 ] + tN[ 2 ] * tN[ 2 ] );
        for ( uint m = 0; m < 3; ++m )
        {
            tJ( m, 2 ) = 0.5 * tT * tN[ m ] / tNn ;
        }

        // cofactor inverse ( same formula as test_EdgeFunctions'
        // test_jacobian, replicated here — this TU has no such helper )
        Matrix< real > tInvJ( 3, 3 );
        const real tDet =
              tJ( 0, 0 ) * ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) )
            - tJ( 0, 1 ) * ( tJ( 1, 0 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 0 ) )
            + tJ( 0, 2 ) * ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) );
        tInvJ( 0, 0 ) = ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) ) / tDet ;
        tInvJ( 0, 1 ) = ( tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 ) ) / tDet ;
        tInvJ( 0, 2 ) = ( tJ( 0, 1 ) * tJ( 1, 2 ) - tJ( 0, 2 ) * tJ( 1, 1 ) ) / tDet ;
        tInvJ( 1, 0 ) = ( tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 ) ) / tDet ;
        tInvJ( 1, 1 ) = ( tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 ) ) / tDet ;
        tInvJ( 1, 2 ) = ( tJ( 0, 2 ) * tJ( 1, 0 ) - tJ( 0, 0 ) * tJ( 1, 2 ) ) / tDet ;
        tInvJ( 2, 0 ) = ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) ) / tDet ;
        tInvJ( 2, 1 ) = ( tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 ) ) / tDet ;
        tInvJ( 2, 2 ) = ( tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 ) ) / tDet ;

        // 2. FD vs E ( copy minus side first — E aliases mE )
        real tMaxFD = 0.0 ;
        Cell< Matrix< real > > tDE ;
        tDE.set_size( 3, {} );
        for ( uint p = 0; p < 3; ++p )
        {
            Matrix< real > tMinus( tEF.E( 1 + 2 * p ) );
            const Matrix< real > & tPlus = tEF.E( 2 + 2 * p );
            tDE( p ).set_size( 3, 6 );
            for ( uint c = 0; c < 3; ++c )
            {
                for ( uint k = 0; k < 6; ++k )
                {
                    tDE( p )( c, k ) =
                            ( tPlus( c, k ) - tMinus( c, k ) )
                            / ( 2.0 * tDxiG );
                }
            }
        }
        for ( uint k = 0; k < 6; ++k )
        {
            for ( uint c = 0; c < 3; ++c )        // field component
            {
                for ( uint m = 0; m < 3; ++m )    // derivative direction
                {
                    real tFD = 0.0 ;
                    for ( uint p = 0; p < 3; ++p )
                    {
                        tFD += tDE( p )( c, k ) * tInvJ( p, m );
                    }
                    EXPECT_NEAR( tG( m + 3 * c, k ), tFD, tEpsilonFDG )
                            << "G vs FD, derivative " << m
                            << " of component " << c << ", dof " << k ;
                    const real tDev = std::abs( tG( m + 3 * c, k ) - tFD );
                    if ( tDev > tMaxFD ) { tMaxFD = tDev ; }
                }
            }
        }

        // 3. curl tie at the base point ( tau = 0.4 != 0: an F0<->F1 swap
        //    is visible here )
        const Matrix< real > & tC = tEF.C( 0 );
        real tMaxTie = 0.0 ;
        for ( uint k = 0; k < 6; ++k )
        {
            real tDev ;
            tDev = std::abs( tC( 0, k ) - ( tG( 7, k ) - tG( 5, k ) ) );
            if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
            tDev = std::abs( tC( 1, k ) - ( tG( 2, k ) - tG( 6, k ) ) );
            if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
            tDev = std::abs( tC( 2, k ) - ( tG( 3, k ) - tG( 1, k ) ) );
            if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
        }
        EXPECT_LT( tMaxTie, tEpsilon ) << "curl tie" ;

        // 3b. curl tie at a NON-base cluster point ( column 6 = tau + h ):
        //     an implementation that ignores aIndex and always reads the
        //     base column would pass every base-point check ( audit
        //     hardening, TS C.3 round )
        {
            Matrix< real > tG6( tEF.G( 6 ) );
            const Matrix< real > & tC6 = tEF.C( 6 );
            for ( uint k = 0; k < 6; ++k )
            {
                EXPECT_NEAR( tC6( 0, k ), tG6( 7, k ) - tG6( 5, k ),
                        tEpsilon ) << "curl tie x at tau+h, dof " << k ;
                EXPECT_NEAR( tC6( 1, k ), tG6( 2, k ) - tG6( 6, k ),
                        tEpsilon ) << "curl tie y at tau+h, dof " << k ;
                EXPECT_NEAR( tC6( 2, k ), tG6( 3, k ) - tG6( 1, k ),
                        tEpsilon ) << "curl tie z at tau+h, dof " << k ;
            }
            // and the tau-dependent entries must actually move between the
            // two points ( F0/F1 differ by tDxiG/2 in tau )
            real tMaxMove = 0.0 ;
            for ( uint r = 0; r < 9; ++r )
            {
                for ( uint k = 0; k < 6; ++k )
                {
                    const real tV = std::abs( tG6( r, k ) - tG( r, k ) );
                    if ( tV > tMaxMove ) { tMaxMove = tV ; }
                }
            }
            EXPECT_GT( tMaxMove, tEpsilon )
                    << "G must depend on the integration point" ;
        }

        // 4. trace at the base point
        for ( uint k = 0; k < 6; ++k )
        {
            EXPECT_NEAR( tG( 0, k ) + tG( 4, k ) + tG( 8, k ), 0.0,
                    tEpsilon ) << "trace of dof " << k ;
        }

        // 5. kernel mode: phi = x + 2y + 3z, local bottom pairs
        //    { 01, 12, 20 }, q_m = s_m c_m copied to the top slots
        real tS[ 12 ];
        tFem->edge_directions( tS );
        const uint tPairs[ 3 ][ 2 ] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
        real tQ[ 6 ];
        for ( uint m = 0; m < 3; ++m )
        {
            const mesh::Node * tA = tElement->node( tPairs[ m ][ 0 ] );
            const mesh::Node * tB = tElement->node( tPairs[ m ][ 1 ] );
            const real tCm =
                      ( tB->x() - tA->x() )
                + 2.0 * ( tB->y() - tA->y() )
                + 3.0 * ( tB->z() - tA->z() );
            tQ[ m ]     = tS[ m ]     * tCm ;
            tQ[ m + 3 ] = tS[ m + 3 ] * tCm ;
        }
        for ( uint r = 0; r < 9; ++r )
        {
            real tGq = 0.0 ;
            for ( uint e = 0; e < 6; ++e )
            {
                tGq += tG( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tGq, 0.0, tEpsilon )
                    << "kernel mode, row " << r ;
        }
        for ( uint r = 0; r < 3; ++r )
        {
            real tCq = 0.0 ;
            for ( uint e = 0; e < 6; ++e )
            {
                tCq += tC( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tCq, 0.0, tEpsilon )
                    << "kernel mode is curl-free, row " << r ;
        }

        // 6. negative control: all-s dofs = the Ampere mode
        real tMaxNeg = 0.0 ;
        for ( uint r = 0; r < 9; ++r )
        {
            real tGq = 0.0 ;
            for ( uint e = 0; e < 6; ++e )
            {
                tGq += tG( r, e ) * tS[ e ];
            }
            const real tV = std::abs( tGq );
            if ( tV > tMaxNeg ) { tMaxNeg = tV ; }
        }
        EXPECT_GT( tMaxNeg, tEpsilon ) << "negative control ( all-s )" ;

        std::cout << "    [ G battery ] max |G - FD(E)| = " << tMaxFD
                  << ", max curl-tie dev = " << tMaxTie << std::endl ;
    }
}

TEST( InterfaceOrientation, Quad4TsGradient )
{
    for ( real tAngle : { 0.0, 0.7 } )
    {
        quad4ts_gradient_battery( tAngle, -1 );
    }
}

TEST( InterfaceOrientation, Quad4TsGradientFlipped )
{
    // flip row 1 = the top edge of the single layer
    quad4ts_gradient_battery( 0.0, 1 );
}

TEST( InterfaceOrientation, Penta6TsGradient )
{
    penta6ts_gradient_battery( -1 );
}

TEST( InterfaceOrientation, Penta6TsGradientFlippedBottom )
{
    // flip a BOTTOM edge ( 1 ): bottom columns use mS[ e ], top columns
    // mS[ e + 3 ] — a fill that applies the bottom sign to both faces
    // hides behind a top-edge-only flip
    penta6ts_gradient_battery( 1 );
}

//------------------------------------------------------------------------------
// G-operator battery for the final pair ( HEX8TS, HEX8TB )
//------------------------------------------------------------------------------

namespace
{
    /**
     * G-operator battery for the HEX8TS thin-shell prism.
     *
     * aRectangle = true replaces the default mid-surface with a ROTATED
     * RECTANGLE ( theta = 0.35, sides 1.3 x 0.8 ): the in-plane frame is
     * constant and orthogonal, so term1 G is the exact gradient ( FD
     * gate ) and div e == 0 ( trace gate ). The default fixture quad is
     * NOT a parallelogram — the per-point nablas make term1 incomplete
     * there by the documented class scope — so with aRectangle = false
     * only the convention-exact gates run: shape, curl tie ( exact by
     * construction on every geometry ), G-must-move.
     *
     *   1. shape ASSERT ( 9 x 8 )
     *   2. curl tie vs compiled C() at the base point and at cluster
     *      column 6, plus G-must-move between them
     *   3. [ rectangle ] FD vs E on the packed cluster through the
     *      test-side Jacobian built from the element's own corner nodes
     *   4. [ rectangle ] trace == 0
     *   5. [ rectangle ] kernel: IN-PLANE phi = x + 2y ( the 8-dof shell
     *      space has no vertical edges, so only in-plane constant fields
     *      are representable ), q_k = s_k ( phi( B ) - phi( A ) ) along
     *      the canonical edges -> G*q = 0 and C*q = 0
     *   6. [ rectangle ] negative control: all-s is the double-loop
     *      Ampere mode -> max |G*q| > tEpsilon
     */
    void
    hex8ts_gradient_battery( const int aFlipEdge, const bool aRectangle )
    {
        const real tT = 0.1 ;

        // rotated rectangle in the z = 0 plane, counterclockwise from +z
        Matrix< real > tCorners( 2, 4 );
        {
            const real tTheta = 0.35 ;
            const real tP0[ 2 ] = { 0.2, 0.1 };
            const real tA[ 2 ]  = {  1.3 * std::cos( tTheta ),
                                     1.3 * std::sin( tTheta ) };
            const real tB[ 2 ]  = { -0.8 * std::sin( tTheta ),
                                     0.8 * std::cos( tTheta ) };
            for ( uint i = 0; i < 2; ++i )
            {
                tCorners( i, 0 ) = tP0[ i ];
                tCorners( i, 1 ) = tP0[ i ] + tA[ i ];
                tCorners( i, 2 ) = tP0[ i ] + tA[ i ] + tB[ i ];
                tCorners( i, 3 ) = tP0[ i ] + tB[ i ];
            }
        }

        fem::test::TS_TestPrism tPrism( ElementType::HEX8TS, tT, aFlipEdge,
                aRectangle ? &tCorners : nullptr );
        fem::Element * tFem = tPrism.element();
        mesh::Element * tElement = tFem->element();

        // packed FD cluster ( PENTA6TS pattern )
        const real tBasePt[ 3 ] = { 0.3, 0.25, 0.4 };
        Matrix< real > tXi( 3, 7 );
        for ( uint i = 0; i < 3; ++i ) { tXi( i, 0 ) = tBasePt[ i ]; }
        uint tCount = 1 ;
        for ( uint p = 0; p < 3; ++p )
        {
            for ( uint s = 0; s < 2; ++s )
            {
                for ( uint i = 0; i < 3; ++i )
                {
                    tXi( i, tCount ) = tBasePt[ i ];
                }
                tXi( p, tCount ) += ( s == 0 ? -tDxiG : tDxiG );
                ++tCount ;
            }
        }

        fem::EF_HEX8TS tEF ;
        tEF.precompute( tXi );
        tEF.link( tFem );

        // 1. shape
        Matrix< real > tG( tEF.G( 0 ) );
        ASSERT_EQ( tG.n_rows(), ( uint ) 9 ) << "G row count" ;
        ASSERT_EQ( tG.n_cols(), ( uint ) 8 ) << "G column count" ;

        // 2. curl tie at the base point and at cluster column 6, plus
        //    G-must-move ( columns differ in zeta -> the F factors move )
        {
            Matrix< real > tC0( tEF.C( 0 ) );
            for ( uint k = 0; k < 8; ++k )
            {
                EXPECT_NEAR( tC0( 0, k ), tG( 7, k ) - tG( 5, k ), tEpsilon )
                        << "curl tie x at base, dof " << k ;
                EXPECT_NEAR( tC0( 1, k ), tG( 2, k ) - tG( 6, k ), tEpsilon )
                        << "curl tie y at base, dof " << k ;
                EXPECT_NEAR( tC0( 2, k ), tG( 3, k ) - tG( 1, k ), tEpsilon )
                        << "curl tie z at base, dof " << k ;
            }

            Matrix< real > tG6( tEF.G( 6 ) );
            const Matrix< real > & tC6 = tEF.C( 6 );
            for ( uint k = 0; k < 8; ++k )
            {
                EXPECT_NEAR( tC6( 0, k ), tG6( 7, k ) - tG6( 5, k ), tEpsilon )
                        << "curl tie x at zeta+h, dof " << k ;
                EXPECT_NEAR( tC6( 1, k ), tG6( 2, k ) - tG6( 6, k ), tEpsilon )
                        << "curl tie y at zeta+h, dof " << k ;
                EXPECT_NEAR( tC6( 2, k ), tG6( 3, k ) - tG6( 1, k ), tEpsilon )
                        << "curl tie z at zeta+h, dof " << k ;
            }

            real tMaxMove = 0.0 ;
            for ( uint r = 0; r < 9; ++r )
            {
                for ( uint k = 0; k < 8; ++k )
                {
                    const real tV = std::abs( tG6( r, k ) - tG( r, k ) );
                    if ( tV > tMaxMove ) { tMaxMove = tV ; }
                }
            }
            EXPECT_GT( tMaxMove, tEpsilon )
                    << "G must depend on the integration point" ;
        }

        if ( ! aRectangle )
        {
            return ;
        }

        // test-side Jacobian from the element's own corner nodes:
        // affine parallelogram map, dx/dxi = ( n1 - n0 )/2,
        // dx/deta = ( n3 - n0 )/2, dx/dzeta = ( t/2 ) n
        Matrix< real > tJ( 3, 3 );
        tJ( 0, 0 ) = 0.5 * ( tElement->node( 1 )->x() - tElement->node( 0 )->x() );
        tJ( 1, 0 ) = 0.5 * ( tElement->node( 1 )->y() - tElement->node( 0 )->y() );
        tJ( 2, 0 ) = 0.5 * ( tElement->node( 1 )->z() - tElement->node( 0 )->z() );
        tJ( 0, 1 ) = 0.5 * ( tElement->node( 3 )->x() - tElement->node( 0 )->x() );
        tJ( 1, 1 ) = 0.5 * ( tElement->node( 3 )->y() - tElement->node( 0 )->y() );
        tJ( 2, 1 ) = 0.5 * ( tElement->node( 3 )->z() - tElement->node( 0 )->z() );
        real tN[ 3 ] = {
            tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 2, 0 ) * tJ( 1, 1 ),
            tJ( 2, 0 ) * tJ( 0, 1 ) - tJ( 0, 0 ) * tJ( 2, 1 ),
            tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 1, 0 ) * tJ( 0, 1 ) };
        const real tNn = std::sqrt(
                tN[ 0 ] * tN[ 0 ] + tN[ 1 ] * tN[ 1 ] + tN[ 2 ] * tN[ 2 ] );
        for ( uint m = 0; m < 3; ++m )
        {
            tJ( m, 2 ) = 0.5 * tT * tN[ m ] / tNn ;
        }

        Matrix< real > tInvJ( 3, 3 );
        const real tDet =
              tJ( 0, 0 ) * ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) )
            - tJ( 0, 1 ) * ( tJ( 1, 0 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 0 ) )
            + tJ( 0, 2 ) * ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) );
        tInvJ( 0, 0 ) = ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) ) / tDet ;
        tInvJ( 0, 1 ) = ( tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 ) ) / tDet ;
        tInvJ( 0, 2 ) = ( tJ( 0, 1 ) * tJ( 1, 2 ) - tJ( 0, 2 ) * tJ( 1, 1 ) ) / tDet ;
        tInvJ( 1, 0 ) = ( tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 ) ) / tDet ;
        tInvJ( 1, 1 ) = ( tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 ) ) / tDet ;
        tInvJ( 1, 2 ) = ( tJ( 0, 2 ) * tJ( 1, 0 ) - tJ( 0, 0 ) * tJ( 1, 2 ) ) / tDet ;
        tInvJ( 2, 0 ) = ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) ) / tDet ;
        tInvJ( 2, 1 ) = ( tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 ) ) / tDet ;
        tInvJ( 2, 2 ) = ( tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 ) ) / tDet ;

        // 3. FD vs E ( copy minus side first — E aliases mE )
        real tMaxFD = 0.0 ;
        Cell< Matrix< real > > tDE ;
        tDE.set_size( 3, {} );
        for ( uint p = 0; p < 3; ++p )
        {
            Matrix< real > tMinus( tEF.E( 1 + 2 * p ) );
            const Matrix< real > & tPlus = tEF.E( 2 + 2 * p );
            tDE( p ).set_size( 3, 8 );
            for ( uint c = 0; c < 3; ++c )
            {
                for ( uint k = 0; k < 8; ++k )
                {
                    tDE( p )( c, k ) =
                            ( tPlus( c, k ) - tMinus( c, k ) )
                            / ( 2.0 * tDxiG );
                }
            }
        }
        for ( uint k = 0; k < 8; ++k )
        {
            for ( uint c = 0; c < 3; ++c )        // field component
            {
                for ( uint m = 0; m < 3; ++m )    // derivative direction
                {
                    real tFD = 0.0 ;
                    for ( uint p = 0; p < 3; ++p )
                    {
                        tFD += tDE( p )( c, k ) * tInvJ( p, m );
                    }
                    EXPECT_NEAR( tG( m + 3 * c, k ), tFD, tEpsilonFDG )
                            << "G vs FD, derivative " << m
                            << " of component " << c << ", dof " << k ;
                    const real tDev = std::abs( tG( m + 3 * c, k ) - tFD );
                    if ( tDev > tMaxFD ) { tMaxFD = tDev ; }
                }
            }
        }

        // 4. trace ( orthogonal constant frame )
        for ( uint k = 0; k < 8; ++k )
        {
            EXPECT_NEAR( tG( 0, k ) + tG( 4, k ) + tG( 8, k ), 0.0,
                    tEpsilon ) << "trace of dof " << k ;
        }

        // 5. kernel: in-plane phi = x + 2y, canonical bottom pairs
        //    { 01, 12, 23, 30 }, top slots copy the bottom values ( phi is
        //    z-free and the prism is extruded straight up )
        real tS[ 12 ];
        tFem->edge_directions( tS );
        const uint tPairs[ 4 ][ 2 ] =
                { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
        real tQ[ 8 ];
        for ( uint m = 0; m < 4; ++m )
        {
            const mesh::Node * tNa = tElement->node( tPairs[ m ][ 0 ] );
            const mesh::Node * tNb = tElement->node( tPairs[ m ][ 1 ] );
            const real tCm =
                      ( tNb->x() - tNa->x() )
                + 2.0 * ( tNb->y() - tNa->y() );
            tQ[ m ]     = tS[ m ]     * tCm ;
            tQ[ m + 4 ] = tS[ m + 4 ] * tCm ;
        }
        const Matrix< real > & tC = tEF.C( 0 );
        for ( uint r = 0; r < 9; ++r )
        {
            real tGq = 0.0 ;
            for ( uint e = 0; e < 8; ++e )
            {
                tGq += tG( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tGq, 0.0, tEpsilon ) << "kernel mode, row " << r ;
        }
        for ( uint r = 0; r < 3; ++r )
        {
            real tCq = 0.0 ;
            for ( uint e = 0; e < 8; ++e )
            {
                tCq += tC( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tCq, 0.0, tEpsilon )
                    << "kernel mode is curl-free, row " << r ;
        }

        // 6. negative control: all-s = the double-loop Ampere mode
        real tMaxNeg = 0.0 ;
        for ( uint r = 0; r < 9; ++r )
        {
            real tGq = 0.0 ;
            for ( uint e = 0; e < 8; ++e )
            {
                tGq += tG( r, e ) * tS[ e ];
            }
            const real tV = std::abs( tGq );
            if ( tV > tMaxNeg ) { tMaxNeg = tV ; }
        }
        EXPECT_GT( tMaxNeg, tEpsilon ) << "negative control ( all-s )" ;

        std::cout << "    [ G battery ] max |G - FD(E)| = " << tMaxFD
                  << std::endl ;
    }

//------------------------------------------------------------------------------

    /**
     * G-operator battery for the HEX8TB side-connector wall
     * ( TS_TestWall, imposed exact-cuboid frame — affine, so every gate
     * is exact ):
     *   1. shape ASSERT ( 9 x 4 )
     *   2. FD vs E on the packed cluster through the test-side Jacobian
     *      from the wall's own corner nodes
     *   3. curl tie vs compiled C() at the base point and at cluster
     *      column 6 ( zeta + h ), plus G-must-move between them —
     *      xi-shifted points would NOT move G ( F has no xi dependence )
     *   4. trace == 0 at the base point ( orthonormal frame )
     *   5. STOKES on the bottom face ( zeta = -1, normal +n ): per-dof
     *      E-circulation around the loop ( edge 0 forward, edge 1
     *      backward ) vs the 2x2 Gauss face integral of curl . n, taken
     *      from BOTH the compiled C and antisym( G ) — per-dof, because
     *      the coefficient vector ( 1, 1, 0, 0 ) telescopes to the
     *      vacuous 0 = 0. This is the element's first compiled curl-sign
     *      gate.
     *   6. kernel: phi = x + 2y + 3z, q_k = s_k ( phi( B ) - phi( A ) )
     *      along the canonical edges — all four longitudinal edges have
     *      the same edge vector, so E*q is the constant tangential
     *      projection of grad phi -> G*q = 0 and C*q = 0
     *   7. negative control: q = ( s0, 0, 0, 0 ) -> max |G*q| > tEpsilon
     */
    void
    hex8tb_gradient_battery( const int aFlipEdge )
    {
        fem::test::TS_TestWall tWall( 2.0, 0.05, 0.1, aFlipEdge );
        fem::Element * tFem = tWall.element();
        mesh::Element * tElement = tFem->element();

        // packed FD cluster
        const real tBasePt[ 3 ] = { 0.3, 0.25, 0.4 };
        Matrix< real > tXi( 3, 7 );
        for ( uint i = 0; i < 3; ++i ) { tXi( i, 0 ) = tBasePt[ i ]; }
        uint tCount = 1 ;
        for ( uint p = 0; p < 3; ++p )
        {
            for ( uint s = 0; s < 2; ++s )
            {
                for ( uint i = 0; i < 3; ++i )
                {
                    tXi( i, tCount ) = tBasePt[ i ];
                }
                tXi( p, tCount ) += ( s == 0 ? -tDxiG : tDxiG );
                ++tCount ;
            }
        }

        fem::EF_HEX8TB tEF ;
        tEF.precompute( tXi );
        tEF.link( tFem );

        // 1. shape
        Matrix< real > tG( tEF.G( 0 ) );
        ASSERT_EQ( tG.n_rows(), ( uint ) 9 ) << "G row count" ;
        ASSERT_EQ( tG.n_cols(), ( uint ) 4 ) << "G column count" ;

        // test-side Jacobian from the wall's own corner nodes ( cuboid ):
        // dx/dxi = ( n1 - n0 )/2, dx/deta = ( n3 - n0 )/2,
        // dx/dzeta = ( n4 - n0 )/2
        Matrix< real > tJ( 3, 3 );
        tJ( 0, 0 ) = 0.5 * ( tElement->node( 1 )->x() - tElement->node( 0 )->x() );
        tJ( 1, 0 ) = 0.5 * ( tElement->node( 1 )->y() - tElement->node( 0 )->y() );
        tJ( 2, 0 ) = 0.5 * ( tElement->node( 1 )->z() - tElement->node( 0 )->z() );
        tJ( 0, 1 ) = 0.5 * ( tElement->node( 3 )->x() - tElement->node( 0 )->x() );
        tJ( 1, 1 ) = 0.5 * ( tElement->node( 3 )->y() - tElement->node( 0 )->y() );
        tJ( 2, 1 ) = 0.5 * ( tElement->node( 3 )->z() - tElement->node( 0 )->z() );
        tJ( 0, 2 ) = 0.5 * ( tElement->node( 4 )->x() - tElement->node( 0 )->x() );
        tJ( 1, 2 ) = 0.5 * ( tElement->node( 4 )->y() - tElement->node( 0 )->y() );
        tJ( 2, 2 ) = 0.5 * ( tElement->node( 4 )->z() - tElement->node( 0 )->z() );

        Matrix< real > tInvJ( 3, 3 );
        const real tDet =
              tJ( 0, 0 ) * ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) )
            - tJ( 0, 1 ) * ( tJ( 1, 0 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 0 ) )
            + tJ( 0, 2 ) * ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) );
        tInvJ( 0, 0 ) = ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) ) / tDet ;
        tInvJ( 0, 1 ) = ( tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 ) ) / tDet ;
        tInvJ( 0, 2 ) = ( tJ( 0, 1 ) * tJ( 1, 2 ) - tJ( 0, 2 ) * tJ( 1, 1 ) ) / tDet ;
        tInvJ( 1, 0 ) = ( tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 ) ) / tDet ;
        tInvJ( 1, 1 ) = ( tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 ) ) / tDet ;
        tInvJ( 1, 2 ) = ( tJ( 0, 2 ) * tJ( 1, 0 ) - tJ( 0, 0 ) * tJ( 1, 2 ) ) / tDet ;
        tInvJ( 2, 0 ) = ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) ) / tDet ;
        tInvJ( 2, 1 ) = ( tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 ) ) / tDet ;
        tInvJ( 2, 2 ) = ( tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 ) ) / tDet ;

        // 2. FD vs E
        real tMaxFD = 0.0 ;
        Cell< Matrix< real > > tDE ;
        tDE.set_size( 3, {} );
        for ( uint p = 0; p < 3; ++p )
        {
            Matrix< real > tMinus( tEF.E( 1 + 2 * p ) );
            const Matrix< real > & tPlus = tEF.E( 2 + 2 * p );
            tDE( p ).set_size( 3, 4 );
            for ( uint c = 0; c < 3; ++c )
            {
                for ( uint k = 0; k < 4; ++k )
                {
                    tDE( p )( c, k ) =
                            ( tPlus( c, k ) - tMinus( c, k ) )
                            / ( 2.0 * tDxiG );
                }
            }
        }
        for ( uint k = 0; k < 4; ++k )
        {
            for ( uint c = 0; c < 3; ++c )        // field component
            {
                for ( uint m = 0; m < 3; ++m )    // derivative direction
                {
                    real tFD = 0.0 ;
                    for ( uint p = 0; p < 3; ++p )
                    {
                        tFD += tDE( p )( c, k ) * tInvJ( p, m );
                    }
                    EXPECT_NEAR( tG( m + 3 * c, k ), tFD, tEpsilonFDG )
                            << "G vs FD, derivative " << m
                            << " of component " << c << ", dof " << k ;
                    const real tDev = std::abs( tG( m + 3 * c, k ) - tFD );
                    if ( tDev > tMaxFD ) { tMaxFD = tDev ; }
                }
            }
        }

        // 3. curl tie at the base point and at zeta + h, plus G-must-move
        {
            Matrix< real > tC0( tEF.C( 0 ) );
            for ( uint k = 0; k < 4; ++k )
            {
                EXPECT_NEAR( tC0( 0, k ), tG( 7, k ) - tG( 5, k ), tEpsilon )
                        << "curl tie x at base, dof " << k ;
                EXPECT_NEAR( tC0( 1, k ), tG( 2, k ) - tG( 6, k ), tEpsilon )
                        << "curl tie y at base, dof " << k ;
                EXPECT_NEAR( tC0( 2, k ), tG( 3, k ) - tG( 1, k ), tEpsilon )
                        << "curl tie z at base, dof " << k ;
            }

            Matrix< real > tG6( tEF.G( 6 ) );
            const Matrix< real > & tC6 = tEF.C( 6 );
            for ( uint k = 0; k < 4; ++k )
            {
                EXPECT_NEAR( tC6( 0, k ), tG6( 7, k ) - tG6( 5, k ), tEpsilon )
                        << "curl tie x at zeta+h, dof " << k ;
                EXPECT_NEAR( tC6( 1, k ), tG6( 2, k ) - tG6( 6, k ), tEpsilon )
                        << "curl tie y at zeta+h, dof " << k ;
                EXPECT_NEAR( tC6( 2, k ), tG6( 3, k ) - tG6( 1, k ), tEpsilon )
                        << "curl tie z at zeta+h, dof " << k ;
            }

            real tMaxMove = 0.0 ;
            for ( uint r = 0; r < 9; ++r )
            {
                for ( uint k = 0; k < 4; ++k )
                {
                    const real tV = std::abs( tG6( r, k ) - tG( r, k ) );
                    if ( tV > tMaxMove ) { tMaxMove = tV ; }
                }
            }
            EXPECT_GT( tMaxMove, tEpsilon )
                    << "G must depend on the integration point" ;
        }

        // 4. trace ( orthonormal frame )
        for ( uint k = 0; k < 4; ++k )
        {
            EXPECT_NEAR( tG( 0, k ) + tG( 4, k ) + tG( 8, k ), 0.0,
                    tEpsilon ) << "trace of dof " << k ;
        }

        // 5. Stokes on the bottom face, per dof: LHS from the compiled
        //    E-circulations ( canonical edges 0 and 1 both run +xi, so
        //    the CCW loop wrt the face normal is edge 0 minus edge 1 );
        //    RHS from 2x2 Gauss on zeta = -1, both C . n and
        //    antisym( G ) . n
        {
            Matrix< real > tCirc ;
            prism_circulation_matrix< fem::EF_HEX8TB >(
                    ElementType::HEX8TS, tFem, tCirc );

            // face normal and area from the corner nodes ( rectangle )
            real tU[ 3 ] = {
                tElement->node( 1 )->x() - tElement->node( 0 )->x(),
                tElement->node( 1 )->y() - tElement->node( 0 )->y(),
                tElement->node( 1 )->z() - tElement->node( 0 )->z() };
            real tV[ 3 ] = {
                tElement->node( 3 )->x() - tElement->node( 0 )->x(),
                tElement->node( 3 )->y() - tElement->node( 0 )->y(),
                tElement->node( 3 )->z() - tElement->node( 0 )->z() };
            real tNf[ 3 ] = {
                tU[ 1 ] * tV[ 2 ] - tU[ 2 ] * tV[ 1 ],
                tU[ 2 ] * tV[ 0 ] - tU[ 0 ] * tV[ 2 ],
                tU[ 0 ] * tV[ 1 ] - tU[ 1 ] * tV[ 0 ] };
            const real tArea = std::sqrt(
                    tNf[ 0 ] * tNf[ 0 ] + tNf[ 1 ] * tNf[ 1 ]
                  + tNf[ 2 ] * tNf[ 2 ] );
            for ( uint i = 0; i < 3; ++i ) { tNf[ i ] /= tArea ; }

            // 2x2 Gauss on the face
            const real tGp = 1.0 / std::sqrt( 3.0 );
            Matrix< real > tXiF( 3, 4 );
            const real tSgn[ 4 ][ 2 ] =
                    { { -1., -1. }, { 1., -1. }, { 1., 1. }, { -1., 1. } };
            for ( uint g = 0; g < 4; ++g )
            {
                tXiF( 0, g ) = tSgn[ g ][ 0 ] * tGp ;
                tXiF( 1, g ) = tSgn[ g ][ 1 ] * tGp ;
                tXiF( 2, g ) = -1.0 ;
            }

            fem::EF_HEX8TB tEFF ;
            tEFF.precompute( tXiF );
            tEFF.link( tFem );

            real tIntC[ 4 ] = { 0., 0., 0., 0. };
            real tIntG[ 4 ] = { 0., 0., 0., 0. };
            for ( uint g = 0; g < 4; ++g )
            {
                Matrix< real > tCf( tEFF.C( g ) );
                const Matrix< real > & tGf = tEFF.G( g );
                for ( uint k = 0; k < 4; ++k )
                {
                    tIntC[ k ] +=
                            ( tCf( 0, k ) * tNf[ 0 ]
                            + tCf( 1, k ) * tNf[ 1 ]
                            + tCf( 2, k ) * tNf[ 2 ] ) * 0.25 * tArea ;
                    const real tCurlG[ 3 ] = {
                            tGf( 7, k ) - tGf( 5, k ),
                            tGf( 2, k ) - tGf( 6, k ),
                            tGf( 3, k ) - tGf( 1, k ) };
                    tIntG[ k ] +=
                            ( tCurlG[ 0 ] * tNf[ 0 ]
                            + tCurlG[ 1 ] * tNf[ 1 ]
                            + tCurlG[ 2 ] * tNf[ 2 ] ) * 0.25 * tArea ;
                }
            }
            // the sign content sits in dofs 0 and 1 ( loop = +-1 ); the
            // top dofs 2, 3 are 0 = 0 on this face ( F and dF/deta both
            // vanish at zeta = -1 ) — their columns are gated by FD and
            // the curl tie instead
            for ( uint k = 0; k < 4; ++k )
            {
                const real tLoop = tCirc( 0, k ) - tCirc( 1, k );
                EXPECT_NEAR( tIntC[ k ], tLoop, tEpsilon )
                        << "Stokes vs C, dof " << k ;
                EXPECT_NEAR( tIntG[ k ], tLoop, tEpsilon )
                        << "Stokes vs antisym( G ), dof " << k ;
            }
        }

        // 6. kernel: phi = x + 2y + 3z along the canonical edges
        real tS[ 12 ];
        tFem->edge_directions( tS );
        real tQ[ 4 ];
        Cell< mesh::Node * > tEdgeNodes ;
        for ( uint m = 0; m < 4; ++m )
        {
            tElement->get_nodes_of_edge( m, tEdgeNodes );
            const mesh::Node * tNa = tEdgeNodes( 0 );
            const mesh::Node * tNb = tEdgeNodes( 1 );
            const real tCm =
                      ( tNb->x() - tNa->x() )
                + 2.0 * ( tNb->y() - tNa->y() )
                + 3.0 * ( tNb->z() - tNa->z() );
            tQ[ m ] = tS[ m ] * tCm ;
        }
        const Matrix< real > & tC = tEF.C( 0 );
        for ( uint r = 0; r < 9; ++r )
        {
            real tGq = 0.0 ;
            for ( uint e = 0; e < 4; ++e )
            {
                tGq += tG( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tGq, 0.0, tEpsilon ) << "kernel mode, row " << r ;
        }
        for ( uint r = 0; r < 3; ++r )
        {
            real tCq = 0.0 ;
            for ( uint e = 0; e < 4; ++e )
            {
                tCq += tC( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tCq, 0.0, tEpsilon )
                    << "kernel mode is curl-free, row " << r ;
        }

        // 7. negative control
        real tMaxNeg = 0.0 ;
        for ( uint r = 0; r < 9; ++r )
        {
            const real tV = std::abs( tG( r, 0 ) * tS[ 0 ] );
            if ( tV > tMaxNeg ) { tMaxNeg = tV ; }
        }
        EXPECT_GT( tMaxNeg, tEpsilon ) << "negative control ( s0, 0, 0, 0 )" ;

        std::cout << "    [ G battery ] max |G - FD(E)| = " << tMaxFD
                  << std::endl ;
    }
}

TEST( InterfaceOrientation, Hex8TsGradient )
{
    // rotated rectangle: full battery ( FD, trace, kernel, controls )
    hex8ts_gradient_battery( -1, true );
}

TEST( InterfaceOrientation, Hex8TsGradientGeneralQuad )
{
    // default general-quad fixture: convention-exact gates only ( the
    // per-point nablas make term1 G incomplete here by documented scope )
    hex8ts_gradient_battery( -1, false );
}

TEST( InterfaceOrientation, Hex8TsGradientFlippedBottom )
{
    // flip a BOTTOM edge ( 1 ): bottom columns use mS[ e ], top columns
    // mS[ e + 4 ]
    hex8ts_gradient_battery( 1, true );
}

TEST( InterfaceOrientation, Hex8TsGradientFlippedTop )
{
    // flip a TOP edge ( 5 ): with only the bottom flip, a fill that
    // drops mS from the top columns alone would survive the battery
    // ( C.3 audit round, both auditors )
    hex8ts_gradient_battery( 5, true );
}

TEST( InterfaceOrientation, Hex8TbGradient )
{
    hex8tb_gradient_battery( -1 );
}

TEST( InterfaceOrientation, Hex8TbGradientFlipped )
{
    // flip edge 1 ( the second bottom longitudinal edge, canonically
    // node 3 -> node 2 ): its column must carry the sign through E, C,
    // G and both sides of the Stokes gate
    hex8tb_gradient_battery( 1 );
}

//------------------------------------------------------------------------------
// Thin-shell normal-field recovery next to an h-conductor
//------------------------------------------------------------------------------

// compute_hn recovers the shell's normal field as the average of the two
// volume traces, projected onto the facet normal ( cl_FEM_Calculator.hpp,
// compute_h_trace ). A volume side that is a phi-region contributes
// -grad phi; a DomainType::Conductor side has no potential — its nodal
// phi is bookkeeping — and contributes its own Nedelec trace E * q as the
// integration-weighted mean over the facet rule. The fixtures set the
// facet's master/slave explicitly, so all three kind combinations run:
// conductor/air ( the corc_solder case ), conductor/conductor ( the only
// production route to a conductor slave ), air/air ( the original path ).
// Fields are seeded directly into the edge_h and phi mesh fields, which
// is the live storage the solve reads.
namespace
{
    real
    dot3( const Vector< real > & aA, const real aB[ 3 ] )
    {
        return aA( 0 ) * aB[ 0 ] + aA( 1 ) * aB[ 1 ] + aA( 2 ) * aB[ 2 ];
    }

    /**
     * conductor master / air slave ( corc_solder ): uniform h = a on both
     * sides, a deliberately nonphysical phi on the conductor-only node
     */
    void
    conductor_air_battery( const int aFlipEdge )
    {
        fem::test::TS_TestConductorShell tStack(
                DomainType::Conductor, DomainType::Air, aFlipEdge );

        const real a[ 3 ] = { 0.3, -0.2, 0.7 };
        const real b[ 3 ] = { 0.0,  0.0, 0.0 };

        tStack.seed_edge_field( tStack.tet_a(), a, b );
        tStack.seed_phi_uniform( tStack.tet_b(), a );

        // the solder-only node of corc_solder: phi is never solved there
        tStack.set_phi( tStack.apex_a(), 1.0e3 );

        fem::Calculator * tCalc = tStack.layer_calculator();
        const Vector< real > & hn = fem::compute_hn( tCalc, 0 );

        const real n[ 3 ] = { 0.0, 0.0, 1.0 };

        // hn = ( a . n ) n, independent of the sign of n
        EXPECT_NEAR( hn( 0 ), 0.0, 1.0e-10 );
        EXPECT_NEAR( hn( 1 ), 0.0, 1.0e-10 );
        EXPECT_NEAR( hn( 2 ), a[ 2 ], 1.0e-10 );

        // O1 agreement gate: the two full traces agree ( the h-side gives
        // a exactly for a field in the Nedelec space, the phi-side gives
        // -grad phi = a ), so a broken h/phi coupling would show here
        const Vector< real > & hm = tCalc->vector( "hm" );
        const Vector< real > & hs = tCalc->vector( "hs" );
        for ( uint i = 0; i < 3; ++i )
        {
            EXPECT_NEAR( hm( i ), a[ i ], 1.0e-10 );
            EXPECT_NEAR( hs( i ), a[ i ], 1.0e-10 );
        }

        // negative control: the old route, -grad phi of the conductor, is
        // dominated by the nonphysical node ( 1e3 over the tet height )
        Vector< real > tPhiA( 4 );
        const Vector< real > & tPhi = tStack.mesh()->field_data( "phi" );
        for ( uint k = 0; k < 4; ++k )
        {
            tPhiA( k ) = tPhi( tStack.tet_a()->node( k )->index() );
        }
        bool tMasterIsConductor = false ;
        bool tSlaveIsConductor  = false ;
        Vector< real > tPhiM ;
        Vector< real > tPhiS ;
        fem::Calculator * tNormalCalc = tCalc->get_normal_calculator(
                tPhiM, tPhiS, tMasterIsConductor, tSlaveIsConductor );
        EXPECT_TRUE( tMasterIsConductor );
        EXPECT_FALSE( tSlaveIsConductor );

        Vector< real > tOld ;
        tOld = -1.0 * tNormalCalc->Bm( 0 ) * tPhiA ;
        EXPECT_GT( std::abs( dot3( tOld, n ) - a[ 2 ] ), 1.0e2 )
                << "the phi-gradient route must fail on the conductor side" ;
    }

    /**
     * conductor / conductor: h = a + b x r on both sides ( in the Nedelec
     * space, normal trace linear over the facet ), nonphysical phi on
     * every node. Checks the weighted mean against the centroid value,
     * the O4 gate ( in 3-D point 0 of the default 7-point rule is the
     * centroid ), and that the trace actually varies over the facet
     */
    void
    conductor_conductor_battery( const int aFlipEdge )
    {
        fem::test::TS_TestConductorShell tStack(
                DomainType::Conductor, DomainType::Conductor, aFlipEdge );

        const real a[ 3 ] = { 0.3, -0.2, 0.7 };
        const real b[ 3 ] = { 0.5, -0.4, 0.2 };

        tStack.seed_edge_field( tStack.tet_a(), a, b );
        tStack.seed_edge_field( tStack.tet_b(), a, b );

        // no phi anywhere on a conductor/conductor shell
        for ( mesh::Node * tNode : tStack.mesh()->nodes() )
        {
            tStack.set_phi( tNode, 1.0e3 );
        }

        // n . h at the facet centroid, n = z
        real c[ 3 ];
        tStack.facet_centroid( c );
        const real tExpected = a[ 2 ] + ( b[ 0 ] * c[ 1 ] - b[ 1 ] * c[ 0 ] );

        fem::Calculator * tCalc = tStack.layer_calculator();
        const Vector< real > & hn = fem::compute_hn( tCalc, 0 );

        EXPECT_NEAR( hn( 0 ), 0.0, 1.0e-10 );
        EXPECT_NEAR( hn( 1 ), 0.0, 1.0e-10 );
        EXPECT_NEAR( hn( 2 ), tExpected, 1.0e-10 );

        // both traces come from the edge field now
        bool tMasterIsConductor = false ;
        bool tSlaveIsConductor  = false ;
        Vector< real > tPhiM ;
        Vector< real > tPhiS ;
        fem::Calculator * tNormalCalc = tCalc->get_normal_calculator(
                tPhiM, tPhiS, tMasterIsConductor, tSlaveIsConductor );
        EXPECT_TRUE( tMasterIsConductor );
        EXPECT_TRUE( tSlaveIsConductor );

        // per-point normal traces of the master side
        const Vector< real > & q = tNormalCalc->nedelec_data_master_h();
        const Vector< real > & w = tNormalCalc->integration()->weights();
        const uint tNumPoints = tNormalCalc->num_intpoints();
        const Vector< real > & tNormal = tNormalCalc->normal( 0 );
        const real tSign = tNormal( 2 ) > 0.0 ? 1.0 : -1.0 ;

        real tMean = 0.0 ;
        real tSum  = 0.0 ;
        real tMaxDev = 0.0 ;
        Vector< real > tH( 3 );
        Vector< real > tV( tNumPoints );
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tH = tNormalCalc->Em( k ) * q ;
            tV( k ) = tSign * ( tNormal( 0 ) * tH( 0 ) + tNormal( 1 ) * tH( 1 ) + tNormal( 2 ) * tH( 2 ) );
            tMean += w( k ) * tV( k );
            tSum  += w( k );
        }
        tMean /= tSum ;
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tMaxDev = std::max( tMaxDev, std::abs( tV( k ) - tMean ) );
        }

        // the weighted mean is the centroid value
        EXPECT_NEAR( tMean, tExpected, 1.0e-10 );

        // O4: in 3-D, k = 0 is the centroid and equals the mean
        EXPECT_NEAR( tV( 0 ), tMean, 1.0e-10 );

        // and the trace is not constant over the facet ( b != 0 )
        EXPECT_GT( tMaxDev, 1.0e-3 );
    }
}

TEST( InterfaceOrientation, ThinShellNormalRecoveryConductorAir )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }
    conductor_air_battery( -1 );
}

TEST( InterfaceOrientation, ThinShellNormalRecoveryConductorAirFlipped )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }
    // reverse local edge 2 of the tets: the dof sign must travel through
    // E and the seeding alike
    conductor_air_battery( 2 );
}

TEST( InterfaceOrientation, ThinShellNormalRecoveryConductorConductor )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }
    conductor_conductor_battery( -1 );
    conductor_conductor_battery( 4 );
}

TEST( InterfaceOrientation, ThinShellNormalRecoveryAirAir )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // the original phi/phi construction, unchanged
    fem::test::TS_TestConductorShell tStack( DomainType::Air, DomainType::Air );

    const real a[ 3 ] = { 0.3, -0.2, 0.7 };
    tStack.seed_phi_uniform( tStack.tet_a(), a );
    tStack.seed_phi_uniform( tStack.tet_b(), a );

    fem::Calculator * tCalc = tStack.layer_calculator();
    const Vector< real > & hn = fem::compute_hn( tCalc, 0 );

    EXPECT_NEAR( hn( 0 ), 0.0, 1.0e-10 );
    EXPECT_NEAR( hn( 1 ), 0.0, 1.0e-10 );
    EXPECT_NEAR( hn( 2 ), a[ 2 ], 1.0e-10 );

    bool tMasterIsConductor = true ;
    bool tSlaveIsConductor  = true ;
    Vector< real > tPhiM ;
    Vector< real > tPhiS ;
    tCalc->get_normal_calculator( tPhiM, tPhiS, tMasterIsConductor, tSlaveIsConductor );
    EXPECT_FALSE( tMasterIsConductor );
    EXPECT_FALSE( tSlaveIsConductor );
}

namespace
{
    /**
     * 2-D conductor / conductor with h = a + b ( -y, x ): the normal trace
     * a_y + b x is linear along the facet and the default line rule has
     * no midpoint, so a single-point sample would be off — the weighted
     * mean is the midpoint value
     */
    void
    conductor_conductor_battery_2d( const int aFlipEdge )
    {
        fem::test::TS_TestConductorShell2D tStack(
                DomainType::Conductor, DomainType::Conductor, aFlipEdge );

        const real a[ 2 ] = { 0.3, 0.7 };
        const real b      = 0.4 ;

        tStack.seed_edge_field( tStack.tri_a(), a, b );
        tStack.seed_edge_field( tStack.tri_b(), a, b );

        for ( mesh::Node * tNode : tStack.mesh()->nodes() )
        {
            tStack.set_phi( tNode, 1.0e3 );
        }

        // n . h at the facet midpoint ( x = L / 2, y = 0 ), n = y
        const real tExpected = a[ 1 ] + b * 0.5 * tStack.length();

        fem::Calculator * tCalc = tStack.layer_calculator();
        const Vector< real > & hn = fem::compute_hn( tCalc, 0 );

        EXPECT_NEAR( hn( 0 ), 0.0, 1.0e-10 );
        EXPECT_NEAR( hn( 1 ), tExpected, 1.0e-10 );

        bool tMasterIsConductor = false ;
        bool tSlaveIsConductor  = false ;
        Vector< real > tPhiM ;
        Vector< real > tPhiS ;
        fem::Calculator * tNormalCalc = tCalc->get_normal_calculator(
                tPhiM, tPhiS, tMasterIsConductor, tSlaveIsConductor );
        EXPECT_TRUE( tMasterIsConductor );
        EXPECT_TRUE( tSlaveIsConductor );

        const Vector< real > & q = tNormalCalc->nedelec_data_slave_h();
        const Vector< real > & w = tNormalCalc->integration()->weights();
        const uint tNumPoints = tNormalCalc->num_intpoints();
        const Vector< real > & tNormal = tNormalCalc->normal( 0 );
        const real tSign = tNormal( 1 ) > 0.0 ? 1.0 : -1.0 ;

        real tMean = 0.0 ;
        real tSum  = 0.0 ;
        real tMaxDev = 0.0 ;
        Vector< real > tH( 2 );
        Vector< real > tV( tNumPoints );
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tH = tNormalCalc->Es( k ) * q ;
            tV( k ) = tSign * ( tNormal( 0 ) * tH( 0 ) + tNormal( 1 ) * tH( 1 ) );
            tMean += w( k ) * tV( k );
            tSum  += w( k );
        }
        tMean /= tSum ;
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tMaxDev = std::max( tMaxDev, std::abs( tV( k ) - tMean ) );
        }

        EXPECT_NEAR( tMean, tExpected, 1.0e-10 );

        // some point of the rule is off the midpoint value: a single-point
        // sample would not do in 2-D
        EXPECT_GT( tMaxDev, 1.0e-3 );
    }
}

TEST( InterfaceOrientation, ThinShellNormalRecoveryConductorConductor2D )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }
    conductor_conductor_battery_2d( -1 );
    conductor_conductor_battery_2d( 1 );
}

TEST( InterfaceOrientation, ThinShellNormalRecoveryConductorAir2D )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::test::TS_TestConductorShell2D tStack( DomainType::Conductor, DomainType::Air );

    const real a[ 2 ] = { 0.3, 0.7 };

    tStack.seed_edge_field( tStack.tri_a(), a, 0.0 );
    tStack.seed_phi_uniform( tStack.tri_b(), a );
    tStack.set_phi( tStack.apex_a(), 1.0e3 );

    fem::Calculator * tCalc = tStack.layer_calculator();
    const Vector< real > & hn = fem::compute_hn( tCalc, 0 );

    EXPECT_NEAR( hn( 0 ), 0.0, 1.0e-10 );
    EXPECT_NEAR( hn( 1 ), a[ 1 ], 1.0e-10 );

    const Vector< real > & hm = tCalc->vector( "hm" );
    const Vector< real > & hs = tCalc->vector( "hs" );
    for ( uint i = 0; i < 2; ++i )
    {
        EXPECT_NEAR( hm( i ), a[ i ], 1.0e-10 );
        EXPECT_NEAR( hs( i ), a[ i ], 1.0e-10 );
    }
}
