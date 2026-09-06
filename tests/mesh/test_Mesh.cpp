/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Mesh class: tensor mesh construction, accessors, map
 * consistency, connectivity on structured grids, scale_mesh, block basics.
 * See: tests_10_mesh.md §4–§9
 *
 * BUG-ME1: Mesh::node(i,j,k) asserts dimension==2 instead of ==3.
 *          Regression test included as DISABLED_.
 *
 * Deferred: I/O, distribution, periodicity, thin shells, curves, order conversion.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "Mesh_Enums.hpp"

namespace
{
    const belfem::real tEps = 1e-12;
}

// =============================================================================
// §4.1  2D Tensor Mesh  [semantic]
// =============================================================================

TEST( TensorMesh2D, Quad4Construction )
{
    // order=1 → QUAD4, 3×3 nodes, step=1.0 → 2×2=4 elements, 9 nodes
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    EXPECT_EQ( tMesh.number_of_nodes(), 9u );
    EXPECT_EQ( tMesh.number_of_elements(), 4u );
    EXPECT_EQ( tMesh.number_of_dimensions(), 2u );
    EXPECT_TRUE( tMesh.is_tensormesh() );
}

TEST( TensorMesh2D, Quad4NodeCoords )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    // node(i,j) — grid order is x-major (i loops over x, j over y)
    // corner (0,0) should be at origin
    belfem::mesh::Node * tN00 = tMesh.node( 0, 0 );
    EXPECT_NEAR( tN00->x(), 0.0, tEps );
    EXPECT_NEAR( tN00->y(), 0.0, tEps );

    // corner (2,2) should be at (2*step, 2*step) = (2,2)
    belfem::mesh::Node * tN22 = tMesh.node( 2, 2 );
    EXPECT_NEAR( tN22->x(), 2.0, tEps );
    EXPECT_NEAR( tN22->y(), 2.0, tEps );
}

TEST( TensorMesh2D, Quad9Construction )
{
    // order=2 → QUAD9, 5×5 nodes, step=0.5 → 2×2=4 elements, 25 nodes
    belfem::Mesh tMesh( 2,
        belfem::Vector< belfem::index_t >( { 5, 5 } ),
        belfem::Vector< belfem::real >( { 0.5, 0.5 } ) );

    EXPECT_EQ( tMesh.number_of_nodes(), 25u );
    EXPECT_EQ( tMesh.number_of_elements(), 4u );
}

TEST( TensorMesh2D, Quad16Construction )
{
    // order=3 → QUAD16, 7×7 nodes, step=1.0 → 2×2=4 elements, 49 nodes
    belfem::Mesh tMesh( 3,
        belfem::Vector< belfem::index_t >( { 7, 7 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    EXPECT_EQ( tMesh.number_of_nodes(), 49u );
    EXPECT_EQ( tMesh.number_of_elements(), 4u );
}

TEST( TensorMesh2D, AllNodesHaveZeroZ )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    for( belfem::index_t k = 0; k < tMesh.number_of_nodes(); ++k )
    {
        EXPECT_NEAR( tMesh.nodes()( k )->z(), 0.0, tEps );
    }
}

// =============================================================================
// §4.2  3D Tensor Mesh  [semantic]
// =============================================================================

TEST( TensorMesh3D, Hex8Construction )
{
    // order=1 → HEX8, 3×3×3 nodes, step=1.0 → 2×2×2=8 elements, 27 nodes
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0, 1.0 } ) );

    EXPECT_EQ( tMesh.number_of_nodes(), 27u );
    EXPECT_EQ( tMesh.number_of_elements(), 8u );
    EXPECT_EQ( tMesh.number_of_dimensions(), 3u );
    EXPECT_TRUE( tMesh.is_tensormesh() );
}

TEST( TensorMesh3D, Hex27Construction )
{
    // order=2 → HEX27, 5×5×5 nodes, step=0.5 → 2×2×2=8 elements, 125 nodes
    belfem::Mesh tMesh( 2,
        belfem::Vector< belfem::index_t >( { 5, 5, 5 } ),
        belfem::Vector< belfem::real >( { 0.5, 0.5, 0.5 } ) );

    EXPECT_EQ( tMesh.number_of_nodes(), 125u );
    EXPECT_EQ( tMesh.number_of_elements(), 8u );
}

// BUG-ME1 regression: node(i,j,k) asserts dimension==2 instead of ==3
// Remove DISABLED_ prefix after fixing the assert in cl_Mesh.hpp line ~1992
TEST( TensorMesh3D, RegressionBugME1_NodeIJKAccessor )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0, 1.0 } ) );

    // After fix: node(0,0,0) should return the origin node
    belfem::mesh::Node * tN = tMesh.node( 0, 0, 0 );
    EXPECT_NEAR( tN->x(), 0.0, tEps );
    EXPECT_NEAR( tN->y(), 0.0, tEps );
    EXPECT_NEAR( tN->z(), 0.0, tEps );

    // far corner
    belfem::mesh::Node * tNfar = tMesh.node( 2, 2, 2 );
    EXPECT_NEAR( tNfar->x(), 2.0, tEps );
    EXPECT_NEAR( tNfar->y(), 2.0, tEps );
    EXPECT_NEAR( tNfar->z(), 2.0, tEps );
}

// =============================================================================
// §4.3  Mesh with Origin Offset  [semantic]
// =============================================================================

TEST( TensorMesh2D, WithOriginOffset )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ),
        belfem::Vector< belfem::real >( { 5.0, 10.0 } ) );

    belfem::mesh::Node * tN00 = tMesh.node( 0, 0 );
    EXPECT_NEAR( tN00->x(), 5.0, tEps );
    EXPECT_NEAR( tN00->y(), 10.0, tEps );

    belfem::mesh::Node * tN22 = tMesh.node( 2, 2 );
    EXPECT_NEAR( tN22->x(), 7.0, tEps );
    EXPECT_NEAR( tN22->y(), 12.0, tEps );
}

// =============================================================================
// §4.4  Mesh Accessors  [semantic]
// =============================================================================

TEST( MeshAccessors, NodeByID )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    // node IDs are 1-based
    belfem::mesh::Node * tN = tMesh.node( 1 );
    EXPECT_NE( tN, nullptr );
    EXPECT_EQ( tN->id(), 1u );
}

TEST( MeshAccessors, ElementByID )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    // element IDs are 1-based
    belfem::mesh::Element * tE = tMesh.element( 1 );
    EXPECT_NE( tE, nullptr );
    EXPECT_EQ( tE->id(), 1u );
}

TEST( MeshAccessors, BlockCountAtLeastOne )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    EXPECT_GE( tMesh.number_of_blocks(), 1u );
}

// =============================================================================
// §5  Map and Index Consistency  [semantic]
// =============================================================================

TEST( MeshMaps, NodeMapConsistent )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    for( belfem::index_t k = 0; k < tMesh.number_of_nodes(); ++k )
    {
        belfem::mesh::Node * tN = tMesh.nodes()( k );
        EXPECT_EQ( tMesh.node( tN->id() ), tN );
    }
}

TEST( MeshMaps, ElementMapConsistent )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    for( belfem::index_t k = 0; k < tMesh.number_of_elements(); ++k )
    {
        belfem::mesh::Element * tE = tMesh.elements()( k );
        EXPECT_EQ( tMesh.element( tE->id() ), tE );
    }
}

TEST( MeshMaps, BlockMapConsistent )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    for( belfem::uint k = 0; k < tMesh.number_of_blocks(); ++k )
    {
        belfem::mesh::Block * tB = tMesh.blocks()( k );
        EXPECT_EQ( tMesh.block( tB->id() ), tB );
    }
}

// =============================================================================
// §6  Connectivity on Tensor Mesh  [semantic]
// The tensor mesh constructor calls finalize() but does NOT set the
// Connectivity::Compute flag. Without it, finalize() skips node-to-element
// connectivity. We must explicitly enable it.
// (Caught by Gemini and ChatGPT during review)
// =============================================================================

TEST( MeshConnectivity, CornerNodeConnectivity )
{
    // 2D QUAD4 grid 4×4 nodes → 3×3=9 elements
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 4, 4 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    // enable connectivity and re-finalize
    tMesh.unfinalize();
    tMesh.set_connectivity( belfem::Connectivity::Compute );
    tMesh.finalize();

    // corner node (0,0) should be connected to exactly 1 element
    belfem::mesh::Node * tCorner = tMesh.node( 0, 0 );
    EXPECT_EQ( tCorner->number_of_elements(), 1u );
}

TEST( MeshConnectivity, InteriorNodeConnectivity )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 4, 4 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    tMesh.unfinalize();
    tMesh.set_connectivity( belfem::Connectivity::Compute );
    tMesh.finalize();

    // interior node (1,1) should be connected to 4 elements
    belfem::mesh::Node * tInterior = tMesh.node( 1, 1 );
    EXPECT_EQ( tInterior->number_of_elements(), 4u );
}

TEST( MeshConnectivity, EdgeNodeConnectivity )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 4, 4 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    tMesh.unfinalize();
    tMesh.set_connectivity( belfem::Connectivity::Compute );
    tMesh.finalize();

    // edge node (1,0) — not corner, on boundary → 2 elements
    belfem::mesh::Node * tEdge = tMesh.node( 1, 0 );
    EXPECT_EQ( tEdge->number_of_elements(), 2u );
}

TEST( MeshConnectivity, EdgesDoNotExistAfterConnectivity )
{
    // Edges are lazy — finalize() with connectivity does NOT create them
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 4, 4 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    tMesh.unfinalize();
    tMesh.set_connectivity( belfem::Connectivity::Compute );
    tMesh.finalize();

    EXPECT_FALSE( tMesh.edges_exist() );
}

TEST( MeshEdges, EdgesExistAfterCreateEdges )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 4, 4 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    tMesh.unfinalize();
    tMesh.set_connectivity( belfem::Connectivity::Compute );
    tMesh.finalize();

    EXPECT_FALSE( tMesh.edges_exist() );

    tMesh.create_edges();

    EXPECT_TRUE( tMesh.edges_exist() );
    EXPECT_GT( tMesh.number_of_edges(), 0u );
}

// =============================================================================
// §7  Scale Mesh  [semantic]
// =============================================================================

TEST( MeshScale, ScaleMeshFactor )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    // verify pre-scale: far corner at (2,2)
    belfem::mesh::Node * tFar = tMesh.node( 2, 2 );
    EXPECT_NEAR( tFar->x(), 2.0, tEps );

    tMesh.scale_mesh( 0.001 );

    EXPECT_NEAR( tFar->x(), 0.002, tEps );
    EXPECT_NEAR( tFar->y(), 0.002, tEps );
}

TEST( MeshScale, ScaleMeshPreservesTopology )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    belfem::index_t tNodesBefore = tMesh.number_of_nodes();
    belfem::index_t tElemsBefore = tMesh.number_of_elements();

    tMesh.scale_mesh( 0.001 );

    EXPECT_EQ( tMesh.number_of_nodes(), tNodesBefore );
    EXPECT_EQ( tMesh.number_of_elements(), tElemsBefore );
}

// =============================================================================
// §8  Block Basics  [semantic]
// =============================================================================

TEST( MeshBlock, BlockElementType )
{
    // (idea from Gemini) — verify block reports correct element type
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    belfem::mesh::Block * tBlock = tMesh.blocks()( 0 );
    EXPECT_EQ( tBlock->element_type(), belfem::ElementType::QUAD4 );
}

TEST( MeshBlock, BlockElementCount )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    // single-block tensor mesh: block element count == mesh element count
    belfem::mesh::Block * tBlock = tMesh.blocks()( 0 );
    EXPECT_EQ( tBlock->number_of_elements(),
               static_cast< belfem::index_t >( tMesh.number_of_elements() ) );
}

TEST( MeshBlock, BlockElementAccess )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    belfem::mesh::Block * tBlock = tMesh.blocks()( 0 );
    belfem::mesh::Element * tElem = tBlock->element( 0 );
    EXPECT_NE( tElem, nullptr );
}

// =============================================================================
// §8b  Edges On Demand  [semantic]
// (idea from ChatGPT) — edges are lazy; create_edges() materializes them
// =============================================================================

TEST( MeshEdges, EdgesDoNotExistInitially )
{
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    EXPECT_FALSE( tMesh.edges_exist() );
    EXPECT_EQ( tMesh.number_of_edges(), 0u );
}

TEST( MeshEdges, EdgesCreatedOnDemand )
{
    // (idea from ChatGPT)
    belfem::Mesh tMesh( 1,
        belfem::Vector< belfem::index_t >( { 3, 3 } ),
        belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

    tMesh.create_edges();

    EXPECT_TRUE( tMesh.edges_exist() );
    // 3×3 QUAD4 grid: 2 horizontal edges per row × 3 rows
    //               + 2 vertical edges per col × 3 cols = 12
    EXPECT_EQ( tMesh.number_of_edges(), 12u );
}

// =============================================================================
// §9  Memory and Cleanup  [valgrind]
// =============================================================================

TEST( MeshMemory, MeshDestructorCleansUp )
{
    // (idea from ChatGPT) — Valgrind target
    {
        belfem::Mesh tMesh( 1,
            belfem::Vector< belfem::index_t >( { 8, 8 } ),
            belfem::Vector< belfem::real >( { 1.0, 1.0 } ) );

        EXPECT_EQ( tMesh.number_of_nodes(), 64u );
    }
    // mesh goes out of scope — Valgrind should report zero leaks
    SUCCEED();
}
