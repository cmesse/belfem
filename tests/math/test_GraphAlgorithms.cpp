/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for graph algorithms: BFS, DFS, find_connected_partitions,
 * find_pseudo_peripheral_vertex/node, symrcm.
 * See: tests_07_graph.md §2–§6
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <set>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "fn_Graph_bfs.hpp"
#include "fn_Graph_dfs.hpp"
#include "fn_Graph_clear.hpp"
#include "fn_Graph_find_connected_partitions.hpp"
#include "fn_Graph_find_pseudo_peripheral_node.hpp"
#include "fn_Graph_find_pseudo_peripheral_vertex.hpp"
#include "fn_Graph_symrcm.hpp"

// =============================================================================
// Test Helpers
// =============================================================================

namespace
{
    // Build a graph from symmetric adjacency lists.
    // aAdjacency(i) is a Cell<index_t> of neighbor indices for vertex i.
    belfem::Graph
    build_test_graph( const belfem::Cell< belfem::Cell< belfem::index_t > > & aAdjacency )
    {
        belfem::index_t tN = aAdjacency.size();
        belfem::Graph tGraph( tN, nullptr );

        for( belfem::index_t i = 0; i < tN; ++i )
        {
            tGraph( i ) = new belfem::graph::Vertex();
            tGraph( i )->set_id( i );
            tGraph( i )->set_index( i );
        }

        // count neighbors
        for( belfem::index_t i = 0; i < tN; ++i )
        {
            for( belfem::index_t j = 0; j < aAdjacency( i ).size(); ++j )
            {
                tGraph( i )->increment_vertex_counter();
            }
        }

        // allocate and fill
        for( belfem::index_t i = 0; i < tN; ++i )
        {
            tGraph( i )->init_vertex_container();
            for( belfem::index_t j = 0; j < aAdjacency( i ).size(); ++j )
            {
                tGraph( i )->insert_vertex( tGraph( aAdjacency( i )( j ) ) );
            }
        }

        return tGraph;
    }

    // Path₅:  0—1—2—3—4
    belfem::Graph make_path5()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 5, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 1 };
        tAdj( 1 ) = { 0, 2 };
        tAdj( 2 ) = { 1, 3 };
        tAdj( 3 ) = { 2, 4 };
        tAdj( 4 ) = { 3 };
        return build_test_graph( tAdj );
    }

    // Star₅:  0 connected to {1,2,3,4}
    belfem::Graph make_star5()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 5, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 1, 2, 3, 4 };
        tAdj( 1 ) = { 0 };
        tAdj( 2 ) = { 0 };
        tAdj( 3 ) = { 0 };
        tAdj( 4 ) = { 0 };
        return build_test_graph( tAdj );
    }

    // K₄: complete graph on 4 vertices
    belfem::Graph make_k4()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 4, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 1, 2, 3 };
        tAdj( 1 ) = { 0, 2, 3 };
        tAdj( 2 ) = { 0, 1, 3 };
        tAdj( 3 ) = { 0, 1, 2 };
        return build_test_graph( tAdj );
    }

    // BinaryTree₇:  0→{1,2}, 1→{3,4}, 2→{5,6} (symmetric)
    belfem::Graph make_binary_tree7()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 7, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 1, 2 };
        tAdj( 1 ) = { 0, 3, 4 };
        tAdj( 2 ) = { 0, 5, 6 };
        tAdj( 3 ) = { 1 };
        tAdj( 4 ) = { 1 };
        tAdj( 5 ) = { 2 };
        tAdj( 6 ) = { 2 };
        return build_test_graph( tAdj );
    }

    // Disconnected: two triangles {0,1,2} and {3,4,5}
    belfem::Graph make_disconnected()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 6, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 1, 2 };
        tAdj( 1 ) = { 0, 2 };
        tAdj( 2 ) = { 0, 1 };
        tAdj( 3 ) = { 4, 5 };
        tAdj( 4 ) = { 3, 5 };
        tAdj( 5 ) = { 3, 4 };
        return build_test_graph( tAdj );
    }

    // Single vertex, no edges
    belfem::Graph make_single()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 1, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = {};
        return build_test_graph( tAdj );
    }

    // Compute bandwidth: max |index(v) - index(neighbor)| over all edges
    belfem::index_t compute_bandwidth( const belfem::Graph & aGraph )
    {
        belfem::index_t tBw = 0;
        for( belfem::index_t i = 0; i < aGraph.size(); ++i )
        {
            for( belfem::uint k = 0; k < aGraph( i )->number_of_vertices(); ++k )
            {
                belfem::index_t tDiff = std::abs(
                    ( long ) aGraph( i )->index()
                  - ( long ) aGraph( i )->vertex( k )->index() );
                tBw = std::max( tBw, tDiff );
            }
        }
        return tBw;
    }

    // Collect the set of neighbor IDs for a vertex
    std::set< belfem::id_t > neighbor_ids( belfem::graph::Vertex * aV )
    {
        std::set< belfem::id_t > tIds;
        for( belfem::uint k = 0; k < aV->number_of_vertices(); ++k )
        {
            tIds.insert( aV->vertex( k )->id() );
        }
        return tIds;
    }
}

// =============================================================================
// §2  BFS  [semantic]
// =============================================================================

TEST( GraphBfs, BfsPath5FromEnd )
{
    belfem::Graph tG = make_path5();

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tG( 0 ) );

    // path: every level has 1 vertex → max width = 1
    EXPECT_EQ( tMaxWidth, 1u );

    // verify levels: 0,1,2,3,4 from vertex with id=0
    // NOTE: BFS reindexes vertices, so find by id
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        if( tG( i )->id() == 0 )
        {
            EXPECT_EQ( tG( i )->level(), 0u );
        }
        if( tG( i )->id() == 4 )
        {
            EXPECT_EQ( tG( i )->level(), 4u );
        }
    }

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsPath5FromMiddle )
{
    belfem::Graph tG = make_path5();

    // find vertex with id=2
    belfem::graph::Vertex * tStart = nullptr;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        if( tG( i )->id() == 2 ) tStart = tG( i );
    }

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tStart );

    // levels from middle: 2,1,0,1,2 → max width = 2 (two vertices at level 2)
    EXPECT_EQ( tMaxWidth, 2u );

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsStar5FromCenter )
{
    belfem::Graph tG = make_star5();

    // vertex id=0 is the center
    belfem::graph::Vertex * tCenter = nullptr;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        if( tG( i )->id() == 0 ) tCenter = tG( i );
    }

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tCenter );

    // all 4 leaves at level 1 → max width = 4
    EXPECT_EQ( tMaxWidth, 4u );

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsStar5FromLeaf )
{
    belfem::Graph tG = make_star5();

    belfem::graph::Vertex * tLeaf = nullptr;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        if( tG( i )->id() == 1 ) tLeaf = tG( i );
    }

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tLeaf );

    // leaf(0) → center(1) → 3 other leaves(2) → max width = 3
    EXPECT_EQ( tMaxWidth, 3u );

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsK4 )
{
    belfem::Graph tG = make_k4();

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tG( 0 ) );

    // all 3 others at level 1 → max width = 3
    EXPECT_EQ( tMaxWidth, 3u );

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsBinaryTree )
{
    belfem::Graph tG = make_binary_tree7();

    belfem::graph::Vertex * tRoot = nullptr;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        if( tG( i )->id() == 0 ) tRoot = tG( i );
    }

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tRoot );

    // level 0: root (1), level 1: 2 children (2), level 2: 4 leaves (4)
    // max width = 4
    EXPECT_EQ( tMaxWidth, 4u );

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsSingleVertex )
{
    belfem::Graph tG = make_single();

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, tG( 0 ) );

    EXPECT_EQ( tMaxWidth, 1u );
    EXPECT_EQ( tG( 0 )->level(), 0u );

    belfem::graph::clear( tG );
}

TEST( GraphBfs, BfsEmptyGraph )
{
    belfem::Graph tG;
    belfem::index_t tMaxWidth = belfem::graph::bfs( tG, nullptr );
    EXPECT_EQ( tMaxWidth, 0u );
}

// --- §2.2 BFS disconnected overload ---

TEST( GraphBfs, BfsDisconnectedGraph )
{
    belfem::Graph tG = make_disconnected();

    belfem::index_t tMaxWidth = belfem::graph::bfs( tG );

    // each triangle: level 0(1), level 1(2) → max width per component = 2
    EXPECT_EQ( tMaxWidth, 2u );

    // all vertices should have levels set
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        EXPECT_NE( tG( i )->level(), belfem::gNoIndex );
    }

    belfem::graph::clear( tG );
}

// =============================================================================
// §3  DFS  [semantic]
// =============================================================================

TEST( GraphDfs, DfsConnectedReturnsOne )
{
    belfem::Graph tG = make_path5();

    belfem::proc_t tComponents = belfem::graph::dfs( tG );

    EXPECT_EQ( tComponents, 1 );

    // all vertices should have the same owner
    belfem::proc_t tOwner = tG( 0 )->owner();
    for( belfem::index_t i = 1; i < tG.size(); ++i )
    {
        EXPECT_EQ( tG( i )->owner(), tOwner );
    }

    belfem::graph::clear( tG );
}

TEST( GraphDfs, DfsDisconnectedReturnsTwoComponents )
{
    belfem::Graph tG = make_disconnected();

    belfem::proc_t tComponents = belfem::graph::dfs( tG );

    EXPECT_EQ( tComponents, 2 );

    // vertices within each triangle share owner
    // find owners by id
    std::set< belfem::proc_t > tOwners;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tOwners.insert( tG( i )->owner() );
    }
    EXPECT_EQ( tOwners.size(), 2u );

    belfem::graph::clear( tG );
}

TEST( GraphDfs, DfsK4ReturnsOne )
{
    belfem::Graph tG = make_k4();
    EXPECT_EQ( belfem::graph::dfs( tG ), 1 );
    belfem::graph::clear( tG );
}

TEST( GraphDfs, DfsSingleVertex )
{
    belfem::Graph tG = make_single();
    EXPECT_EQ( belfem::graph::dfs( tG ), 1 );
    belfem::graph::clear( tG );
}

TEST( GraphDfs, DfsEmptyGraph )
{
    belfem::Graph tG;
    EXPECT_EQ( belfem::graph::dfs( tG ), 0 );
}

TEST( GraphDfs, DfsAllVerticesFlagged )
{
    belfem::Graph tG = make_path5();
    belfem::graph::dfs( tG );

    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        EXPECT_TRUE( tG( i )->is_flagged() );
    }

    belfem::graph::clear( tG );
}

TEST( GraphDfs, DfsLevelsSet )
{
    // After DFS, all vertices in a connected graph should have level != gNoIndex
    belfem::Graph tG = make_path5();
    belfem::graph::dfs( tG );

    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        EXPECT_NE( tG( i )->level(), belfem::gNoIndex );
    }

    belfem::graph::clear( tG );
}

// =============================================================================
// §4  Connected Partitions  [semantic]
// =============================================================================

TEST( GraphPartitions, FindConnectedPartitionsOneComponent )
{
    belfem::Graph tG = make_path5();

    belfem::index_t tLargestSize = belfem::graph::find_connected_partitions( tG );

    // 1 partition of 5 vertices → returns 5 (size of largest)
    EXPECT_EQ( tLargestSize, 5u );

    // verify all vertices share the same owner (1 partition)
    std::set< belfem::proc_t > tOwners;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tOwners.insert( tG( i )->owner() );
    }
    EXPECT_EQ( tOwners.size(), 1u );

    belfem::graph::clear( tG );
}

TEST( GraphPartitions, FindConnectedPartitionsTwoComponents )
{
    belfem::Graph tG = make_disconnected();

    belfem::index_t tLargestSize = belfem::graph::find_connected_partitions( tG );

    // 2 partitions of 3 each → largest = 3
    EXPECT_EQ( tLargestSize, 3u );

    std::set< belfem::proc_t > tOwners;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tOwners.insert( tG( i )->owner() );
    }
    EXPECT_EQ( tOwners.size(), 2u );

    belfem::graph::clear( tG );
}

TEST( GraphPartitions, FindConnectedPartitionsThreeIsolated )
{
    // 3 isolated vertices → 3 partitions of size 1
    belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 3, belfem::Cell< belfem::index_t >() );
    tAdj( 0 ) = {};
    tAdj( 1 ) = {};
    tAdj( 2 ) = {};
    belfem::Graph tG = build_test_graph( tAdj );

    belfem::index_t tLargestSize = belfem::graph::find_connected_partitions( tG );

    EXPECT_EQ( tLargestSize, 1u );

    std::set< belfem::proc_t > tOwners;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tOwners.insert( tG( i )->owner() );
    }
    EXPECT_EQ( tOwners.size(), 3u );

    belfem::graph::clear( tG );
}

TEST( GraphPartitions, FindConnectedPartitionsEmpty )
{
    // Source was fixed: early return for empty graph
    belfem::Graph tG;
    belfem::index_t tResult = belfem::graph::find_connected_partitions( tG );
    EXPECT_EQ( tResult, 0u );
}

// =============================================================================
// §5  Pseudo-Peripheral Vertex  [semantic]
// =============================================================================

TEST( GraphPseudoPeripheral, PseudoPeripheralPath5 )
{
    belfem::Graph tG = make_path5();

    belfem::graph::Vertex * tV =
        belfem::graph::find_pseudo_peripheral_vertex( tG );

    // should return an endpoint (id 0 or 4)
    ASSERT_NE( tV, nullptr );
    EXPECT_TRUE( tV->id() == 0 || tV->id() == 4 );

    belfem::graph::clear( tG );
}

TEST( GraphPseudoPeripheral, PseudoPeripheralBinaryTree )
{
    belfem::Graph tG = make_binary_tree7();

    belfem::graph::Vertex * tV =
        belfem::graph::find_pseudo_peripheral_vertex( tG );

    ASSERT_NE( tV, nullptr );
    // should be a leaf (id 3,4,5,6) — max eccentricity
    EXPECT_TRUE( tV->id() >= 3 && tV->id() <= 6 );

    belfem::graph::clear( tG );
}

TEST( GraphPseudoPeripheral, PseudoPeripheralK4 )
{
    belfem::Graph tG = make_k4();

    belfem::graph::Vertex * tV =
        belfem::graph::find_pseudo_peripheral_vertex( tG );

    ASSERT_NE( tV, nullptr );
    // all equivalent — just verify no crash

    belfem::graph::clear( tG );
}

TEST( GraphPseudoPeripheral, PseudoPeripheralSingleVertex )
{
    belfem::Graph tG = make_single();

    belfem::graph::Vertex * tV =
        belfem::graph::find_pseudo_peripheral_vertex( tG );

    ASSERT_NE( tV, nullptr );

    belfem::graph::clear( tG );
}

TEST( GraphPseudoPeripheral, PseudoPeripheralEmptyGraph )
{
    belfem::Graph tG;

    belfem::graph::Vertex * tV =
        belfem::graph::find_pseudo_peripheral_vertex( tG );

    EXPECT_EQ( tV, nullptr );
}

TEST( GraphPseudoPeripheral, PseudoPeripheralNodePath5 )
{
    // BUG-G2 regression: find_pseudo_peripheral_node() confused BFS width
    // with depth. On Path₅, it returned a near neighbor instead of an endpoint.
    // After fix, it should return an endpoint (id 0 or 4, degree 1).
    belfem::Graph tG = make_path5();

    belfem::graph::Vertex * tV =
        belfem::graph::find_pseudo_peripheral_node( tG, tG( 0 ) );

    ASSERT_NE( tV, nullptr );

    // on Path₅, the pseudo-peripheral node must be an endpoint (degree 1)
    EXPECT_EQ( tV->number_of_vertices(), 1u );
    // and its id must be 0 or 4 (the two endpoints)
    EXPECT_TRUE( tV->id() == 0 || tV->id() == 4 );

    belfem::graph::clear( tG );
}

// =============================================================================
// §6  SymRCM  [semantic]
// =============================================================================

TEST( GraphSymrcm, SymrcmPath5 )
{
    belfem::Graph tG = make_path5();

    // save neighbor IDs before reordering
    std::vector< std::set< belfem::id_t > > tOrigAdj( tG.size() );
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tOrigAdj[ tG( i )->id() ] = neighbor_ids( tG( i ) );
    }

    belfem::index_t tBwBefore = compute_bandwidth( tG );

    belfem::graph::symrcm( tG );

    belfem::index_t tBwAfter = compute_bandwidth( tG );

    // bandwidth should not increase
    EXPECT_LE( tBwAfter, tBwBefore );

    // verify permutation validity: all indices 0..N-1 appear exactly once
    std::set< belfem::index_t > tIndices;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tIndices.insert( tG( i )->index() );
    }
    EXPECT_EQ( tIndices.size(), tG.size() );

    // verify adjacency is preserved after reordering (Codex finding)
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        EXPECT_EQ( neighbor_ids( tG( i ) ), tOrigAdj[ tG( i )->id() ] );
    }

    belfem::graph::clear( tG );
}

TEST( GraphSymrcm, SymrcmStar5 )
{
    belfem::Graph tG = make_star5();

    belfem::index_t tBwBefore = compute_bandwidth( tG );

    belfem::graph::symrcm( tG );

    belfem::index_t tBwAfter = compute_bandwidth( tG );

    // star graph has high bandwidth (center connects to all); symrcm should reduce it
    EXPECT_LE( tBwAfter, tBwBefore );

    // verify permutation validity
    std::set< belfem::index_t > tIndices;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tIndices.insert( tG( i )->index() );
    }
    EXPECT_EQ( tIndices.size(), tG.size() );

    belfem::graph::clear( tG );
}

TEST( GraphSymrcm, SymrcmBinaryTree )
{
    belfem::Graph tG = make_binary_tree7();

    belfem::index_t tBwBefore = compute_bandwidth( tG );

    belfem::graph::symrcm( tG );

    belfem::index_t tBwAfter = compute_bandwidth( tG );

    EXPECT_LE( tBwAfter, tBwBefore );

    belfem::graph::clear( tG );
}

TEST( GraphSymrcm, SymrcmDisconnected )
{
    belfem::Graph tG = make_disconnected();

    // should handle both components without crash
    belfem::graph::symrcm( tG );

    // verify all indices unique
    std::set< belfem::index_t > tIndices;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tIndices.insert( tG( i )->index() );
    }
    EXPECT_EQ( tIndices.size(), tG.size() );

    belfem::graph::clear( tG );
}

TEST( GraphSymrcm, SymrcmSingleVertex )
{
    belfem::Graph tG = make_single();
    belfem::graph::symrcm( tG );
    EXPECT_EQ( tG( 0 )->index(), 0u );
    belfem::graph::clear( tG );
}

TEST( GraphSymrcm, SymrcmEmptyGraph )
{
    belfem::Graph tG;
    belfem::graph::symrcm( tG );
    EXPECT_EQ( tG.size(), 0u );
}

TEST( GraphSymrcm, SymrcmPermutationIsValid )
{
    belfem::Graph tG = make_star5();

    belfem::graph::symrcm( tG );

    std::set< belfem::index_t > tIndices;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tIndices.insert( tG( i )->index() );
    }
    // every index 0..4 appears exactly once
    EXPECT_EQ( tIndices.size(), 5u );
    EXPECT_EQ( *tIndices.begin(), 0u );
    EXPECT_EQ( *tIndices.rbegin(), 4u );

    belfem::graph::clear( tG );
}

// =============================================================================
// §9  Memory Management  [semantic]
// =============================================================================

TEST( GraphMemory, GraphClearDeletesVertices )
{
    belfem::Graph tG = make_path5();

    belfem::graph::clear( tG );

    EXPECT_EQ( tG.size(), 0u );
}
