/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for graph utilities: CSR export, graph::sort (BUG-G1),
 * apply_graph_permutation, comparison operators.
 * See: tests_07_graph.md §7
 */

#include <gtest/gtest.h>
#include <set>
#include <algorithm>
#include <stdexcept>

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Vector.hpp"
#include "fn_Graph_clear.hpp"
#include "fn_Graph_sort.hpp"
#include "graphtools.hpp"
#include "commtools.hpp"
#include "op_Graph_Vertex_Degree.hpp"
#include "op_Graph_Vertex_ID.hpp"
#include "op_Graph_Vertex_Index.hpp"
#include "op_Graph_Vertex_Level.hpp"
#include "op_Graph_Vertex_Owner.hpp"

#ifdef BELFEM_METIS
#include "fn_Graph_METIS.hpp"
#endif

// =============================================================================
// Test Helpers (same builder as test_GraphAlgorithms.cpp)
// =============================================================================

namespace
{
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

        for( belfem::index_t i = 0; i < tN; ++i )
        {
            for( belfem::index_t j = 0; j < aAdjacency( i ).size(); ++j )
            {
                tGraph( i )->increment_vertex_counter();
            }
        }

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

    // path5 with a self-loop on every vertex: the shape of a graph built from
    // a sparse matrix, whose diagonal becomes a loop
    belfem::Graph make_path5_with_self_loops()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 5, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 0, 1 };
        tAdj( 1 ) = { 0, 1, 2 };
        tAdj( 2 ) = { 1, 2, 3 };
        tAdj( 3 ) = { 2, 3, 4 };
        tAdj( 4 ) = { 3, 4 };
        return build_test_graph( tAdj );
    }

    // three isolated vertices, each with only a self-loop: after the loop drop
    // the CSR has vertices but no edges
    belfem::Graph make_self_loops_only()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 3, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 0 };
        tAdj( 1 ) = { 1 };
        tAdj( 2 ) = { 2 };
        return build_test_graph( tAdj );
    }

    belfem::Graph make_k4()
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 4, belfem::Cell< belfem::index_t >() );
        tAdj( 0 ) = { 1, 2, 3 };
        tAdj( 1 ) = { 0, 2, 3 };
        tAdj( 2 ) = { 0, 1, 3 };
        tAdj( 3 ) = { 0, 1, 2 };
        return build_test_graph( tAdj );
    }
}

// =============================================================================
// §7.1  CSR Export  [semantic]
// =============================================================================

TEST( GraphTools, BuildAdjacencyPath5 )
{
    belfem::Graph tG = make_path5();
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    // 5 vertices → tVertices has 6 entries (N+1)
    EXPECT_EQ( tVertices.length(), 6u );

    // total edges: 4 undirected edges × 2 directions = 8
    EXPECT_EQ( tVertices( 5 ), 8 );

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildAdjacencyK4 )
{
    belfem::Graph tG = make_k4();
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    // 4 vertices → tVertices has 5 entries
    EXPECT_EQ( tVertices.length(), 5u );

    // 6 undirected edges × 2 = 12
    EXPECT_EQ( tVertices( 4 ), 12 );

    // every vertex has 3 neighbors
    for( int i = 0; i < 4; ++i )
    {
        EXPECT_EQ( tVertices( i + 1 ) - tVertices( i ), 3 );
    }

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildAdjacencyEmpty )
{
    belfem::Graph tG;
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    // N=0 → tVertices has 1 entry (value 0)
    EXPECT_EQ( tVertices.length(), 1u );
    EXPECT_EQ( tVertices( 0 ), 0 );
}

TEST( GraphTools, BuildAdjacencySingleVertex )
{
    belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( 1, belfem::Cell< belfem::index_t >() );
    tAdj( 0 ) = {};
    belfem::Graph tG = build_test_graph( tAdj );

    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    EXPECT_EQ( tVertices.length(), 2u );
    EXPECT_EQ( tVertices( 0 ), 0 );
    EXPECT_EQ( tVertices( 1 ), 0 );

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildAdjacencyVerticesMonotone )
{
    belfem::Graph tG = make_path5();
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    for( belfem::index_t i = 1; i < tVertices.length(); ++i )
    {
        EXPECT_LE( tVertices( i - 1 ), tVertices( i ) );
    }

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildAdjacencyDegreeConsistent )
{
    belfem::Graph tG = make_path5();
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        int tDegree = tVertices( i + 1 ) - tVertices( i );
        EXPECT_EQ( static_cast< belfem::uint >( tDegree ),
                   tG( i )->number_of_vertices() );
    }

    belfem::graph::clear( tG );
}

// self-loops never reach the CSR the graph libraries see: METIS_NodeNDP
// corrupts its heap on them ( reproduced 2026-09-04 on the DistMatrix graph )
TEST( GraphTools, BuildAdjacencyDropsSelfLoops )
{
    belfem::Graph tG = make_path5_with_self_loops();
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    ASSERT_EQ( tVertices.length(), 6u );
    EXPECT_EQ( tVertices( 5 ), 8 );            // the 5 loops are gone, 8 directed entries stay
    EXPECT_EQ( tEdges.length(), 8u );
    for ( int i = 0; i < 5; ++i )
    {
        for ( int e = tVertices( i ); e < tVertices( i + 1 ); ++e )
        {
            EXPECT_NE( tEdges( e ), i ) << "row " << i;
        }
    }

    // the Graph itself must KEEP its loops: DistMatrix builds the permuted
    // SpMatrix from the Graph and needs the diagonal
    for ( belfem::index_t i = 0; i < 5; ++i )
    {
        bool tHasLoop = false;
        for ( belfem::uint k = 0; k < tG( i )->number_of_vertices(); ++k )
        {
            if ( tG( i )->vertex( k ) == tG( i ) ) tHasLoop = true;
        }
        EXPECT_TRUE( tHasLoop ) << "vertex " << i << " lost its self-loop";
    }

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildParAdjacencyDropsSelfLoops )
{
    belfem::Graph tG = make_path5_with_self_loops();
    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    ASSERT_EQ( tVertices( 0 ).length(), 6u );
    EXPECT_EQ( tVertices( 0 )( 5 ), 8 );
    EXPECT_EQ( tEdges( 0 ).length(), 8u );
    for ( int i = 0; i < 5; ++i )
    {
        for ( int e = tVertices( 0 )( i ); e < tVertices( 0 )( i + 1 ); ++e )
        {
            EXPECT_NE( tEdges( 0 )( e ), i ) << "row " << i;
        }
    }

    belfem::graph::clear( tG );
}

#ifdef BELFEM_METIS
// the crash itself: METIS_NodeNDP on a graph with self-loops corrupted its heap
// ( munmap_chunk(): invalid pointer ) before the builders dropped them. This is
// the serial gate for D5; the matrix-graph reproduction lives in the devlog
TEST( GraphTools, MetisNdpOrdersGraphWithSelfLoops )
{
    belfem::Graph tG = make_path5_with_self_loops();

    belfem::graph::metis_ndp( tG, 2 );

    belfem::Cell< belfem::uint > tSeen( 5, 0 );
    for ( belfem::graph::Vertex * tV : tG )
    {
        ASSERT_LT( tV->index(), 5u );
        ++tSeen( tV->index() );
    }
    for ( belfem::index_t k = 0; k < 5; ++k )
    {
        EXPECT_EQ( tSeen( k ), 1u ) << "index " << k;
    }

    belfem::graph::clear( tG );
}
#endif

// vertices but no edges once the loops are gone: the CSR terminal is 0 and the
// edge buffer keeps its one-element placeholder ( the libraries take a raw
// pointer, and an empty Vector's data() may be null under Armadillo )
TEST( GraphTools, BuildAdjacencySelfLoopsOnlyKeepsPlaceholder )
{
    belfem::Graph tG = make_self_loops_only();
    belfem::Vector< int > tVertices;
    belfem::Vector< int > tEdges;

    belfem::graph::build_graph_adjacency( tG, tVertices, tEdges );

    ASSERT_EQ( tVertices.length(), 4u );
    EXPECT_EQ( tVertices( 3 ), 0 );
    EXPECT_EQ( tEdges.length(), 1u );

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildParAdjacencySelfLoopsOnlyKeepsPlaceholder )
{
    belfem::Graph tG = make_self_loops_only();
    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    ASSERT_EQ( tVertices( 0 ).length(), 4u );   // 3 vertices on root, N+1 offsets
    EXPECT_EQ( tVertices( 0 )( 3 ), 0 );
    EXPECT_EQ( tEdges( 0 ).length(), 1u );

    belfem::graph::clear( tG );
}

// -----------------------------------------------------------------------------
// §7.1b  Parallel CSR export: owner sentinels  [semantic]
// build_pargraph_adjacency() indexes per-rank counters by vertex owner. The
// fixture never calls set_owner(), so every vertex carries the gNoOwner default:
// this is the path that once wrote out of bounds in release builds.
// -----------------------------------------------------------------------------

TEST( GraphTools, BuildParAdjacencyUnownedGoesToRoot )
{
    belfem::Graph tG = make_path5();
    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    const belfem::proc_t tSize = belfem::comm_size();
    ASSERT_EQ( tDistribution.length(), static_cast< size_t >( tSize + 1 ) );
    EXPECT_EQ( tDistribution( 0 ), 0 );
    EXPECT_EQ( tDistribution( tSize ), 5 );

    // all five vertices land on root: N+1 offsets, 4 undirected edges × 2
    ASSERT_EQ( tVertices( 0 ).length(), 6u );
    EXPECT_EQ( tVertices( 0 )( 5 ), 8 );

    for ( belfem::graph::Vertex * tV : tG )
    {
        EXPECT_EQ( tV->owner(), 0 );
    }

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildParAdjacencyRejectsNegativeOwner )
{
    belfem::Graph tG = make_path5();
    tG( 2 )->set_owner( -1 );

    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    EXPECT_THROW( belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges ),
                  std::runtime_error );

    // the throw happens in the first pass: vertices before the bad one are
    // already unflagged and remapped, but the reordered graph has not been
    // built, so tG still owns every vertex
    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildParAdjacencyCommSizeMarkerGoesToRoot )
{
    // the Kernel marks unassigned entities with comm_size() itself; that
    // sentinel must keep mapping to root exactly like the gNoOwner default
    belfem::Graph tG = make_path5();
    for ( belfem::graph::Vertex * tV : tG )
    {
        tV->set_owner( belfem::comm_size() );
    }

    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    const belfem::proc_t tSize = belfem::comm_size();
    ASSERT_EQ( tDistribution.length(), static_cast< size_t >( tSize + 1 ) );
    EXPECT_EQ( tDistribution( tSize ), 5 );
    ASSERT_EQ( tVertices( 0 ).length(), 6u );

    for ( belfem::graph::Vertex * tV : tG )
    {
        EXPECT_EQ( tV->owner(), 0 );
    }

    belfem::graph::clear( tG );
}

TEST( GraphTools, BuildParAdjacencyRejectsOwnerAboveCommSize )
{
    belfem::Graph tG = make_path5();

    // neither sentinel: comm_size() itself is the Kernel marker, gNoOwner the default
    const belfem::proc_t tBad = belfem::comm_size() + 1;
    ASSERT_NE( tBad, belfem::gNoOwner );
    tG( 2 )->set_owner( tBad );

    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    EXPECT_THROW( belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges ),
                  std::runtime_error );

    belfem::graph::clear( tG );
}

// =============================================================================
// §7.2  Graph Sort — BUG-G1 Regression  [semantic]
// =============================================================================

TEST( GraphTools, GraphSortRegressionBugG1 )
{
    // BUG-G1 regression: graph::sort() was a no-op (begin,begin → empty range).
    // After fix, graph should be sorted by vertex index.
    belfem::Graph tG = make_path5();

    // shuffle indices: assign 4,2,0,3,1
    tG( 0 )->set_index( 4 );
    tG( 1 )->set_index( 2 );
    tG( 2 )->set_index( 0 );
    tG( 3 )->set_index( 3 );
    tG( 4 )->set_index( 1 );

    belfem::graph::sort( tG );

    // after fix: graph should be sorted by index
    EXPECT_EQ( tG( 0 )->index(), 0u );
    EXPECT_EQ( tG( 1 )->index(), 1u );
    EXPECT_EQ( tG( 2 )->index(), 2u );
    EXPECT_EQ( tG( 3 )->index(), 3u );
    EXPECT_EQ( tG( 4 )->index(), 4u );

    belfem::graph::clear( tG );
}

// =============================================================================
// §7.3  Apply Graph Permutation  [semantic]
// =============================================================================

TEST( GraphTools, ApplyPermutationIdentity )
{
    belfem::Graph tG = make_path5();

    belfem::Vector< belfem::index_t > tPerm( 5 );
    for( belfem::index_t i = 0; i < 5; ++i )
    {
        tPerm( i ) = i;
    }

    belfem::graph::apply_graph_permutation( tG, tPerm );

    // identity permutation → indices unchanged
    for( belfem::index_t i = 0; i < 5; ++i )
    {
        EXPECT_EQ( tG( i )->index(), i );
    }

    belfem::graph::clear( tG );
}

TEST( GraphTools, ApplyPermutationReverse )
{
    belfem::Graph tG = make_path5();

    // save original neighbor IDs
    std::vector< std::set< belfem::id_t > > tOrigAdj( 5 );
    for( belfem::index_t i = 0; i < 5; ++i )
    {
        tOrigAdj[ tG( i )->id() ] = std::set< belfem::id_t >();
        for( belfem::uint k = 0; k < tG( i )->number_of_vertices(); ++k )
        {
            tOrigAdj[ tG( i )->id() ].insert( tG( i )->vertex( k )->id() );
        }
    }

    // reverse permutation: {4, 3, 2, 1, 0}
    belfem::Vector< belfem::index_t > tPerm = { 4, 3, 2, 1, 0 };

    belfem::graph::apply_graph_permutation( tG, tPerm );

    // all indices 0..4 should still appear exactly once
    std::set< belfem::index_t > tIndices;
    for( belfem::index_t i = 0; i < 5; ++i )
    {
        tIndices.insert( tG( i )->index() );
    }
    EXPECT_EQ( tIndices.size(), 5u );

    // verify adjacency is preserved (Codex finding: was missing)
    for( belfem::index_t i = 0; i < 5; ++i )
    {
        belfem::id_t tId = tG( i )->id();
        std::set< belfem::id_t > tNewAdj;
        for( belfem::uint k = 0; k < tG( i )->number_of_vertices(); ++k )
        {
            tNewAdj.insert( tG( i )->vertex( k )->id() );
        }
        EXPECT_EQ( tNewAdj, tOrigAdj[ tId ] );
    }

    belfem::graph::clear( tG );
}

// =============================================================================
// §7.4  Comparison Operators  [semantic]
// =============================================================================

TEST( GraphOperators, OpVertexDegree )
{
    belfem::graph::Vertex tA, tB;

    // give tA 1 neighbor, tB 3 neighbors (using explicit init)
    belfem::graph::Vertex tN1, tN2, tN3;
    tA.init_vertex_container( 1 );
    tA.insert_vertex( &tN1 );

    tB.init_vertex_container( 3 );
    tB.insert_vertex( &tN1 );
    tB.insert_vertex( &tN2 );
    tB.insert_vertex( &tN3 );

    // opVertexDegree: aA->number_of_vertices() < aB->number_of_vertices()
    EXPECT_TRUE( belfem::opVertexDegree( &tA, &tB ) );
    EXPECT_FALSE( belfem::opVertexDegree( &tB, &tA ) );
}

TEST( GraphOperators, OpVertexID )
{
    belfem::graph::Vertex tA, tB;
    tA.set_id( 5 );
    tB.set_id( 10 );

    EXPECT_TRUE( belfem::opVertexID( &tA, &tB ) );
    EXPECT_FALSE( belfem::opVertexID( &tB, &tA ) );
}

TEST( GraphOperators, OpVertexIndex )
{
    belfem::graph::Vertex tA, tB;
    tA.set_index( 2 );
    tB.set_index( 8 );

    EXPECT_TRUE( belfem::opVertexIndex( &tA, &tB ) );
    EXPECT_FALSE( belfem::opVertexIndex( &tB, &tA ) );
}

TEST( GraphOperators, OpVertexLevel )
{
    belfem::graph::Vertex tA, tB;
    tA.set_level( 1 );
    tB.set_level( 4 );

    EXPECT_TRUE( belfem::opVertexLevel( &tA, &tB ) );
    EXPECT_FALSE( belfem::opVertexLevel( &tB, &tA ) );
}

TEST( GraphOperators, OpVertexOwner )
{
    belfem::graph::Vertex tA, tB;
    tA.set_owner( 0 );
    tB.set_owner( 3 );

    EXPECT_TRUE( belfem::opVertexOwner( &tA, &tB ) );
    EXPECT_FALSE( belfem::opVertexOwner( &tB, &tA ) );
}

// =============================================================================
// §8  External Library Wrappers (Conditional)
// =============================================================================

#ifdef BELFEM_METIS

TEST( GraphMetis, MetisNdPath5 )
{
    belfem::Graph tG = make_path5();

    belfem::graph::metis_nd( tG );

    // verify valid permutation: all indices 0..4 unique
    std::set< belfem::index_t > tIndices;
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        tIndices.insert( tG( i )->index() );
    }
    EXPECT_EQ( tIndices.size(), 5u );

    belfem::graph::clear( tG );
}

TEST( GraphMetis, MetisPartition2 )
{
    belfem::Graph tG = make_path5();

    belfem::graph::metis_partition( tG, 2 );

    // all owners should be 0 or 1
    for( belfem::index_t i = 0; i < tG.size(); ++i )
    {
        EXPECT_GE( tG( i )->owner(), 0 );
        EXPECT_LE( tG( i )->owner(), 1 );
    }

    belfem::graph::clear( tG );
}

#endif // BELFEM_METIS
