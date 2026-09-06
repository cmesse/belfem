/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for graph::Vertex class: construction, properties, flags,
 * adjacency container lifecycle, sort, reverse, disabled element container.
 * See: tests_07_graph.md §1
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"

// =============================================================================
// §1.1  Construction & Properties  [semantic]
// =============================================================================

TEST( GraphVertex, DefaultConstruction )
{
    belfem::graph::Vertex tV;
    EXPECT_EQ( ( belfem::id_t ) tV.id(), belfem::gNoID );
    EXPECT_EQ( ( belfem::index_t ) tV.index(), belfem::gNoIndex );
    EXPECT_EQ( tV.owner(), belfem::gNoOwner );   // proc_t on both sides, no cast
    EXPECT_EQ( tV.level(), 0u );
}

TEST( GraphVertex, SetGetId )
{
    belfem::graph::Vertex tV;
    tV.set_id( 42 );
    EXPECT_EQ( tV.id(), 42u );
}

TEST( GraphVertex, SetGetIndex )
{
    belfem::graph::Vertex tV;
    tV.set_index( 7 );
    EXPECT_EQ( tV.index(), 7u );
}

TEST( GraphVertex, SetGetOwner )
{
    belfem::graph::Vertex tV;
    tV.set_owner( 3 );
    EXPECT_EQ( tV.owner(), 3 );
}

TEST( GraphVertex, SetGetLevel )
{
    belfem::graph::Vertex tV;
    tV.set_level( 5 );
    EXPECT_EQ( tV.level(), 5u );
}

// =============================================================================
// §1.2  Flag System  [semantic]
// =============================================================================

TEST( GraphVertexFlag, FlagDefault )
{
    belfem::graph::Vertex tV;
    for( uint8_t i = 0; i < 8; ++i )
    {
        EXPECT_FALSE( tV.is_flagged( i ) );
    }
}

TEST( GraphVertexFlag, FlagAndTest )
{
    belfem::graph::Vertex tV;
    tV.flag( 0 );
    EXPECT_TRUE( tV.is_flagged( 0 ) );
    EXPECT_FALSE( tV.is_flagged( 1 ) );
}

TEST( GraphVertexFlag, UnflagAndTest )
{
    belfem::graph::Vertex tV;
    tV.flag( 3 );
    EXPECT_TRUE( tV.is_flagged( 3 ) );
    tV.unflag( 3 );
    EXPECT_FALSE( tV.is_flagged( 3 ) );
}

TEST( GraphVertexFlag, MultipleFlagsIndependent )
{
    belfem::graph::Vertex tV;
    tV.flag( 0 );
    tV.flag( 5 );
    EXPECT_TRUE( tV.is_flagged( 0 ) );
    EXPECT_TRUE( tV.is_flagged( 5 ) );
    EXPECT_FALSE( tV.is_flagged( 1 ) );
    EXPECT_FALSE( tV.is_flagged( 7 ) );
}

TEST( GraphVertexFlag, AllEightFlags )
{
    belfem::graph::Vertex tV;
    for( uint8_t i = 0; i < 8; ++i )
    {
        tV.flag( i );
    }
    for( uint8_t i = 0; i < 8; ++i )
    {
        EXPECT_TRUE( tV.is_flagged( i ) );
    }
}

// --- §1.3 Flag System [debug] ---

#ifndef NDEBUG

TEST( GraphVertexDebug, FlagIndexOutOfBoundsThrows )
{
    belfem::graph::Vertex tV;
    EXPECT_THROW( tV.flag( 8 ), std::runtime_error );
}

TEST( GraphVertexDebug, UnflagIndexOutOfBoundsThrows )
{
    belfem::graph::Vertex tV;
    EXPECT_THROW( tV.unflag( 8 ), std::runtime_error );
}

TEST( GraphVertexDebug, IsFlaggedIndexOutOfBoundsThrows )
{
    belfem::graph::Vertex tV;
    EXPECT_THROW( tV.is_flagged( 8 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §1.4  Adjacency Container Lifecycle  [semantic]
// =============================================================================

TEST( GraphVertexAdj, CountAllocFillPattern )
{
    // The standard pattern: increment N times → init → insert N times
    belfem::graph::Vertex tCenter;
    belfem::graph::Vertex tA, tB, tC;
    tA.set_id( 10 );
    tB.set_id( 20 );
    tC.set_id( 30 );

    tCenter.increment_vertex_counter();
    tCenter.increment_vertex_counter();
    tCenter.increment_vertex_counter();
    tCenter.init_vertex_container();

    tCenter.insert_vertex( &tA );
    tCenter.insert_vertex( &tB );
    tCenter.insert_vertex( &tC );

    EXPECT_EQ( tCenter.number_of_vertices(), 3u );
    EXPECT_EQ( tCenter.vertex( 0 )->id(), 10u );
    EXPECT_EQ( tCenter.vertex( 1 )->id(), 20u );
    EXPECT_EQ( tCenter.vertex( 2 )->id(), 30u );
}

TEST( GraphVertexAdj, InitWithExplicitSize )
{
    belfem::graph::Vertex tV;
    belfem::graph::Vertex tN1, tN2;

    tV.init_vertex_container( 5 );
    tV.insert_vertex( &tN1 );
    tV.insert_vertex( &tN2 );

    EXPECT_EQ( tV.number_of_vertices(), 2u );
}

TEST( GraphVertexAdj, ResetClearsContainer )
{
    belfem::graph::Vertex tV;
    belfem::graph::Vertex tN1;

    tV.increment_vertex_counter();
    tV.init_vertex_container();
    tV.insert_vertex( &tN1 );
    EXPECT_EQ( tV.number_of_vertices(), 1u );

    tV.reset_vertex_container();
    EXPECT_EQ( tV.number_of_vertices(), 0u );
}

TEST( GraphVertexAdj, ReInitAfterFill )
{
    belfem::graph::Vertex tV;
    belfem::graph::Vertex tN1, tN2, tN3;

    // first fill with 3
    tV.init_vertex_container( 3 );
    tV.insert_vertex( &tN1 );
    tV.insert_vertex( &tN2 );
    tV.insert_vertex( &tN3 );
    EXPECT_EQ( tV.number_of_vertices(), 3u );

    // re-init with 2 — old buffer freed, new buffer allocated
    tV.init_vertex_container( 2 );
    tV.insert_vertex( &tN1 );
    tV.insert_vertex( &tN2 );
    EXPECT_EQ( tV.number_of_vertices(), 2u );
}

TEST( GraphVertexAdj, EmptyContainer )
{
    belfem::graph::Vertex tV;
    // counter is 0 → init should produce empty container
    tV.init_vertex_container();
    EXPECT_EQ( tV.number_of_vertices(), 0u );
}

// --- §1.5 Adjacency Access [debug] ---

#ifndef NDEBUG

TEST( GraphVertexDebug, VertexAccessOutOfBoundsThrows )
{
    belfem::graph::Vertex tV;
    belfem::graph::Vertex tN1, tN2;

    tV.init_vertex_container( 2 );
    tV.insert_vertex( &tN1 );
    tV.insert_vertex( &tN2 );

    EXPECT_THROW( tV.vertex( 3 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §1.6  Sort and Reverse  [semantic]
// =============================================================================

TEST( GraphVertexSort, SortVerticesByIndex )
{
    belfem::graph::Vertex tV;
    belfem::graph::Vertex tA, tB, tC;
    tA.set_index( 5 );
    tB.set_index( 1 );
    tC.set_index( 3 );

    tV.init_vertex_container( 3 );
    tV.insert_vertex( &tA );
    tV.insert_vertex( &tB );
    tV.insert_vertex( &tC );

    tV.sort_vertices();

    EXPECT_EQ( tV.vertex( 0 )->index(), 1u );
    EXPECT_EQ( tV.vertex( 1 )->index(), 3u );
    EXPECT_EQ( tV.vertex( 2 )->index(), 5u );
}

TEST( GraphVertexSort, ReverseVertices )
{
    belfem::graph::Vertex tV;
    belfem::graph::Vertex tA, tB, tC;
    tA.set_id( 1 );
    tB.set_id( 2 );
    tC.set_id( 3 );

    tV.init_vertex_container( 3 );
    tV.insert_vertex( &tA );
    tV.insert_vertex( &tB );
    tV.insert_vertex( &tC );

    tV.reverse_vertices();

    EXPECT_EQ( tV.vertex( 0 )->id(), 3u );
    EXPECT_EQ( tV.vertex( 1 )->id(), 2u );
    EXPECT_EQ( tV.vertex( 2 )->id(), 1u );
}

TEST( GraphVertexSort, ReverseEmpty )
{
    belfem::graph::Vertex tV;
    tV.init_vertex_container( 0 );
    // should not crash
    tV.reverse_vertices();
    EXPECT_EQ( tV.number_of_vertices(), 0u );
}

// =============================================================================
// §1.7  Element Container (Disabled)  [semantic]
// =============================================================================

// BELFEM_ERROR is always active — no #ifndef NDEBUG guard needed.

TEST( GraphVertexElement, InitElementContainerThrows )
{
    belfem::graph::Vertex tV;
    EXPECT_THROW( tV.init_element_container(), std::runtime_error );
}

TEST( GraphVertexElement, ResetElementContainerThrows )
{
    belfem::graph::Vertex tV;
    EXPECT_THROW( tV.reset_element_container(), std::runtime_error );
}
