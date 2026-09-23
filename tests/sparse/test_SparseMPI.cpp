/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Tier 2 MPI tests for the parallel nested-dissection wiring:
 * block_distribution + build_pargraph_adjacency, parmetis_nd ( including its
 * empty-rank fallback ), ptscotch_nd ( including an empty rank ), and the
 * PETSc solve that reaches them through DistMatrix.
 *
 * Requires: mpirun -np N ( 2 or 4 ). Registered through the TESTRANKS path in
 * config/scripts/Add_Test.cmake, never by hand. main(), gComm/gLog and the MPI
 * lifecycle live in test_sparsempi_main.cpp, which carries the launcher
 * sentinel and the verdict fold.
 *
 * Every test runs its collective calls on EVERY rank. A test that hangs is a
 * failure this binary is meant to surface; the CTest TIMEOUT bounds it.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Graph_Vertex.hpp"
#include "fn_Graph_clear.hpp"
#include "graphtools.hpp"
#include "fn_Graph_ParMETIS.hpp"
#include "fn_Graph_PTSCOTCH.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Solver.hpp"
#include "cl_SolverParameters.hpp"
#include "en_SolverEnums.hpp"

extern belfem::Communicator gComm;
extern belfem::Logger       gLog;

namespace
{
    // -------------------------------------------------------------------------
    // fixtures
    // -------------------------------------------------------------------------

    belfem::Graph
    build_graph( const belfem::Cell< belfem::Cell< belfem::index_t > > & aAdjacency )
    {
        const belfem::index_t tN = aAdjacency.size();
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

    // m x m grid, 4-neighbour stencil, no self loops: 2 m ( m - 1 ) undirected edges
    belfem::Graph
    make_grid( const belfem::index_t aM )
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( aM * aM, belfem::Cell< belfem::index_t >() );
        for( belfem::index_t i = 0; i < aM; ++i )
        {
            for( belfem::index_t j = 0; j < aM; ++j )
            {
                const belfem::index_t k = i * aM + j;
                if( i > 0 )      tAdj( k ).push( k - aM );
                if( i < aM - 1 ) tAdj( k ).push( k + aM );
                if( j > 0 )      tAdj( k ).push( k - 1 );
                if( j < aM - 1 ) tAdj( k ).push( k + 1 );
            }
        }
        return build_graph( tAdj );
    }

    // path on aN vertices
    belfem::Graph
    make_path( const belfem::index_t aN )
    {
        belfem::Cell< belfem::Cell< belfem::index_t > > tAdj( aN, belfem::Cell< belfem::index_t >() );
        for( belfem::index_t k = 0; k < aN; ++k )
        {
            if( k > 0 )      tAdj( k ).push( k - 1 );
            if( k < aN - 1 ) tAdj( k ).push( k + 1 );
        }
        return build_graph( tAdj );
    }

    // the split block_distribution promises: N/P on the first P-(N mod P) ranks
    belfem::index_t
    expected_block_size( const belfem::index_t aN, const belfem::proc_t aP, const belfem::proc_t aRank )
    {
        const belfem::index_t tP     = static_cast< belfem::index_t >( aP );
        const belfem::index_t tDiv   = aN / tP;
        const belfem::index_t tSplit = tP - ( aN % tP );
        return static_cast< belfem::index_t >( aRank ) < tSplit ? tDiv : tDiv + 1;
    }

    // -------------------------------------------------------------------------
    // oracles ( root only )
    // -------------------------------------------------------------------------

    // index() is a bijection onto 0..n-1, the id set is untouched, and the
    // degree sum is unchanged - what an ordering may and may not do to a graph
    void
    expect_valid_ordering( const belfem::Graph & aGraph, const belfem::index_t aN, const belfem::index_t aDegreeSum )
    {
        ASSERT_EQ( aGraph.size(), aN );

        belfem::Cell< belfem::uint > tIndexSeen( aN, 0 );
        belfem::Cell< belfem::uint > tIdSeen( aN, 0 );
        belfem::index_t tDegreeSum = 0;

        for( belfem::graph::Vertex * tV : aGraph )
        {
            ASSERT_LT( tV->index(), aN );
            ASSERT_LT( tV->id(), aN );
            ++tIndexSeen( tV->index() );
            ++tIdSeen( tV->id() );
            tDegreeSum += tV->number_of_vertices();
        }
        for( belfem::index_t k = 0; k < aN; ++k )
        {
            EXPECT_EQ( tIndexSeen( k ), 1u ) << "index " << k;
            EXPECT_EQ( tIdSeen( k ), 1u )    << "id " << k;
        }
        EXPECT_EQ( tDegreeSum, aDegreeSum );
    }

    const belfem::index_t tM       = 8;                              // grid side
    const belfem::index_t tGridN   = tM * tM;                        // 64 vertices
    const belfem::index_t tGridDeg = 2 * ( 2 * tM * ( tM - 1 ) );    // directed entries
}

// =============================================================================
// §1  block_distribution + build_pargraph_adjacency  ( root only, no library )
// =============================================================================

TEST( ParallelND, BlockDistributionSlices )
{
    if( belfem::comm_rank() != 0 )
    {
        return;   // the builder is root-only by contract; nothing collective here
    }

    const belfem::proc_t tP = belfem::comm_size();
    belfem::Graph tG = make_grid( tM );

    belfem::graph::block_distribution( tG, tP );

    // owners in graph order are non-decreasing and sized by the split
    belfem::Cell< belfem::index_t > tCount( tP, 0 );
    belfem::proc_t tLast = 0;
    for( belfem::graph::Vertex * tV : tG )
    {
        ASSERT_GE( tV->owner(), tLast );
        ASSERT_LT( tV->owner(), tP );
        tLast = tV->owner();
        ++tCount( tV->owner() );
    }
    for( belfem::proc_t p = 0; p < tP; ++p )
    {
        EXPECT_EQ( tCount( p ), expected_block_size( tGridN, tP, p ) ) << "rank " << p;
    }

    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;

    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    ASSERT_EQ( tDistribution.length(), static_cast< size_t >( tP + 1 ) );
    EXPECT_EQ( tDistribution( 0 ), 0 );
    EXPECT_EQ( tDistribution( tP ), static_cast< int >( tGridN ) );

    int tTotalEdges = 0;
    for( belfem::proc_t p = 0; p < tP; ++p )
    {
        const int tSize = tDistribution( p + 1 ) - tDistribution( p );
        EXPECT_EQ( tSize, static_cast< int >( expected_block_size( tGridN, tP, p ) ) );
        ASSERT_EQ( tVertices( p ).length(), static_cast< size_t >( tSize + 1 ) );
        EXPECT_EQ( tVertices( p )( 0 ), 0 );

        // the CSR terminal is the slice's degree sum, and it is what the
        // buffer holds ( no empty owner here, so no placeholder )
        EXPECT_EQ( static_cast< size_t >( tVertices( p )( tSize ) ), tEdges( p ).length() );
        tTotalEdges += tVertices( p )( tSize );

        for( size_t e = 0; e < tEdges( p ).length(); ++e )
        {
            EXPECT_GE( tEdges( p )( e ), 0 );
            EXPECT_LT( tEdges( p )( e ), static_cast< int >( tGridN ) );
        }
    }
    EXPECT_EQ( tTotalEdges, static_cast< int >( tGridDeg ) );

    belfem::graph::clear( tG );
}

TEST( ParallelND, BlockDistributionEmptyRankWhenFewerVerticesThanRanks )
{
    if( belfem::comm_rank() != 0 )
    {
        return;
    }
    const belfem::proc_t tP = belfem::comm_size();
    const belfem::index_t tN = static_cast< belfem::index_t >( tP ) - 1;   // one rank must stay empty
    belfem::Graph tG = make_path( tN );

    belfem::graph::block_distribution( tG, tP );   // tDiv = 0: the two-loop split must not divide by it

    belfem::Vector< int > tDistribution;
    belfem::Cell< belfem::Vector< int > > tVertices;
    belfem::Cell< belfem::Vector< int > > tEdges;
    belfem::graph::build_pargraph_adjacency( tG, tDistribution, tVertices, tEdges );

    // rank 0 is the empty one ( remainder goes to the LAST ranks )
    EXPECT_EQ( tDistribution( 1 ) - tDistribution( 0 ), 0 );
    ASSERT_EQ( tVertices( 0 ).length(), 1u );
    EXPECT_EQ( tVertices( 0 )( 0 ), 0 );
    EXPECT_EQ( tEdges( 0 ).length(), 1u );   // the placeholder D3 is about

    belfem::graph::clear( tG );
}

// =============================================================================
// §2  parmetis_nd  ( collective )
// =============================================================================

#ifdef BELFEM_PARMETIS

TEST( ParallelND, ParmetisOrdersBlockDistributedGraph )
{
    belfem::Graph tG;
    if( belfem::comm_rank() == 0 )
    {
        tG = make_grid( tM );
        belfem::graph::block_distribution( tG, belfem::comm_size() );
    }

    belfem::graph::parmetis_nd( tG );   // every rank

    if( belfem::comm_rank() == 0 )
    {
        expect_valid_ordering( tG, tGridN, tGridDeg );
        belfem::graph::clear( tG );
    }
    else
    {
        EXPECT_EQ( tG.size(), 0u );
    }
}

// D4, production trigger: N < P leaves a rank empty under the block split.
// The verdict must reach every rank and every rank must return; root falls
// back to serial METIS and still delivers a bijection. A hang here is the
// failure this test exists for.
TEST( ParallelND, ParmetisFallsBackWhenFewerVerticesThanRanks )
{
    const belfem::index_t tN = static_cast< belfem::index_t >( belfem::comm_size() ) - 1;
    belfem::Graph tG;
    if( belfem::comm_rank() == 0 )
    {
        tG = make_path( tN );
        belfem::graph::block_distribution( tG, belfem::comm_size() );
    }

    belfem::graph::parmetis_nd( tG );

    if( belfem::comm_rank() == 0 )
    {
        expect_valid_ordering( tG, tN, 2 * ( tN - 1 ) );
        belfem::graph::clear( tG );
    }
    else
    {
        EXPECT_EQ( tG.size(), 0u );
    }
}

// D4, synthetic trigger: a caller that hands every vertex to rank 0
TEST( ParallelND, ParmetisFallsBackWhenAllOwnersOnRoot )
{
    belfem::Graph tG;
    if( belfem::comm_rank() == 0 )
    {
        tG = make_grid( tM );
        for( belfem::graph::Vertex * tV : tG )
        {
            tV->set_owner( 0 );
        }
    }

    belfem::graph::parmetis_nd( tG );

    if( belfem::comm_rank() == 0 )
    {
        expect_valid_ordering( tG, tGridN, tGridDeg );
        belfem::graph::clear( tG );
    }
}

#endif // BELFEM_PARMETIS

// =============================================================================
// §3  ptscotch_nd  ( collective; empty ranks are legal, D3 )
// =============================================================================

#ifdef BELFEM_PTSCOTCH

TEST( ParallelND, PtscotchOrdersBlockDistributedGraph )
{
    belfem::Graph tG;
    if( belfem::comm_rank() == 0 )
    {
        tG = make_grid( tM );
        belfem::graph::block_distribution( tG, belfem::comm_size() );
    }

    belfem::graph::ptscotch_nd( tG );

    if( belfem::comm_rank() == 0 )
    {
        expect_valid_ordering( tG, tGridN, tGridDeg );
        belfem::graph::clear( tG );
    }
    else
    {
        EXPECT_EQ( tG.size(), 0u );
    }
}

TEST( ParallelND, PtscotchAcceptsEmptyRank )
{
    belfem::Graph tG;
    if( belfem::comm_rank() == 0 )
    {
        tG = make_grid( tM );
        for( belfem::graph::Vertex * tV : tG )
        {
            tV->set_owner( 0 );   // every other rank gets the placeholder slice
        }
    }

    belfem::graph::ptscotch_nd( tG );

    if( belfem::comm_rank() == 0 )
    {
        expect_valid_ordering( tG, tGridN, tGridDeg );
        belfem::graph::clear( tG );
    }
}

#endif // BELFEM_PTSCOTCH

// =============================================================================
// §4  end to end: the call site. PETSc at np > 1 builds a DistMatrix, whose
//     create_matrix() runs the ordering every rank must take part in.
// =============================================================================

#ifdef BELFEM_PETSC

namespace
{
    // 2-D five-point Laplacian on an m x m grid, SPD, dense assembly is fine at 64
    belfem::Matrix< belfem::real >
    make_poisson( const belfem::index_t aM )
    {
        const belfem::index_t tN = aM * aM;
        belfem::Matrix< belfem::real > tK( tN, tN, 0.0 );
        for( belfem::index_t i = 0; i < aM; ++i )
        {
            for( belfem::index_t j = 0; j < aM; ++j )
            {
                const belfem::index_t k = i * aM + j;
                tK( k, k ) = 4.0;
                if( i > 0 )      tK( k, k - aM ) = -1.0;
                if( i < aM - 1 ) tK( k, k + aM ) = -1.0;
                if( j > 0 )      tK( k, k - 1 )  = -1.0;
                if( j < aM - 1 ) tK( k, k + 1 )  = -1.0;
            }
        }
        return tK;
    }

    // CG + Jacobi, pinned: Jacobi is invariant under a symmetric permutation,
    // so the ordering cannot change the iteration and "same solution" is a
    // statement about the permutation plumbing, not about the preconditioner.
    // aM = grid side; aM * aM rows
    void
    solve_poisson_and_verify( const belfem::ReorderingMethod aMethod, const belfem::index_t aM )
    {
        // the production contract: the matrix and the vectors exist on
        // rank 0 only; every other rank hands in an empty matrix and empty
        // vectors, and PETSC::solve broadcasts the lengths it needs. Building
        // the full matrix on every rank would let a non-root dereference of
        // aMatrix stay green. SpMatrix is neither copyable nor movable, so the
        // root matrix is built through a pointer
        const bool tRoot = belfem::comm_rank() == 0;
        const belfem::index_t tN = aM * aM;

        belfem::Matrix< belfem::real > tK;
        belfem::Vector< belfem::real > tXknown;
        belfem::Vector< belfem::real > tB;
        belfem::Vector< belfem::real > tX;
        belfem::SpMatrix * tA = nullptr;

        if( tRoot )
        {
            tK = make_poisson( aM );
            tA = new belfem::SpMatrix( tK, belfem::SpMatrixType::CSR );

            tXknown.set_size( tN );
            for( belfem::index_t k = 0; k < tN; ++k )
            {
                tXknown( k ) = 1.0 + 0.01 * static_cast< belfem::real >( k );
            }
            tB = tK * tXknown;
            tX.set_size( tN, 0.0 );
        }
        else
        {
            tA = new belfem::SpMatrix();
        }

        belfem::SolverParameters tParams( belfem::SolverType::PETSc );
        tParams.set_reordering_method( aMethod );
        tParams.set_krylov_method( belfem::KrylovMethod::CG );
        tParams.set_preconditioner( belfem::Preconditioner::JACOBI );
        tParams.set_relative_tolerance( 1e-12 );

        belfem::Solver tSolver( tParams );
        tSolver.solve( *tA, tX, tB );

        if( tRoot )
        {
            for( belfem::index_t k = 0; k < tN; ++k )
            {
                EXPECT_NEAR( tX( k ), tXknown( k ), 1e-8 ) << "row " << k;
            }
        }

        delete tA;
    }
}

// the control: today's path, so a failure below is attributable to the wiring
#ifdef BELFEM_METIS
TEST( ParallelND, PetscSolvesWithMetisOrdering )
{
    solve_poisson_and_verify( belfem::ReorderingMethod::METIS, tM );
}
#endif

#ifdef BELFEM_PARMETIS
TEST( ParallelND, PetscSolvesWithParmetisOrdering )
{
    solve_poisson_and_verify( belfem::ReorderingMethod::PARMETIS, tM );
}

// the D4 fallback through the whole DistMatrix tail: a 1 x 1 grid is one row,
// fewer than any Tier 2 rank count, so parmetis_nd falls back on every rank
// and create_matrix must still deliver a solvable permuted matrix
TEST( ParallelND, PetscSolvesWithParmetisFallbackOnFewerRowsThanRanks )
{
    solve_poisson_and_verify( belfem::ReorderingMethod::PARMETIS, 1 );
}
#endif

#ifdef BELFEM_PTSCOTCH
TEST( ParallelND, PetscSolvesWithPtscotchOrdering )
{
    solve_poisson_and_verify( belfem::ReorderingMethod::PTSCOTCH, tM );
}
#endif

#endif // BELFEM_PETSC

// =============================================================================
// §5  SolverParameters reach every rank
// =============================================================================

// rank 0 states every member that has a setter away from its default ( the
// METIS NodeNDP flag has none and stays at its default ), the other ranks keep
// the defaults, and synchronize() must leave every rank with rank 0's record,
// have-flags included. The payload is positional, so a member dropped or
// shifted in pack or unpack desyncs the members behind it
TEST( ParallelND, SolverParametersSynchronizeAcrossRanks )
{
    belfem::SolverParameters tParams( belfem::SolverType::PETSc );

    if( belfem::comm_rank() == 0 )
    {
        tParams.set_distributed_matrix_type( belfem::DistributedMatrixType::CSC );
        tParams.set_preconditioner( belfem::Preconditioner::GAMG );
        tParams.set_krylov_method( belfem::KrylovMethod::GMRES );
        tParams.set_reordering_method( belfem::ReorderingMethod::NATURAL );
        tParams.set_compression_method( belfem::CompressionMethod::BLR );
        tParams.set_use_initial_guess( true );
        tParams.set_matrix_matching( false );
        tParams.set_relative_tolerance( 1e-6 );
        tParams.set_compression_cutoff( 1e-13 );
        tParams.set_absolute_tolerance( 1e-12 );
        tParams.set_memory_budget( 4000 );
        tParams.set_max_iterations( 7 );
    }

    tParams.synchronize();

    EXPECT_EQ( tParams.distributed_matrix_type(), belfem::DistributedMatrixType::CSC );
    EXPECT_EQ( tParams.preconditioner(), belfem::Preconditioner::GAMG );
    EXPECT_EQ( tParams.krylov_method(), belfem::KrylovMethod::GMRES );
    EXPECT_EQ( tParams.reordering_method(), belfem::ReorderingMethod::NATURAL );
    EXPECT_EQ( tParams.compression_method(), belfem::CompressionMethod::BLR );
    EXPECT_TRUE( tParams.use_initial_guess() );
    EXPECT_FALSE( tParams.use_matrix_matching() );
    EXPECT_TRUE( tParams.have_relative_tolerance() );
    EXPECT_EQ( tParams.relative_tolerance(), 1e-6 );
    EXPECT_TRUE( tParams.have_compression_cutoff() );
    EXPECT_EQ( tParams.compression_cutoff(), 1e-13 );
    EXPECT_TRUE( tParams.have_absolute_tolerance() );
    EXPECT_EQ( tParams.absolute_tolerance(), 1e-12 );
    EXPECT_TRUE( tParams.have_memory_budget() );
    EXPECT_EQ( tParams.memory_budget(), 4000u );
    EXPECT_EQ( tParams.max_iterations(), 7u );
}

// =============================================================================
// §6  STRUMPACK honors the caller's initial iterate in parallel
// =============================================================================

#ifdef BELFEM_STRUMPACK

namespace
{
    // unpreconditioned GMRES on the n x n 1-D Laplacian ( condition number
    // ~ 4 n^2 / pi^2 ) cannot reach the exit test within STRUMPACK's 50
    // iterations from a zero start, so the answer depends on whether the
    // guess reached the solver. Returns the max-norm error on rank 0, 0
    // elsewhere
    belfem::real
    strumpack_error_from_guess( const bool aUseGuess )
    {
        const belfem::index_t tN = 200;
        const bool tRoot = belfem::comm_rank() == 0;

        belfem::Vector< belfem::real > tXknown;
        belfem::Vector< belfem::real > tB;
        belfem::Vector< belfem::real > tX;
        belfem::SpMatrix * tA = nullptr;

        if( tRoot )
        {
            belfem::Matrix< belfem::real > tK( tN, tN, 0.0 );
            tXknown.set_size( tN );
            for( belfem::index_t k = 0; k < tN; ++k )
            {
                tK( k, k ) = 2.0;
                if( k > 0 )      tK( k, k - 1 ) = -1.0;
                if( k < tN - 1 ) tK( k, k + 1 ) = -1.0;
                tXknown( k ) = std::sin( 0.1 * static_cast< belfem::real >( k ) ) + 2.0;
            }
            tA = new belfem::SpMatrix( tK, belfem::SpMatrixType::CSR );
            tB = tK * tXknown;

            // the guess IS the solution
            tX = tXknown;
        }
        else
        {
            tA = new belfem::SpMatrix();
        }

        belfem::SolverParameters tParams( belfem::SolverType::STRUMPACK );
        tParams.set_krylov_method( belfem::KrylovMethod::GMRES );
        tParams.set_preconditioner( belfem::Preconditioner::NONE );
        tParams.set_use_initial_guess( aUseGuess );

        belfem::Solver tSolver( tParams );
        tSolver.solve( *tA, tX, tB );

        belfem::real tError = 0.0;
        if( tRoot )
        {
            // written so that a NaN entry makes the error NaN: std::max
            // would drop it, and an undistributed guess is uninitialized
            // memory that can hold one
            for( belfem::index_t k = 0; k < tN; ++k )
            {
                const belfem::real tDiff = std::abs( tX( k ) - tXknown( k ) );
                if( ! ( tDiff <= tError ) )
                {
                    tError = tDiff;
                }
            }
        }

        delete tA;
        return tError;
    }
}

TEST( ParallelND, StrumpackHonorsTheDistributedInitialIterate )
{
    const belfem::real tWithGuess    = strumpack_error_from_guess( true );
    const belfem::real tWithoutGuess = strumpack_error_from_guess( false );

    if( belfem::comm_rank() == 0 )
    {
        std::printf( "STRUMPACK max error: with guess %.3e, without %.3e\n", tWithGuess, tWithoutGuess );

        EXPECT_LT( tWithGuess, 1e-10 );

        // the control: without the guess the same solve lands far off, so
        // the first expectation is a statement about the guess
        EXPECT_GT( tWithoutGuess, 1e-3 );
    }
}

#endif // BELFEM_STRUMPACK
