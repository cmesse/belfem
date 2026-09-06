/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the difference-constraint feasibility solver
 * graph::spfa_difference_constraints -- the SPFA rectifier behind the
 * non-unit thin-cut cleaning path.
 *
 * The solver answers: does a theta exist with
 *     theta( head( a ) ) - theta( tail( a ) ) <= weight( a )   for every arc a
 * and if not, produce a negative cycle as proof.
 *
 * Every test here checks the ANSWER independently rather than trusting the
 * solver's word for it: a feasible verdict is verified by testing all
 * constraints against the returned theta, and an infeasible verdict by
 * verifying the returned certificate really is a closed walk of strictly
 * negative total weight. That matters because the property under test is a
 * proof obligation, not a numeric value.
 */

#include <gtest/gtest.h>
#include <cstdint>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "fn_Graph_spfa.hpp"

using namespace belfem;

// =============================================================================
// Test helpers
// =============================================================================

namespace
{
    struct System
    {
        index_t          mNumVertices = 0 ;
        Cell< index_t >  mTails ;
        Cell< index_t >  mHeads ;
        Cell< int64_t >  mWeights ;

        void
        arc( const index_t aTail, const index_t aHead, const int64_t aWeight )
        {
            mTails.push( aTail );
            mHeads.push( aHead );
            mWeights.push( aWeight );
        }
    };

//------------------------------------------------------------------------------

    // Solve, require feasibility, and verify EVERY constraint against the
    // returned theta. A solver that returned true with a bad theta would pass
    // a bare EXPECT_TRUE; it cannot pass this.
    void
    expect_feasible( const System & aSystem )
    {
        Cell< int64_t > tTheta ;
        Cell< index_t > tCycle ;

        const bool tFeasible = graph::spfa_difference_constraints(
                aSystem.mNumVertices, aSystem.mTails, aSystem.mHeads,
                aSystem.mWeights, tTheta, tCycle );

        ASSERT_TRUE( tFeasible ) << "solver reported infeasible on a feasible system" ;

        ASSERT_EQ( tTheta.size(), aSystem.mNumVertices );

        for ( index_t a = 0 ; a < aSystem.mTails.size() ; ++a )
        {
            const int64_t tLhs = tTheta( aSystem.mHeads( a ) )
                               - tTheta( aSystem.mTails( a ) );

            EXPECT_LE( tLhs, aSystem.mWeights( a ) )
                << "arc " << a << " ( " << aSystem.mTails( a ) << " -> "
                << aSystem.mHeads( a ) << ", w = " << aSystem.mWeights( a )
                << " ) violated by the returned theta" ;
        }
    }

//------------------------------------------------------------------------------

    // Solve, require infeasibility, and verify the certificate is a genuine
    // closed walk whose weights sum strictly negative.
    //
    // Strictly negative matters: a zero-weight parent cycle is NOT a proof of
    // infeasibility ( it is a dead-end walk in the predecessor graph ), and an
    // earlier version of the solver accepted one.
    void
    expect_infeasible( const System & aSystem )
    {
        Cell< int64_t > tTheta ;
        Cell< index_t > tCycle ;

        const bool tFeasible = graph::spfa_difference_constraints(
                aSystem.mNumVertices, aSystem.mTails, aSystem.mHeads,
                aSystem.mWeights, tTheta, tCycle );

        ASSERT_FALSE( tFeasible ) << "solver reported feasible on an infeasible system" ;

        ASSERT_GT( tCycle.size(), 0u ) << "no certificate returned" ;

        int64_t tSum = 0 ;

        for ( index_t k = 0 ; k < tCycle.size() ; ++k )
        {
            const index_t tArc = tCycle( k );

            ASSERT_LT( tArc, aSystem.mTails.size() ) << "certificate names arc "
                << tArc << ", which does not exist" ;

            tSum += aSystem.mWeights( tArc );

            // the head of this arc must be the tail of the next, wrapping at
            // the end -- otherwise it is not a cycle
            const index_t tNext = tCycle( ( k + 1 ) % tCycle.size() );

            EXPECT_EQ( aSystem.mHeads( tArc ), aSystem.mTails( tNext ) )
                << "certificate is not a closed walk at position " << k ;
        }

        EXPECT_LT( tSum, 0 ) << "certificate cycle does not have negative weight" ;
    }

//------------------------------------------------------------------------------

    // A deterministic small PRNG, so a failure is reproducible. std::rand is
    // not used because its sequence is implementation-defined.
    struct Rng
    {
        uint64_t mState ;

        explicit Rng( const uint64_t aSeed ) : mState( aSeed ) {}

        uint64_t
        next()
        {
            // xorshift64*
            mState ^= mState >> 12 ;
            mState ^= mState << 25 ;
            mState ^= mState >> 27 ;
            return mState * 0x2545F4914F6CDD1DULL ;
        }

        // in [ 0, aBound )
        index_t
        below( const index_t aBound )
        {
            return static_cast< index_t >( this->next() % aBound );
        }

        int64_t
        between( const int64_t aLow, const int64_t aHigh )
        {
            return aLow + static_cast< int64_t >(
                    this->next() % static_cast< uint64_t >( aHigh - aLow + 1 ) );
        }
    };
}

// =============================================================================
// Degenerate and structural cases
// =============================================================================

TEST( GraphSpfa, NoArcsIsFeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 5 ;

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

// Vertices that appear in no arc are documented to stay at theta = 0.
TEST( GraphSpfa, IsolatedVerticesStayAtZero )
{
    System tSystem ;
    tSystem.mNumVertices = 6 ;
    tSystem.arc( 0, 1, 3 );
    tSystem.arc( 1, 2, 3 );

    Cell< int64_t > tTheta ;
    Cell< index_t > tCycle ;

    ASSERT_TRUE( graph::spfa_difference_constraints(
            tSystem.mNumVertices, tSystem.mTails, tSystem.mHeads,
            tSystem.mWeights, tTheta, tCycle ) );

    // 3, 4 and 5 are in no arc
    EXPECT_EQ( tTheta( 3 ), 0 );
    EXPECT_EQ( tTheta( 4 ), 0 );
    EXPECT_EQ( tTheta( 5 ), 0 );
}

//------------------------------------------------------------------------------

TEST( GraphSpfa, SimpleChainIsFeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 4 ;
    tSystem.arc( 0, 1, 2 );
    tSystem.arc( 1, 2, 2 );
    tSystem.arc( 2, 3, 2 );

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

// The production shape: a chain of unit constraints that must rectify. Both
// directions are present, which pins every difference exactly.
TEST( GraphSpfa, TightUnitChainIsFeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 5 ;

    for ( index_t k = 0 ; k + 1 < 5 ; ++k )
    {
        tSystem.arc( k, k + 1, 1 );
        tSystem.arc( k + 1, k, -1 );
    }

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

TEST( GraphSpfa, MultipleComponentsSolvedIndependently )
{
    System tSystem ;
    tSystem.mNumVertices = 7 ;

    // component A
    tSystem.arc( 0, 1, 1 );
    tSystem.arc( 1, 2, 1 );

    // component B, disjoint
    tSystem.arc( 4, 5, -2 );
    tSystem.arc( 5, 6, 3 );

    expect_feasible( tSystem );
}

// =============================================================================
// Self loops
// =============================================================================

// theta( v ) - theta( v ) = 0 <= w holds for any non-negative w.
TEST( GraphSpfa, NonNegativeSelfLoopIsFeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 3 ;
    tSystem.arc( 1, 1, 0 );
    tSystem.arc( 2, 2, 5 );

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

// 0 <= -1 is unsatisfiable, and the one-arc cycle is the certificate.
TEST( GraphSpfa, NegativeSelfLoopIsInfeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 3 ;
    tSystem.arc( 1, 1, -1 );

    expect_infeasible( tSystem );
}

// =============================================================================
// Cycles: the distinction the solver exists to get right
// =============================================================================

// A zero-weight cycle is FEASIBLE -- it only forces every theta on the cycle
// to be equal. Reporting it as infeasible was a real defect found in audit,
// so this test is a regression guard, not a formality.
TEST( GraphSpfa, ZeroWeightCycleIsFeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 3 ;
    tSystem.arc( 0, 1, 0 );
    tSystem.arc( 1, 2, 0 );
    tSystem.arc( 2, 0, 0 );

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

// A cycle whose weights sum to zero but individually differ is still feasible.
TEST( GraphSpfa, BalancedCycleIsFeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 3 ;
    tSystem.arc( 0, 1, 5 );
    tSystem.arc( 1, 2, -3 );
    tSystem.arc( 2, 0, -2 );

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

TEST( GraphSpfa, NegativeTriangleIsInfeasible )
{
    System tSystem ;
    tSystem.mNumVertices = 3 ;
    tSystem.arc( 0, 1, -1 );
    tSystem.arc( 1, 2, -1 );
    tSystem.arc( 2, 0, -1 );

    expect_infeasible( tSystem );
}

//------------------------------------------------------------------------------

// The greedy-hang case that motivated the whole rectifier: a 3-cycle in which
// every edge carries coefficient c = 2, which in production maps to a weight
// of 1 - c = -1. Greedy sweeps loop on this for ever; the solver must instead
// return a chained negative-weight certificate naming the throat.
TEST( GraphSpfa, GreedyHangThreeCycleReturnsCertificate )
{
    System tSystem ;
    tSystem.mNumVertices = 3 ;

    const int64_t tC = 2 ;

    tSystem.arc( 0, 1, 1 - tC );
    tSystem.arc( 1, 2, 1 - tC );
    tSystem.arc( 2, 0, 1 - tC );

    expect_infeasible( tSystem );
}

//------------------------------------------------------------------------------

// A negative cycle buried inside a larger feasible graph must still be found,
// and the certificate must name only the cycle arcs.
TEST( GraphSpfa, NegativeCycleInsideLargerGraphIsFound )
{
    System tSystem ;
    tSystem.mNumVertices = 8 ;

    // feasible tail
    tSystem.arc( 0, 1, 4 );
    tSystem.arc( 1, 2, 4 );

    // the offending cycle
    tSystem.arc( 3, 4, -2 );
    tSystem.arc( 4, 5, -2 );
    tSystem.arc( 5, 3, -2 );

    // feasible spur hanging off it
    tSystem.arc( 5, 6, 7 );
    tSystem.arc( 6, 7, 7 );

    expect_infeasible( tSystem );
}

// =============================================================================
// Randomized and stress
// =============================================================================

// Feasible BY CONSTRUCTION: draw potentials first, then give every arc a
// weight no smaller than the difference it must admit. Any "infeasible"
// verdict here is a false negative in the solver.
//
// This is also the regression guard for the spurious update-count trigger:
// the classic SPFA cnt > n test can fire before the predecessor graph holds a
// cycle, and an early version aborted on exactly this class of instance.
TEST( GraphSpfa, RandomFeasibleInstances )
{
    Rng tRng( 0x5EEDu );

    for ( uint tCase = 0 ; tCase < 200 ; ++tCase )
    {
        const index_t tNumVertices = 3 + tRng.below( 30 );

        // the hidden solution
        Cell< int64_t > tTruth ;
        for ( index_t v = 0 ; v < tNumVertices ; ++v )
        {
            tTruth.push( tRng.between( -50, 50 ) );
        }

        System tSystem ;
        tSystem.mNumVertices = tNumVertices ;

        const index_t tNumArcs = tNumVertices + tRng.below( 3 * tNumVertices );

        for ( index_t a = 0 ; a < tNumArcs ; ++a )
        {
            const index_t tTail = tRng.below( tNumVertices );
            const index_t tHead = tRng.below( tNumVertices );

            // slack >= 0 keeps the hidden solution admissible
            const int64_t tSlack = tRng.between( 0, 10 );

            tSystem.arc( tTail, tHead,
                         tTruth( tHead ) - tTruth( tTail ) + tSlack );
        }

        expect_feasible( tSystem );
    }
}

//------------------------------------------------------------------------------

// Take a feasible instance and inject one cycle that cannot be satisfied.
// The verdict must flip, and the certificate must verify.
TEST( GraphSpfa, RandomInstancesWithInjectedNegativeCycle )
{
    Rng tRng( 0xC0FFEEu );

    for ( uint tCase = 0 ; tCase < 100 ; ++tCase )
    {
        const index_t tNumVertices = 4 + tRng.below( 20 );

        Cell< int64_t > tTruth ;
        for ( index_t v = 0 ; v < tNumVertices ; ++v )
        {
            tTruth.push( tRng.between( -20, 20 ) );
        }

        System tSystem ;
        tSystem.mNumVertices = tNumVertices ;

        const index_t tNumArcs = tNumVertices + tRng.below( 2 * tNumVertices );

        for ( index_t a = 0 ; a < tNumArcs ; ++a )
        {
            const index_t tTail = tRng.below( tNumVertices );
            const index_t tHead = tRng.below( tNumVertices );

            tSystem.arc( tTail, tHead,
                         tTruth( tHead ) - tTruth( tTail ) + tRng.between( 0, 5 ) );
        }

        // a 3-cycle on distinct vertices whose weights sum to -3
        tSystem.arc( 0, 1, -1 );
        tSystem.arc( 1, 2, -1 );
        tSystem.arc( 2, 0, -1 );

        expect_infeasible( tSystem );
    }
}

//------------------------------------------------------------------------------

// Long chain: exercises the BFS warm start and the queue, and would expose an
// operation cap set below the real O( V * E ) need.
TEST( GraphSpfa, LongChainStress )
{
    const index_t tNumVertices = 2000 ;

    System tSystem ;
    tSystem.mNumVertices = tNumVertices ;

    for ( index_t k = 0 ; k + 1 < tNumVertices ; ++k )
    {
        tSystem.arc( k, k + 1, 1 );
    }

    // one long back arc, loose enough to stay feasible
    tSystem.arc( tNumVertices - 1, 0, static_cast< int64_t >( tNumVertices ) );

    expect_feasible( tSystem );
}

//------------------------------------------------------------------------------

// Same size, but the back arc closes a negative cycle.
TEST( GraphSpfa, LongChainWithNegativeReturnArcIsInfeasible )
{
    const index_t tNumVertices = 500 ;

    System tSystem ;
    tSystem.mNumVertices = tNumVertices ;

    for ( index_t k = 0 ; k + 1 < tNumVertices ; ++k )
    {
        tSystem.arc( k, k + 1, 1 );
    }

    tSystem.arc( tNumVertices - 1, 0,
                 -static_cast< int64_t >( tNumVertices ) );

    expect_infeasible( tSystem );
}
