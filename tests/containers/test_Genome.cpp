/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Genome<B,N> container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §9 (test matrix)
 *
 * Genome is not a container in the strict sense, but exercises Bitset and
 * Vector in a non-trivial way.  Tests are round-trip functional tests only
 * — we do NOT test private encode/decode internals.
 *
 * No BELFEM_ASSERT/BELFEM_ERROR in this class → no debug-only tests.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>    // std::srand

#include "typedefs.hpp"
#include "cl_Genome.hpp"
#include "cl_Vector.hpp"
#include "cl_Bitset.hpp"

// BELFEM_REAL_MAX macro expands to std::numeric_limits<real>::max()
// which needs unqualified 'real' in scope.
using belfem::real;

// =============================================================================
// Test configuration
// =============================================================================
//
// B = 8 bits per chromosome → 2^8 - 1 = 255 quantization levels
// N = 3 chromosomes (parameters)
//
// Maximum round-trip quantization error per parameter (linear scale):
//   (max - min) / (2 * (2^B - 1))  ≈  (max - min) / 510
//
// We use a conservative tolerance of (max - min) / (2^B - 1) which is
// roughly 2× the theoretical max error.

static constexpr size_t B = 8;
static constexpr size_t N = 3;

// Helper: compute per-parameter quantization tolerance for linear scale
static belfem::real quantization_tol( belfem::real aMin, belfem::real aMax )
{
    return ( aMax - aMin ) / static_cast< belfem::real >( ( 1u << B ) - 1 )
           + belfem::BELFEM_EPSILON;
}

// =============================================================================
// §9.1  Tests  [semantic]
// =============================================================================

TEST( Genome, SetValuesGetValuesRoundTrip )
{
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;   // all linear (default: all zero)

    tMin( 0 ) = 0.0;   tMax( 0 ) = 10.0;
    tMin( 1 ) = -5.0;  tMax( 1 ) = 5.0;
    tMin( 2 ) = 100.0; tMax( 2 ) = 200.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // set known values
    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = 3.5;
    tInput( 1 ) = -1.0;
    tInput( 2 ) = 150.0;

    tGenome.set_values( tInput );

    // decode
    belfem::Vector< belfem::real > tOutput( N );
    tGenome.get_values( tOutput );

    // verify round-trip within quantization tolerance
    for( size_t i = 0; i < N; ++i )
    {
        EXPECT_NEAR( tOutput( i ), tInput( i ),
                     quantization_tol( tMin( i ), tMax( i ) ) );
    }
}

TEST( Genome, LinearScaleRoundTrip )
{
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;   // all linear

    tMin( 0 ) = 0.0;    tMax( 0 ) = 1.0;
    tMin( 1 ) = -100.0; tMax( 1 ) = 100.0;
    tMin( 2 ) = 50.0;   tMax( 2 ) = 51.0;   // narrow range

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // test at midpoints
    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = 0.5;
    tInput( 1 ) = 0.0;
    tInput( 2 ) = 50.5;

    tGenome.set_values( tInput );

    belfem::Vector< belfem::real > tOutput( N );
    tGenome.get_values( tOutput );

    for( size_t i = 0; i < N; ++i )
    {
        EXPECT_NEAR( tOutput( i ), tInput( i ),
                     quantization_tol( tMin( i ), tMax( i ) ) );
    }
}

TEST( Genome, LogScaleRoundTrip )
{
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    // all log-scale (positive values required)
    tTypes.set( 0 );
    tTypes.set( 1 );
    tTypes.set( 2 );

    tMin( 0 ) = 1.0;    tMax( 0 ) = 1000.0;
    tMin( 1 ) = 0.001;  tMax( 1 ) = 1.0;
    tMin( 2 ) = 10.0;   tMax( 2 ) = 100.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = 100.0;
    tInput( 1 ) = 0.01;
    tInput( 2 ) = 50.0;

    tGenome.set_values( tInput );

    belfem::Vector< belfem::real > tOutput( N );
    tGenome.get_values( tOutput );

    // log-scale tolerance: quantization step in log space scaled to real space
    for( size_t i = 0; i < N; ++i )
    {
        belfem::real tLogMin = std::log( tMin( i ) );
        belfem::real tLogMax = std::log( tMax( i ) );
        belfem::real tLogTol = ( tLogMax - tLogMin ) / static_cast< belfem::real >( ( 1u << B ) - 1 );

        // in real space, tolerance scales with the value
        belfem::real tExpectedTol = tOutput( i ) * ( std::exp( tLogTol ) - 1.0 ) + belfem::BELFEM_EPSILON;
        EXPECT_NEAR( tOutput( i ), tInput( i ), tExpectedTol );
    }
}

TEST( Genome, ClampingBehavior )
{
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;   // all linear

    tMin( 0 ) = 0.0;  tMax( 0 ) = 10.0;
    tMin( 1 ) = 0.0;  tMax( 1 ) = 10.0;
    tMin( 2 ) = 0.0;  tMax( 2 ) = 10.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // set values outside [min, max]
    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = -5.0;    // below min → clamped to 0
    tInput( 1 ) = 15.0;    // above max → clamped to 10
    tInput( 2 ) = 5.0;     // within range → unchanged

    tGenome.set_values( tInput );

    belfem::Vector< belfem::real > tOutput( N );
    tGenome.get_values( tOutput );

    belfem::real tTol = quantization_tol( 0.0, 10.0 );

    EXPECT_NEAR( tOutput( 0 ), 0.0,  tTol );   // clamped to min
    EXPECT_NEAR( tOutput( 1 ), 10.0, tTol );   // clamped to max
    EXPECT_NEAR( tOutput( 2 ), 5.0,  tTol );   // unchanged
}

TEST( Genome, EqualMinMaxIsPassThrough )
{
    // Only param 1 has min==max; others are normal ranges.
    // More realistic than setting all three equal (idea from ChatGPT).
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;   tMax( 0 ) = 10.0;
    tMin( 1 ) = 42.0;  tMax( 1 ) = 42.0;   // min == max → passthrough
    tMin( 2 ) = -5.0;  tMax( 2 ) = 5.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // set_values with arbitrary input for param 1
    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = 5.0;
    tInput( 1 ) = 99.0;   // should be ignored, returns 42.0
    tInput( 2 ) = -2.5;

    tGenome.set_values( tInput );

    belfem::Vector< belfem::real > tOutput( N );
    tGenome.get_values( tOutput );

    // param 1: passthrough → returns min regardless of input
    EXPECT_DOUBLE_EQ( tOutput( 1 ), 42.0 );

    // other params: normal round-trip
    EXPECT_NEAR( tOutput( 0 ), tInput( 0 ),
                 quantization_tol( tMin( 0 ), tMax( 0 ) ) );
    EXPECT_NEAR( tOutput( 2 ), tInput( 2 ),
                 quantization_tol( tMin( 2 ), tMax( 2 ) ) );
}

TEST( Genome, RandomizeProducesValidValues )
{
    // Include a min==max param to verify passthrough works with randomize
    // (idea from ChatGPT)
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;   tMax( 0 ) = 10.0;
    tMin( 1 ) = -5.0;  tMax( 1 ) = 5.0;
    tMin( 2 ) = 7.0;   tMax( 2 ) = 7.0;   // min == max

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // seed for reproducibility (non-MPI: rand() uses std::rand())
    std::srand( 42 );

    // run randomize several times, verify bounds each time
    for( int tTrial = 0; tTrial < 10; ++tTrial )
    {
        tGenome.randomize();

        belfem::Vector< belfem::real > tOutput( N );
        tGenome.get_values( tOutput );

        EXPECT_GE( tOutput( 0 ), tMin( 0 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 0 ), tMax( 0 ) + belfem::BELFEM_EPSILON );

        EXPECT_GE( tOutput( 1 ), tMin( 1 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 1 ), tMax( 1 ) + belfem::BELFEM_EPSILON );

        // min==max param: should always return 7.0
        EXPECT_DOUBLE_EQ( tOutput( 2 ), 7.0 );
    }
}

TEST( Genome, RandomizeLogScaleBranch )
{
    // Codex finding: randomize() never tested the log-scale branch.
    // tTypes.test(i)==true means log-scale; min/max must be > 0.
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 1.0;     tMax( 0 ) = 1000.0;  tTypes.set( 0 );   // log-scale
    tMin( 1 ) = 0.01;    tMax( 1 ) = 100.0;    tTypes.set( 1 );   // log-scale
    tMin( 2 ) = 0.0;     tMax( 2 ) = 10.0;     // linear (bit 2 stays false)

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    std::srand( 123 );

    for( int tTrial = 0; tTrial < 10; ++tTrial )
    {
        tGenome.randomize();

        belfem::Vector< belfem::real > tOutput( N );
        tGenome.get_values( tOutput );

        EXPECT_GE( tOutput( 0 ), tMin( 0 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 0 ), tMax( 0 ) + belfem::BELFEM_EPSILON );

        EXPECT_GE( tOutput( 1 ), tMin( 1 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 1 ), tMax( 1 ) + belfem::BELFEM_EPSILON );

        EXPECT_GE( tOutput( 2 ), tMin( 2 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 2 ), tMax( 2 ) + belfem::BELFEM_EPSILON );
    }
}

TEST( Genome, InheritProducesValidValues )
{
    // Include a min==max param to verify passthrough works with inherit
    // (idea from ChatGPT)
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;   tMax( 0 ) = 10.0;
    tMin( 1 ) = -5.0;  tMax( 1 ) = 5.0;
    tMin( 2 ) = 9.0;   tMax( 2 ) = 9.0;   // min == max

    belfem::Genome< B, N > tMom( tMin, tMax, tTypes );
    belfem::Genome< B, N > tDad( tMin, tMax, tTypes );
    belfem::Genome< B, N > tChild( tMin, tMax, tTypes );

    std::srand( 123 );

    // give parents known values
    belfem::Vector< belfem::real > tMomVals( N );
    tMomVals( 0 ) = 3.0;
    tMomVals( 1 ) = -2.0;
    tMomVals( 2 ) = 9.0;
    tMom.set_values( tMomVals );

    belfem::Vector< belfem::real > tDadVals( N );
    tDadVals( 0 ) = 7.0;
    tDadVals( 1 ) = 2.0;
    tDadVals( 2 ) = 9.0;
    tDad.set_values( tDadVals );

    // inherit and verify bounds (stochastic due to crossover + mutation)
    for( int tTrial = 0; tTrial < 10; ++tTrial )
    {
        tChild.inherit( &tMom, &tDad );

        belfem::Vector< belfem::real > tOutput( N );
        tChild.get_values( tOutput );

        // child values must be within [min, max] even after mutation
        EXPECT_GE( tOutput( 0 ), tMin( 0 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 0 ), tMax( 0 ) + belfem::BELFEM_EPSILON );

        EXPECT_GE( tOutput( 1 ), tMin( 1 ) - belfem::BELFEM_EPSILON );
        EXPECT_LE( tOutput( 1 ), tMax( 1 ) + belfem::BELFEM_EPSILON );

        // min==max param: should always return 9.0
        EXPECT_DOUBLE_EQ( tOutput( 2 ), 9.0 );
    }
}

TEST( Genome, KillResetsState )
{
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;  tMax( 0 ) = 10.0;
    tMin( 1 ) = 0.0;  tMax( 1 ) = 10.0;
    tMin( 2 ) = 0.0;  tMax( 2 ) = 10.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // make genome alive
    tGenome.set_fitness( 1.0 );
    EXPECT_TRUE( tGenome.is_alive() );

    tGenome.kill();

    EXPECT_FALSE( tGenome.is_alive() );
    EXPECT_DOUBLE_EQ( tGenome.fitness(), BELFEM_REAL_MAX );
}

TEST( Genome, FitnessSetAndGet )
{
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;  tMax( 0 ) = 1.0;
    tMin( 1 ) = 0.0;  tMax( 1 ) = 1.0;
    tMin( 2 ) = 0.0;  tMax( 2 ) = 1.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // initial fitness is BELFEM_REAL_MAX (not alive)
    EXPECT_DOUBLE_EQ( tGenome.fitness(), BELFEM_REAL_MAX );
    EXPECT_FALSE( tGenome.is_alive() );

    // set fitness → alive
    tGenome.set_fitness( 42.5 );
    EXPECT_DOUBLE_EQ( tGenome.fitness(), 42.5 );
    EXPECT_TRUE( tGenome.is_alive() );

    // overwrite
    tGenome.set_fitness( 0.001 );
    EXPECT_DOUBLE_EQ( tGenome.fitness(), 0.001 );
    EXPECT_TRUE( tGenome.is_alive() );
}

TEST( Genome, EncodingDoesNotMakeGenomeAlive )
{
    // set_values() should NOT affect fitness or alive state (idea from ChatGPT)
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;  tMax( 0 ) = 10.0;
    tMin( 1 ) = 0.0;  tMax( 1 ) = 10.0;
    tMin( 2 ) = 0.0;  tMax( 2 ) = 10.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = 5.0;
    tInput( 1 ) = 3.0;
    tInput( 2 ) = 7.0;

    tGenome.set_values( tInput );

    EXPECT_FALSE( tGenome.is_alive() );
    EXPECT_DOUBLE_EQ( tGenome.fitness(), BELFEM_REAL_MAX );
}

// =============================================================================
// Supplemental: boundary values
// =============================================================================

TEST( Genome, RoundTripAtExactBounds )
{
    // encode values exactly at min and max — these should quantize
    // to the integer endpoints 0 and (2^B - 1) with zero error
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;   tMax( 0 ) = 10.0;
    tMin( 1 ) = -5.0;  tMax( 1 ) = 5.0;
    tMin( 2 ) = 100.0; tMax( 2 ) = 200.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    // set values to exact min
    belfem::Vector< belfem::real > tInputMin( N );
    tInputMin( 0 ) = tMin( 0 );
    tInputMin( 1 ) = tMin( 1 );
    tInputMin( 2 ) = tMin( 2 );

    tGenome.set_values( tInputMin );

    belfem::Vector< belfem::real > tOutputMin( N );
    tGenome.get_values( tOutputMin );

    for( size_t i = 0; i < N; ++i )
    {
        EXPECT_DOUBLE_EQ( tOutputMin( i ), tMin( i ) );
    }

    // set values to exact max
    belfem::Vector< belfem::real > tInputMax( N );
    tInputMax( 0 ) = tMax( 0 );
    tInputMax( 1 ) = tMax( 1 );
    tInputMax( 2 ) = tMax( 2 );

    tGenome.set_values( tInputMax );

    belfem::Vector< belfem::real > tOutputMax( N );
    tGenome.get_values( tOutputMax );

    for( size_t i = 0; i < N; ++i )
    {
        EXPECT_DOUBLE_EQ( tOutputMax( i ), tMax( i ) );
    }
}

TEST( Genome, MixedLinearAndLogScale )
{
    // param 0: linear, param 1: log, param 2: linear
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tTypes.set( 1 );   // only param 1 is log-scale

    tMin( 0 ) = 0.0;    tMax( 0 ) = 10.0;
    tMin( 1 ) = 0.01;   tMax( 1 ) = 100.0;    // log-scale, positive
    tMin( 2 ) = -50.0;  tMax( 2 ) = 50.0;

    belfem::Genome< B, N > tGenome( tMin, tMax, tTypes );

    belfem::Vector< belfem::real > tInput( N );
    tInput( 0 ) = 5.0;
    tInput( 1 ) = 1.0;     // log-scale
    tInput( 2 ) = -25.0;

    tGenome.set_values( tInput );

    belfem::Vector< belfem::real > tOutput( N );
    tGenome.get_values( tOutput );

    // linear params: standard tolerance
    EXPECT_NEAR( tOutput( 0 ), tInput( 0 ),
                 quantization_tol( tMin( 0 ), tMax( 0 ) ) );
    EXPECT_NEAR( tOutput( 2 ), tInput( 2 ),
                 quantization_tol( tMin( 2 ), tMax( 2 ) ) );

    // log-scale param: relative tolerance
    belfem::real tLogMin = std::log( tMin( 1 ) );
    belfem::real tLogMax = std::log( tMax( 1 ) );
    belfem::real tLogTol = ( tLogMax - tLogMin ) / static_cast< belfem::real >( ( 1u << B ) - 1 );
    belfem::real tTol = tOutput( 1 ) * ( std::exp( tLogTol ) - 1.0 ) + belfem::BELFEM_EPSILON;
    EXPECT_NEAR( tOutput( 1 ), tInput( 1 ), tTol );
}

// =============================================================================
// Supplemental: opGenomeSort comparator
// =============================================================================

TEST( Genome, GenomeSortComparatorUsesFitness )
{
    // opGenomeSort is a public functor for sorting genomes by fitness
    // (idea from ChatGPT)
    belfem::Vector< belfem::real > tMin( N );
    belfem::Vector< belfem::real > tMax( N );
    belfem::Bitset< N >    tTypes;

    tMin( 0 ) = 0.0;  tMax( 0 ) = 1.0;
    tMin( 1 ) = 0.0;  tMax( 1 ) = 1.0;
    tMin( 2 ) = 0.0;  tMax( 2 ) = 1.0;

    belfem::Genome< B, N > tA( tMin, tMax, tTypes );
    belfem::Genome< B, N > tB( tMin, tMax, tTypes );

    tA.set_fitness( 2.0 );
    tB.set_fitness( 5.0 );

    belfem::opGenomeSort< B, N > tComp;

    // A has lower fitness → should sort before B
    EXPECT_TRUE( tComp( &tA, &tB ) );
    EXPECT_FALSE( tComp( &tB, &tA ) );

    // equal fitness → neither is "less than" the other
    tB.set_fitness( 2.0 );
    EXPECT_FALSE( tComp( &tA, &tB ) );
    EXPECT_FALSE( tComp( &tB, &tA ) );
}
