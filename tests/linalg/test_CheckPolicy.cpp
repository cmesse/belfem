/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * The check policy of the linalg wrapper: the backends' bounds and size-conformance checks
 * follow BELFEM_ASSERTIONS_ACTIVE on Armadillo and Blaze alike, Blaze's shape exceptions stay
 * on in every build, and the dispatch headers assert the contracts the backends stop checking.
 * See src/linalg/doc/linalg_usage_guide.md, pitfall 5.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_eigen.hpp"

// The policy is a compile-time contract: one translation unit sees exactly one backend.
#if defined( BELFEM_ARMADILLO )
#  if BELFEM_ASSERTIONS_ACTIVE
#    ifndef ARMA_CHECK_CONFORMANCE
#      error "debug build, but Armadillo's conformance checks are off"
#    endif
#  else
#    ifdef ARMA_CHECK_CONFORMANCE
#      error "release build, but Armadillo's conformance checks are on"
#    endif
#  endif
#elif defined( BELFEM_BLAZE )
#  if BELFEM_ASSERTIONS_ACTIVE != BLAZE_USER_ASSERTION
#    error "BLAZE_USER_ASSERTION does not follow BELFEM_ASSERTIONS_ACTIVE"
#  endif
#endif

// -----------------------------------------------------------------------------
// every build: what stays fatal on both backends
// -----------------------------------------------------------------------------

#ifdef BELFEM_BLAZE
TEST( CheckPolicy, ShapeMismatchOnBlazeThrowsInEveryBuild )
{
    belfem::Matrix< belfem::real > tA( 2, 2, 1.0 );
    belfem::Matrix< belfem::real > tB( 3, 3, 1.0 );
    EXPECT_THROW( tA += tB, std::invalid_argument );
}
#endif

TEST( CheckPolicy, EigenOfNonFiniteMatrixIsAnErrorInEveryBuild )
{
    // Armadillo's eig_gen reports failure and empties its result; Blaze's eigen computes
    // zero eigenvalues without complaint. The wrapper checks finiteness before either
    belfem::Matrix< belfem::real > tA = {
        { std::numeric_limits< belfem::real >::quiet_NaN(), 0.0 },
        { 0.0, 1.0 } };
    belfem::Vector< belfem::real > tValues;
    EXPECT_THROW( belfem::eigen( tA, tValues ), std::runtime_error );
    EXPECT_THROW( belfem::eigen_sym( tA, tValues ), std::runtime_error );
}

TEST( CheckPolicy, EigenOfNonSquareMatrixIsAnErrorInEveryBuild )
{
    // the backend's own square check is off in release; the wrapper's is not
    belfem::Matrix< belfem::real > tA( 2, 3, 1.0 );
    belfem::Vector< belfem::real > tValues;
    EXPECT_THROW( belfem::eigen( tA, tValues ), std::runtime_error );
    EXPECT_THROW( belfem::eigen_sym( tA, tValues ), std::runtime_error );
}

TEST( CheckPolicy, CrossOfBackendOperandAndVector )
{
    // the mixed overloads forward the backend operand unchanged
    const belfem::Vector< belfem::real > tA = { 1.0, 0.0, 0.0 };   // const: the mixed overloads bind it
    const belfem::Vector< belfem::real > tB = { 0.0, 1.0, 0.0 };
    // assignment, not copy-initialization: the Armadillo Vector has no constructor from a
    // Glue expression yet ( round C of the hardening plan ), operator= takes any expression
    belfem::Vector< belfem::real > tC( 3 );
    belfem::Vector< belfem::real > tD( 3 );
    tC = belfem::cross( tA.vector_data(), tB );
    tD = belfem::cross( tA, tB.vector_data() );
    EXPECT_DOUBLE_EQ( tC( 2 ), 1.0 );
    EXPECT_DOUBLE_EQ( tD( 2 ), 1.0 );
    EXPECT_DOUBLE_EQ( tC( 0 ) + tC( 1 ) + tD( 0 ) + tD( 1 ), 0.0 );
}

// -----------------------------------------------------------------------------
// debug build only: the contracts the backends stop checking in release
// -----------------------------------------------------------------------------

#if BELFEM_ASSERTIONS_ACTIVE

TEST( CheckPolicyDebug, DotOfUnequalLengthsIsAsserted )
{
    belfem::Vector< belfem::real > tA( 3, 1.0 );
    belfem::Vector< belfem::real > tB( 4, 1.0 );
    EXPECT_THROW( belfem::dot( tA, tB ), std::runtime_error );
}

TEST( CheckPolicyDebug, DetOfNonSquareIsAsserted )
{
    belfem::Matrix< belfem::real > tA( 2, 3, 1.0 );
    EXPECT_THROW( belfem::det( tA ), std::runtime_error );
}

TEST( CheckPolicyDebug, InvOfNonSquareIsAsserted )
{
    belfem::Matrix< belfem::real > tA( 2, 3, 1.0 );
    EXPECT_THROW( belfem::inv( tA ), std::runtime_error );
}

#ifdef BELFEM_ARMADILLO
TEST( CheckPolicyDebug, ShapeMismatchOnArmadilloThrowsWhileChecksAreOn )
{
    belfem::Matrix< belfem::real > tA( 2, 2, 1.0 );
    belfem::Matrix< belfem::real > tB( 3, 3, 1.0 );
    EXPECT_THROW( tA += tB, std::logic_error );
}

TEST( CheckPolicyDebug, BackendBoundsOnArmadilloThrowWhileChecksAreOn )
{
    belfem::Matrix< belfem::real > tA( 2, 2, 1.0 );
    EXPECT_THROW( static_cast< void >( tA.matrix_data()( 5, 5 ) ), std::out_of_range );
}
#endif

#endif // BELFEM_ASSERTIONS_ACTIVE
