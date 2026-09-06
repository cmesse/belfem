/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Quaternion<T> class and related free functions
 * See: tests_00_strategy.md, tests_03_quaternion.md
 *
 * Quaternion is a plain value type (inline storage, trivially copyable);
 * the earlier external-buffer mode was removed — see cl_Quaternion.hpp.
 * Notation: q = (a, b, c, d) = (w, x, y, z) = (scalar, vector)
 */

#include <gtest/gtest.h>
#include <cstring>
#include <type_traits>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Quaternion.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_quaternion_to_rotation_matrix.hpp"
#include "fn_quaternion_from_rotation_matrix.hpp"
#include "fn_quaternion_rotate_vector.hpp"

namespace
{
    const belfem::real tEps = 100.0 * belfem::BELFEM_EPSILON;  // unit-norm, identity comparisons
    const belfem::real tTol = 1e-12;                    // rotation round-trips

    // Returns true if q1 and q2 represent the same rotation (q ≈ ±q)
    bool same_rotation(
        const belfem::Quaternion< belfem::real > & aQ1,
        const belfem::Quaternion< belfem::real > & aQ2 )
    {
        belfem::real tDotAbs = std::abs( belfem::dot( aQ1, aQ2 ) );
        return std::abs( tDotAbs - 1.0 ) < tEps;
    }
}

// =============================================================================
// §1.1  Owned-Mode Construction  [semantic]
// =============================================================================

TEST( Quaternion, DefaultConstructorZeroes )
{
    belfem::Quaternion< belfem::real > tQ;
    EXPECT_NEAR( tQ.a(), 0.0, tTol );
    EXPECT_NEAR( tQ.b(), 0.0, tTol );
    EXPECT_NEAR( tQ.c(), 0.0, tTol );
    EXPECT_NEAR( tQ.d(), 0.0, tTol );
    EXPECT_NE( tQ.data(), nullptr );
}

TEST( Quaternion, FourArgConstructor )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    EXPECT_NEAR( tQ.a(), 1.0, tTol );
    EXPECT_NEAR( tQ.b(), 2.0, tTol );
    EXPECT_NEAR( tQ.c(), 3.0, tTol );
    EXPECT_NEAR( tQ.d(), 4.0, tTol );
}

TEST( Quaternion, FromVector3PureQuaternion )
{
    belfem::Vector< belfem::real > tV = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tV );
    EXPECT_NEAR( tQ.a(), 0.0, tTol );
    EXPECT_NEAR( tQ.b(), 1.0, tTol );
    EXPECT_NEAR( tQ.c(), 2.0, tTol );
    EXPECT_NEAR( tQ.d(), 3.0, tTol );
}

TEST( Quaternion, FromAxisAngleProducesUnitQuat )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 2.0 );
    EXPECT_NEAR( tQ.norm(), 1.0, tEps );
}

TEST( Quaternion, FromAxisAngleNormalizesAxis )
{
    // unnormalized axis {2,0,0} should produce same result as {1,0,0}
    belfem::Vector< belfem::real > tAxis2 = { 2.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tAxis1 = { 1.0, 0.0, 0.0 };
    belfem::real tAngle = M_PI / 2.0;

    belfem::Quaternion< belfem::real > tQ1( tAxis1, tAngle );
    belfem::Quaternion< belfem::real > tQ2( tAxis2, tAngle );

    EXPECT_NEAR( tQ1.a(), tQ2.a(), tTol );
    EXPECT_NEAR( tQ1.b(), tQ2.b(), tTol );
    EXPECT_NEAR( tQ1.c(), tQ2.c(), tTol );
    EXPECT_NEAR( tQ1.d(), tQ2.d(), tTol );
}

TEST( Quaternion, IdentityStaticMethod )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    EXPECT_NEAR( tQ.a(), 1.0, tTol );
    EXPECT_NEAR( tQ.b(), 0.0, tTol );
    EXPECT_NEAR( tQ.c(), 0.0, tTol );
    EXPECT_NEAR( tQ.d(), 0.0, tTol );
}

TEST( Quaternion, InitializerListAssignment )
{
    belfem::Quaternion< belfem::real > tQ;
    tQ = { 5.0, 6.0, 7.0, 8.0 };
    EXPECT_NEAR( tQ.a(), 5.0, tTol );
    EXPECT_NEAR( tQ.d(), 8.0, tTol );
}

// =============================================================================
// §1.2  Value-Type Contract  [semantic]
// -----------------------------------------------------------------------------
// The class deliberately dropped the earlier borrowed-buffer (external
// pointer) mode: it is a plain value type, trivially copyable, safe to
// memcpy / MPI-transfer (cl_Quaternion.hpp:29-33). These tests lock that
// contract; the old §1.2 borrowed-mode tests targeted the removed API.
// =============================================================================

TEST( Quaternion, ValueTypeIsTriviallyCopyable )
{
    static_assert(
        std::is_trivially_copyable< belfem::Quaternion< belfem::real > >::value,
        "Quaternion must stay trivially copyable (MPI/memcpy contract)" );

    // memcpy round-trip preserves the value ( the MPI-transfer contract )
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB;
    std::memcpy( tB.data(), tA.data(), 4 * sizeof( belfem::real ) );
    EXPECT_NEAR( tB.a(), 1.0, tTol );
    EXPECT_NEAR( tB.b(), 2.0, tTol );
    EXPECT_NEAR( tB.c(), 3.0, tTol );
    EXPECT_NEAR( tB.d(), 4.0, tTol );
}

TEST( Quaternion, DataIsInlinePerObject )
{
    // each object carries its own inline storage — no aliasing between
    // objects, no external pointers
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( tA );
    EXPECT_NE( tA.data(), tB.data() );

    tB.a() = 99.0;
    EXPECT_NEAR( tA.a(), 1.0, tTol );
}

// =============================================================================
// §1.3  Copy & Move Semantics  [semantic]
// =============================================================================

TEST( Quaternion, CopyConstructorDeepCopies )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( tA );
    tB.a() = 99.0;
    EXPECT_NEAR( tA.a(), 1.0, tTol );
}

TEST( Quaternion, CopyAssignmentCopiesValues )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB;
    tB = tA;
    EXPECT_NE( tB.data(), tA.data() );   // separate buffers
    tB.a() = 99.0;
    EXPECT_NEAR( tA.a(), 1.0, tTol );
}

TEST( Quaternion, SelfCopyAssignment )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    tQ = tQ;
    EXPECT_NEAR( tQ.a(), 1.0, tTol );
    EXPECT_NEAR( tQ.d(), 4.0, tTol );
}

TEST( Quaternion, MoveConstructorIsTrivialCopy )
{
    // rule of zero: move == copy for the inline value type; the source
    // stays valid with its own storage ( the old owned/borrowed pointer
    // semantics are gone with the removed external-buffer mode )
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );

    belfem::Quaternion< belfem::real > tB( std::move( tA ) );

    EXPECT_NEAR( tB.a(), 1.0, tTol );
    EXPECT_NEAR( tB.d(), 4.0, tTol );
    EXPECT_NE( tB.data(), tA.data() );

    // source remains a valid, independent object
    tA.a() = 50.0 ;
    EXPECT_NEAR( tB.a(), 1.0, tTol );
}

TEST( Quaternion, MoveAssignmentCopiesValues )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB;
    tB = std::move( tA );
    EXPECT_NEAR( tB.a(), 1.0, tTol );
    EXPECT_NEAR( tB.d(), 4.0, tTol );
    // move assignment copies values, does NOT transfer ownership
    EXPECT_NO_THROW( tA.a() = 0.0 );
}

TEST( Quaternion, ConjReturnsIndependentValue )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tConj = tQ.conj();
    // result is its own value object and does not alias the source
    EXPECT_NE( tConj.data(), tQ.data() );
    EXPECT_NEAR( tConj.a(),  1.0, tTol );
    EXPECT_NEAR( tConj.b(), -2.0, tTol );
}

TEST( Quaternion, InvReturnsIndependentValue )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tInv = tQ.inv();
    EXPECT_NE( tInv.data(), tQ.data() );
}

TEST( Quaternion, BinaryOperatorReturnsIndependentValue )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 0.0, 0.0, 0.0 );
    belfem::Quaternion< belfem::real > tB( 0.0, 1.0, 0.0, 0.0 );

    auto tSum  = tA + tB;
    auto tDiff = tA - tB;
    auto tProd = tA * tB;
    EXPECT_NE( tSum.data(),  tA.data() );
    EXPECT_NE( tSum.data(),  tB.data() );
    EXPECT_NE( tDiff.data(), tA.data() );
    EXPECT_NE( tProd.data(), tA.data() );
}

// The three tests below targeted the removed external-buffer constructors
// ( FromVector / FromAxisAngle / CopyConstructor WithExternalBuffer ); their
// value-mode equivalents are §1.1 FromVector3PureQuaternion,
// FromAxisAngleProducesUnitQuat, and §1.3 CopyConstructorDeepCopies.
#if 0 // removed external-buffer API — kept for reference until sign-off
TEST( Quaternion, FromVectorWithExternalBuffer )
{
    belfem::real tBuf[4];
    belfem::Vector< belfem::real > tV = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tV, tBuf );
    EXPECT_EQ( tQ.data(), tBuf );
    EXPECT_NEAR( tBuf[0], 0.0, tTol );   // pure quaternion: w=0
    EXPECT_NEAR( tBuf[1], 1.0, tTol );
}

TEST( Quaternion, FromAxisAngleWithExternalBuffer )
{
    belfem::real tBuf[4];
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 2.0, tBuf );
    EXPECT_EQ( tQ.data(), tBuf );
    EXPECT_NEAR( tQ.norm(), 1.0, tEps );
}

TEST( Quaternion, CopyConstructorWithExternalBuffer )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::real tBuf[4];
    belfem::Quaternion< belfem::real > tB( tA, tBuf );
    EXPECT_EQ( tB.data(), tBuf );
    EXPECT_NEAR( tBuf[0], 1.0, tTol );
    EXPECT_NEAR( tBuf[3], 4.0, tTol );
    // modify original, copy should be independent (stored in tBuf)
    tA.a() = 99.0;
    EXPECT_NEAR( tB.a(), 1.0, tTol );
}
#endif // removed external-buffer API

// =============================================================================
// §2  Accessors & Iterators  [semantic]
// =============================================================================

TEST( Quaternion, AccessorsMatchData )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    EXPECT_NEAR( tQ.a(), tQ.data()[0], tTol );
    EXPECT_NEAR( tQ.b(), tQ.data()[1], tTol );
    EXPECT_NEAR( tQ.c(), tQ.data()[2], tTol );
    EXPECT_NEAR( tQ.d(), tQ.data()[3], tTol );
}

TEST( Quaternion, ConstAccessors )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    const belfem::Quaternion< belfem::real > & tRef = tQ;
    EXPECT_NEAR( tRef.a(), 1.0, tTol );
    EXPECT_NE( tRef.data(), nullptr );
}

TEST( Quaternion, IteratorRange )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    EXPECT_EQ( tQ.end() - tQ.begin(), 4 );
}

TEST( Quaternion, RangeBasedForLoop )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    belfem::real tSum = 0.0;
    size_t tCount = 0;
    for( auto tVal : tQ )
    {
        tSum += tVal;
        ++tCount;
    }
    EXPECT_EQ( tCount, 4u );
    EXPECT_NEAR( tSum, 10.0, tTol );
}

// =============================================================================
// §3.1  Compound Operators  [semantic]
// =============================================================================

TEST( Quaternion, PlusEquals )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    tA += tB;
    EXPECT_NEAR( tA.a(), 6.0,  tTol );
    EXPECT_NEAR( tA.b(), 8.0,  tTol );
    EXPECT_NEAR( tA.c(), 10.0, tTol );
    EXPECT_NEAR( tA.d(), 12.0, tTol );
}

TEST( Quaternion, MinusEquals )
{
    belfem::Quaternion< belfem::real > tA( 6.0, 8.0, 10.0, 12.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    tA -= tB;
    EXPECT_NEAR( tA.a(), 1.0, tTol );
    EXPECT_NEAR( tA.d(), 4.0, tTol );
}

TEST( Quaternion, TimesEqualsScalar )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    tQ *= 2.0;
    EXPECT_NEAR( tQ.a(), 2.0, tTol );
    EXPECT_NEAR( tQ.d(), 8.0, tTol );
}

TEST( Quaternion, DivideEqualsScalar )
{
    belfem::Quaternion< belfem::real > tQ( 2.0, 4.0, 6.0, 8.0 );
    tQ /= 2.0;
    EXPECT_NEAR( tQ.a(), 1.0, tTol );
    EXPECT_NEAR( tQ.d(), 4.0, tTol );
}

TEST( Quaternion, HamiltonProductInPlace )
{
    // i*j = k: (0,1,0,0)*(0,0,1,0) = (0,0,0,1)
    belfem::Quaternion< belfem::real > tI( 0.0, 1.0, 0.0, 0.0 );
    belfem::Quaternion< belfem::real > tJ( 0.0, 0.0, 1.0, 0.0 );
    tI *= tJ;
    EXPECT_NEAR( tI.a(), 0.0, tTol );
    EXPECT_NEAR( tI.d(), 1.0, tTol );
}

TEST( Quaternion, HamiltonProductBasisRules )
{
    // test all 6 basis products (idea from ChatGPT)
    belfem::Quaternion< belfem::real > tI( 0.0, 1.0, 0.0, 0.0 );
    belfem::Quaternion< belfem::real > tJ( 0.0, 0.0, 1.0, 0.0 );
    belfem::Quaternion< belfem::real > tK( 0.0, 0.0, 0.0, 1.0 );

    // positive cycle: i*j=k, j*k=i, k*i=j
    auto tIJ = tI * tJ;
    EXPECT_NEAR( tIJ.a(), 0.0, tTol );
    EXPECT_NEAR( tIJ.d(), 1.0, tTol );

    auto tJK = tJ * tK;
    EXPECT_NEAR( tJK.a(), 0.0, tTol );
    EXPECT_NEAR( tJK.b(), 1.0, tTol );

    auto tKI = tK * tI;
    EXPECT_NEAR( tKI.a(), 0.0, tTol );
    EXPECT_NEAR( tKI.c(), 1.0, tTol );

    // negative cycle: j*i=-k, k*j=-i, i*k=-j
    auto tJI = tJ * tI;
    EXPECT_NEAR( tJI.d(), -1.0, tTol );

    auto tKJ = tK * tJ;
    EXPECT_NEAR( tKJ.b(), -1.0, tTol );

    auto tIK = tI * tK;
    EXPECT_NEAR( tIK.c(), -1.0, tTol );
}

TEST( Quaternion, HamiltonProductNotCommutative )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    auto tAB = tA * tB;
    auto tBA = tB * tA;
    bool tDiffers = ( std::abs( tAB.b() - tBA.b() ) > tTol ) ||
                    ( std::abs( tAB.c() - tBA.c() ) > tTol ) ||
                    ( std::abs( tAB.d() - tBA.d() ) > tTol );
    EXPECT_TRUE( tDiffers );
}

TEST( Quaternion, HamiltonProductAssociative )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    belfem::Quaternion< belfem::real > tC( 0.5, -1.0, 2.0, -0.5 );
    auto tAB_C = ( tA * tB ) * tC;
    auto tA_BC = tA * ( tB * tC );
    EXPECT_NEAR( tAB_C.a(), tA_BC.a(), tTol );
    EXPECT_NEAR( tAB_C.b(), tA_BC.b(), tTol );
    EXPECT_NEAR( tAB_C.c(), tA_BC.c(), tTol );
    EXPECT_NEAR( tAB_C.d(), tA_BC.d(), tTol );
}

// =============================================================================
// §3.2  Binary Operators  [semantic]
// =============================================================================

TEST( Quaternion, BinaryAddition )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    auto tC = tA + tB;
    EXPECT_NEAR( tC.a(), 6.0, tTol );
    EXPECT_NEAR( tC.d(), 12.0, tTol );
}

TEST( Quaternion, BinarySubtraction )
{
    belfem::Quaternion< belfem::real > tA( 5.0, 6.0, 7.0, 8.0 );
    belfem::Quaternion< belfem::real > tB( 1.0, 2.0, 3.0, 4.0 );
    auto tC = tA - tB;
    EXPECT_NEAR( tC.a(), 4.0, tTol );
    EXPECT_NEAR( tC.d(), 4.0, tTol );
}

TEST( Quaternion, ScalarMultiplyRight )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tC = tQ * 3.0;
    EXPECT_NEAR( tC.a(), 3.0, tTol );
    EXPECT_NEAR( tC.d(), 12.0, tTol );
}

TEST( Quaternion, ScalarMultiplyLeft )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tC = 3.0 * tQ;
    EXPECT_NEAR( tC.a(), 3.0, tTol );
    EXPECT_NEAR( tC.d(), 12.0, tTol );
}

TEST( Quaternion, ScalarDivision )
{
    belfem::Quaternion< belfem::real > tQ( 2.0, 4.0, 6.0, 8.0 );
    auto tC = tQ / 2.0;
    EXPECT_NEAR( tC.a(), 1.0, tTol );
    EXPECT_NEAR( tC.d(), 4.0, tTol );
}

// =============================================================================
// §3.4  Equality  [semantic]
// =============================================================================

TEST( Quaternion, EqualitySameQuaternion )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    EXPECT_TRUE( tQ == tQ );
}

TEST( Quaternion, EqualityWithinEpsilon )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::real tTiny = belfem::BELFEM_EPSILON * 0.5;
    belfem::Quaternion< belfem::real > tB( 1.0 + tTiny, 2.0, 3.0, 4.0 );
    EXPECT_TRUE( tA == tB );
}

TEST( Quaternion, InequalityBeyondEpsilon )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 1.0, 2.0, 3.0, 5.0 );
    EXPECT_FALSE( tA == tB );
    EXPECT_TRUE( tA != tB );
}

TEST( Quaternion, EqualityIsNotRotationEquivalence )
{
    // q and -q represent same rotation, but operator== is component-wise
    // (idea from ChatGPT — documents intentional behavior)
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tNeg = tQ * -1.0;
    EXPECT_FALSE( tQ == tNeg );
    EXPECT_TRUE( tQ != tNeg );
}

// =============================================================================
// §4.1  Norm  [semantic]
// =============================================================================

TEST( Quaternion, NormIdentity )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    EXPECT_NEAR( tQ.norm(), 1.0, tTol );
}

TEST( Quaternion, NormZero )
{
    belfem::Quaternion< belfem::real > tQ;
    EXPECT_NEAR( tQ.norm(), 0.0, tTol );
}

TEST( Quaternion, NormScaling )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tScaled = tQ * 3.0;
    EXPECT_NEAR( tScaled.norm(), tQ.norm() * 3.0, tTol );
}

TEST( Quaternion, NormMultiplicative )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 0.5, -1.0, 2.0, -0.5 );
    auto tProduct = tA * tB;
    EXPECT_NEAR( tProduct.norm(), tA.norm() * tB.norm(), tTol );
}

// =============================================================================
// §4.2  Conjugate  [semantic]
// =============================================================================

TEST( Quaternion, ConjNegatesImaginary )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tConj = tQ.conj();
    EXPECT_NEAR( tConj.a(), 1.0,  tTol );
    EXPECT_NEAR( tConj.b(), -2.0, tTol );
    EXPECT_NEAR( tConj.c(), -3.0, tTol );
    EXPECT_NEAR( tConj.d(), -4.0, tTol );
}

TEST( Quaternion, DoubleConjIsIdentity )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tDbl = tQ.conj().conj();
    EXPECT_NEAR( tDbl.a(), tQ.a(), tTol );
    EXPECT_NEAR( tDbl.b(), tQ.b(), tTol );
    EXPECT_NEAR( tDbl.c(), tQ.c(), tTol );
    EXPECT_NEAR( tDbl.d(), tQ.d(), tTol );
}

TEST( Quaternion, QTimesConjIsScalar )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    auto tProduct = tQ * tQ.conj();
    belfem::real tNormSq = tQ.norm() * tQ.norm();
    EXPECT_NEAR( tProduct.a(), tNormSq, tTol );
    EXPECT_NEAR( tProduct.b(), 0.0, tTol );
    EXPECT_NEAR( tProduct.c(), 0.0, tTol );
    EXPECT_NEAR( tProduct.d(), 0.0, tTol );
}

// =============================================================================
// §4.3  Inverse  [semantic]
// =============================================================================

TEST( Quaternion, QTimesInvIsIdentity )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    // test both directions (idea from ChatGPT)
    auto tP1 = tQ * tQ.inv();
    auto tP2 = tQ.inv() * tQ;
    EXPECT_NEAR( tP1.a(), 1.0, tTol );
    EXPECT_NEAR( tP1.b(), 0.0, tTol );
    EXPECT_NEAR( tP1.c(), 0.0, tTol );
    EXPECT_NEAR( tP1.d(), 0.0, tTol );
    EXPECT_NEAR( tP2.a(), 1.0, tTol );
    EXPECT_NEAR( tP2.b(), 0.0, tTol );
}

TEST( Quaternion, InvOfUnitIsConj )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 3.0 );
    auto tInv = tQ.inv();
    auto tConj = tQ.conj();
    EXPECT_NEAR( tInv.a(), tConj.a(), tTol );
    EXPECT_NEAR( tInv.b(), tConj.b(), tTol );
    EXPECT_NEAR( tInv.c(), tConj.c(), tTol );
    EXPECT_NEAR( tInv.d(), tConj.d(), tTol );
}

// =============================================================================
// §4.5  Normalize  [semantic]
// =============================================================================

TEST( Quaternion, NormalizeProducesUnitNorm )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    tQ.normalize();
    EXPECT_NEAR( tQ.norm(), 1.0, tTol );
}

TEST( Quaternion, NormalizePreservesDirection )
{
    belfem::Quaternion< belfem::real > tQ( 0.0, 0.0, 0.0, 4.0 );
    tQ.normalize();
    EXPECT_NEAR( tQ.a(), 0.0, tTol );
    EXPECT_NEAR( tQ.d(), 1.0, tTol );
}

TEST( Quaternion, NormalizeReturnsSelf )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > & tRef = tQ.normalize();
    EXPECT_EQ( &tRef, &tQ );
}

// =============================================================================
// §4.6  Dot and Cross  [semantic]
// =============================================================================

TEST( Quaternion, DotProductBasic )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 0.0, 0.0, 0.0 );
    belfem::Quaternion< belfem::real > tB( 0.0, 1.0, 0.0, 0.0 );
    EXPECT_NEAR( belfem::dot( tA, tB ), 0.0, tTol );
}

TEST( Quaternion, DotProductSelf )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    EXPECT_NEAR( belfem::dot( tQ, tQ ), tQ.norm() * tQ.norm(), tTol );
}

TEST( Quaternion, CrossProductIsVectorPartOnly )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    auto tC = belfem::cross( tA, tB );
    EXPECT_NEAR( tC.a(), 0.0, tTol );
}

TEST( Quaternion, CrossProductKnownValues )
{
    // (idea from ChatGPT — verify actual cross product result)
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    auto tC = belfem::cross( tA, tB );
    // cross of (2,3,4) × (6,7,8) = (3*8-4*7, 4*6-2*8, 2*7-3*6) = (-4, 8, -4)
    EXPECT_NEAR( tC.b(), -4.0, tTol );
    EXPECT_NEAR( tC.c(), 8.0,  tTol );
    EXPECT_NEAR( tC.d(), -4.0, tTol );
}

TEST( Quaternion, CrossProductAnticommutative )
{
    belfem::Quaternion< belfem::real > tA( 1.0, 2.0, 3.0, 4.0 );
    belfem::Quaternion< belfem::real > tB( 5.0, 6.0, 7.0, 8.0 );
    auto tAB = belfem::cross( tA, tB );
    auto tBA = belfem::cross( tB, tA );
    EXPECT_NEAR( tAB.b(), -tBA.b(), tTol );
    EXPECT_NEAR( tAB.c(), -tBA.c(), tTol );
    EXPECT_NEAR( tAB.d(), -tBA.d(), tTol );
}

// =============================================================================
// §5.1  Rotation via member function  [semantic]
// =============================================================================

TEST( Quaternion, IdentityRotationIsNoOp )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    belfem::Vector< belfem::real > tV = { 1.0, 2.0, 3.0 };
    auto tR = tQ.rotate( tV );
    EXPECT_NEAR( tR( 0 ), 1.0, tTol );
    EXPECT_NEAR( tR( 1 ), 2.0, tTol );
    EXPECT_NEAR( tR( 2 ), 3.0, tTol );
}

TEST( Quaternion, Rotate90AboutZ )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 2.0 );
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };
    auto tR = tQ.rotate( tV );
    EXPECT_NEAR( tR( 0 ), 0.0, tTol );
    EXPECT_NEAR( tR( 1 ), 1.0, tTol );
    EXPECT_NEAR( tR( 2 ), 0.0, tTol );
}

TEST( Quaternion, Rotate180AboutX )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 0.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI );
    belfem::Vector< belfem::real > tV = { 0.0, 1.0, 0.0 };
    auto tR = tQ.rotate( tV );
    EXPECT_NEAR( tR( 0 ), 0.0,  tTol );
    EXPECT_NEAR( tR( 1 ), -1.0, tTol );
    EXPECT_NEAR( tR( 2 ), 0.0,  tTol );
}

TEST( Quaternion, RotationPreservesNorm )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 1.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 1.23 );
    belfem::Vector< belfem::real > tV = { 3.0, 4.0, 5.0 };
    auto tR = tQ.rotate( tV );
    belfem::real tOrigNorm = std::sqrt( 9.0 + 16.0 + 25.0 );
    belfem::real tRotNorm  = std::sqrt( tR(0)*tR(0) + tR(1)*tR(1) + tR(2)*tR(2) );
    EXPECT_NEAR( tRotNorm, tOrigNorm, tTol );
}

TEST( Quaternion, RotateAndInverseRotateRecovers )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 0.7 );
    belfem::Vector< belfem::real > tV = { 1.0, -2.0, 3.0 };
    auto tRotated = tQ.rotate( tV );
    auto tRecovered = tQ.conj().rotate( tRotated );
    EXPECT_NEAR( tRecovered( 0 ), tV( 0 ), tTol );
    EXPECT_NEAR( tRecovered( 1 ), tV( 1 ), tTol );
    EXPECT_NEAR( tRecovered( 2 ), tV( 2 ), tTol );
}

TEST( Quaternion, CompositionViaMultiplication )
{
    belfem::Vector< belfem::real > tAxis1 = { 1.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tAxis2 = { 0.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis1, 0.5 );
    belfem::Quaternion< belfem::real > tQ2( tAxis2, 0.3 );
    belfem::Vector< belfem::real > tV = { 1.0, 2.0, 3.0 };

    // (q2*q1).rotate(v) should equal q2.rotate(q1.rotate(v))
    auto tComposed = ( tQ2 * tQ1 ).rotate( tV );
    auto tSequential = tQ2.rotate( tQ1.rotate( tV ) );
    EXPECT_NEAR( tComposed( 0 ), tSequential( 0 ), tTol );
    EXPECT_NEAR( tComposed( 1 ), tSequential( 1 ), tTol );
    EXPECT_NEAR( tComposed( 2 ), tSequential( 2 ), tTol );
}

TEST( Quaternion, SignEquivalence )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 1.0 );
    belfem::Quaternion< belfem::real > tNegQ = tQ * -1.0;
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };

    auto tR1 = tQ.rotate( tV );
    auto tR2 = tNegQ.rotate( tV );
    EXPECT_NEAR( tR1( 0 ), tR2( 0 ), tTol );
    EXPECT_NEAR( tR1( 1 ), tR2( 1 ), tTol );
    EXPECT_NEAR( tR1( 2 ), tR2( 2 ), tTol );
}

// =============================================================================
// §5.2  Free function quaternion_rotate_vector  [semantic]
// =============================================================================

TEST( Quaternion, OptimizedMatchesMemberRotate )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 4.0 );
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };

    auto tR1 = tQ.rotate( tV );
    belfem::Vector< belfem::real > tR2( 3, 0.0 );
    belfem::quaternion_rotate_vector( tQ, tV, tR2 );

    EXPECT_NEAR( tR1( 0 ), tR2( 0 ), tTol );
    EXPECT_NEAR( tR1( 1 ), tR2( 1 ), tTol );
    EXPECT_NEAR( tR1( 2 ), tR2( 2 ), tTol );
}

// =============================================================================
// §6  Quaternion ↔ Rotation Matrix Conversion  [semantic]
// =============================================================================

TEST( Quaternion, IdentityQuatGivesIdentityMatrix )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    belfem::Matrix< belfem::real > tR( 3, 3, 0.0 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            EXPECT_NEAR( tR( i, j ), ( i == j ) ? 1.0 : 0.0, tTol );
        }
    }
}

TEST( Quaternion, ResultIsOrthogonal )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 0.7 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );

    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            belfem::real tDot = 0.0;
            for( size_t k = 0; k < 3; ++k )
                tDot += tR( k, i ) * tR( k, j );
            EXPECT_NEAR( tDot, ( i == j ) ? 1.0 : 0.0, tTol );
        }
    }
}

TEST( Quaternion, MatrixRotationMatchesQuaternionRotation )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 3.0 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );

    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };
    auto tViaQ = tQ.rotate( tV );
    belfem::Vector< belfem::real > tViaR( tR * tV );

    EXPECT_NEAR( tViaQ( 0 ), tViaR( 0 ), tTol );
    EXPECT_NEAR( tViaQ( 1 ), tViaR( 1 ), tTol );
    EXPECT_NEAR( tViaQ( 2 ), tViaR( 2 ), tTol );
}

TEST( Quaternion, QuatToMatrixToQuatRoundTrip )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 1.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 1.23 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );
    auto tQ2 = belfem::quaternion_from_rotation_matrix( tR );
    EXPECT_TRUE( same_rotation( tQ, tQ2 ) );
}

TEST( Quaternion, FromRotationMatrixExercisesAllBranches )
{
    // Shepperd's method: 4 branches. Near-180° about each axis makes
    // the corresponding diagonal dominant.

    // Branch 1: trace > 0 (small rotation)
    {
        belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
        belfem::Quaternion< belfem::real > tQ( tAxis, 0.1 );
        belfem::Matrix< belfem::real > tR( 3, 3 );
        belfem::quaternion_to_rotation_matrix( tQ, tR );
        auto tQ2 = belfem::quaternion_from_rotation_matrix( tR );
        EXPECT_TRUE( same_rotation( tQ, tQ2 ) );
    }
    // Branch 2: r00 dominant (near-180° about x)
    {
        belfem::Vector< belfem::real > tAxis = { 1.0, 0.0, 0.0 };
        belfem::Quaternion< belfem::real > tQ( tAxis, M_PI - 0.01 );
        belfem::Matrix< belfem::real > tR( 3, 3 );
        belfem::quaternion_to_rotation_matrix( tQ, tR );
        auto tQ2 = belfem::quaternion_from_rotation_matrix( tR );
        EXPECT_TRUE( same_rotation( tQ, tQ2 ) );
    }
    // Branch 3: r11 dominant (near-180° about y)
    {
        belfem::Vector< belfem::real > tAxis = { 0.0, 1.0, 0.0 };
        belfem::Quaternion< belfem::real > tQ( tAxis, M_PI - 0.01 );
        belfem::Matrix< belfem::real > tR( 3, 3 );
        belfem::quaternion_to_rotation_matrix( tQ, tR );
        auto tQ2 = belfem::quaternion_from_rotation_matrix( tR );
        EXPECT_TRUE( same_rotation( tQ, tQ2 ) );
    }
    // Branch 4: r22 dominant (near-180° about z)
    {
        belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
        belfem::Quaternion< belfem::real > tQ( tAxis, M_PI - 0.01 );
        belfem::Matrix< belfem::real > tR( 3, 3 );
        belfem::quaternion_to_rotation_matrix( tQ, tR );
        auto tQ2 = belfem::quaternion_from_rotation_matrix( tR );
        EXPECT_TRUE( same_rotation( tQ, tQ2 ) );
    }
}

// =============================================================================
// §7  Slerp  [semantic]
// =============================================================================

TEST( Quaternion, SlerpAtZeroReturnsFirst )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis, 0.0 );
    belfem::Quaternion< belfem::real > tQ2( tAxis, M_PI / 2.0 );
    auto tResult = belfem::slerp( tQ1, tQ2, 0.0 );
    EXPECT_TRUE( same_rotation( tResult, tQ1 ) );
}

TEST( Quaternion, SlerpAtOneReturnsSecond )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis, 0.0 );
    belfem::Quaternion< belfem::real > tQ2( tAxis, M_PI / 2.0 );
    auto tResult = belfem::slerp( tQ1, tQ2, 1.0 );
    EXPECT_TRUE( same_rotation( tResult, tQ2 ) );
}

TEST( Quaternion, SlerpMidpoint90Deg )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis, 0.0 );
    belfem::Quaternion< belfem::real > tQ2( tAxis, M_PI / 2.0 );
    auto tMid = belfem::slerp( tQ1, tQ2, 0.5 );

    belfem::real tHalf = M_PI / 8.0;
    EXPECT_NEAR( tMid.a(), std::cos( tHalf ), tTol );
    EXPECT_NEAR( tMid.d(), std::sin( tHalf ), tTol );
}

TEST( Quaternion, SlerpSameQuaternion )
{
    // (from test plan §7.1 — noticed by Grok)
    belfem::Vector< belfem::real > tAxis = { 1.0, 0.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 0.5 );
    auto tResult = belfem::slerp( tQ, tQ, 0.5 );
    EXPECT_TRUE( same_rotation( tResult, tQ ) );
}

TEST( Quaternion, SlerpResultIsUnit )
{
    belfem::Vector< belfem::real > tAxis1 = { 1.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tAxis2 = { 0.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis1, 0.3 );
    belfem::Quaternion< belfem::real > tQ2( tAxis2, 1.7 );
    for( int i = 0; i <= 10; ++i )
    {
        belfem::real tT = static_cast< belfem::real >( i ) / 10.0;
        auto tResult = belfem::slerp( tQ1, tQ2, tT );
        EXPECT_NEAR( tResult.norm(), 1.0, tTol );
    }
}

TEST( Quaternion, SlerpShortestPath )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis, 0.0 );
    belfem::Quaternion< belfem::real > tQ2( tAxis, M_PI / 2.0 );
    belfem::Quaternion< belfem::real > tNegQ2 = tQ2 * -1.0;

    auto tMid1 = belfem::slerp( tQ1, tQ2, 0.5 );
    auto tMid2 = belfem::slerp( tQ1, tNegQ2, 0.5 );

    EXPECT_TRUE( same_rotation( tMid1, tMid2 ) );
}

TEST( Quaternion, SlerpNearlyParallelFallback )
{
    // when q1 ≈ q2, slerp uses linear interpolation fallback
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ1( tAxis, 0.0 );
    belfem::Quaternion< belfem::real > tQ2( tAxis, 1e-14 );  // nearly identical
    auto tMid = belfem::slerp( tQ1, tQ2, 0.5 );
    EXPECT_NEAR( tMid.norm(), 1.0, tTol );
    EXPECT_TRUE( same_rotation( tMid, tQ1 ) );
}

// =============================================================================
// Missing coverage rows (Codex finding #4)
// =============================================================================

TEST( Quaternion, Rotate120AboutDiagonal )
{
    // 120° about (1,1,1)/√3: cyclic permutation x→y→z→x
    belfem::Vector< belfem::real > tAxis = { 1.0, 1.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 2.0 * M_PI / 3.0 );
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };   // x-axis
    auto tR = tQ.rotate( tV );
    EXPECT_NEAR( tR( 0 ), 0.0, tTol );
    EXPECT_NEAR( tR( 1 ), 1.0, tTol );
    EXPECT_NEAR( tR( 2 ), 0.0, tTol );
}

TEST( Quaternion, OptimizedPreservesNorm )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 1.23 );
    belfem::Vector< belfem::real > tV = { 3.0, 4.0, 5.0 };
    belfem::Vector< belfem::real > tOut( 3, 0.0 );
    belfem::quaternion_rotate_vector( tQ, tV, tOut );
    belfem::real tOrigNorm = std::sqrt( 9.0 + 16.0 + 25.0 );
    belfem::real tRotNorm  = std::sqrt( tOut(0)*tOut(0) + tOut(1)*tOut(1) + tOut(2)*tOut(2) );
    EXPECT_NEAR( tRotNorm, tOrigNorm, tTol );
}

TEST( Quaternion, ResultHasDetPlusOne )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 0.7 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );

    // det(R) = +1 for proper rotation
    belfem::real tDet = tR(0,0) * ( tR(1,1)*tR(2,2) - tR(1,2)*tR(2,1) )
                      - tR(0,1) * ( tR(1,0)*tR(2,2) - tR(1,2)*tR(2,0) )
                      + tR(0,2) * ( tR(1,0)*tR(2,1) - tR(1,1)*tR(2,0) );
    EXPECT_NEAR( tDet, 1.0, tTol );
}

TEST( Quaternion, NinetyDegAboutZMatrix )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 1.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 2.0 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );

    // expected: [[0,-1,0],[1,0,0],[0,0,1]]
    EXPECT_NEAR( tR( 0, 0 ),  0.0, tTol );
    EXPECT_NEAR( tR( 0, 1 ), -1.0, tTol );
    EXPECT_NEAR( tR( 1, 0 ),  1.0, tTol );
    EXPECT_NEAR( tR( 1, 1 ),  0.0, tTol );
    EXPECT_NEAR( tR( 2, 2 ),  1.0, tTol );
}

TEST( Quaternion, IdentityMatrixGivesIdentityQuat )
{
    belfem::Matrix< belfem::real > tR( 3, 3, 0.0 );
    tR( 0, 0 ) = 1.0;
    tR( 1, 1 ) = 1.0;
    tR( 2, 2 ) = 1.0;
    auto tQ = belfem::quaternion_from_rotation_matrix( tR );
    auto tId = belfem::Quaternion< belfem::real >::identity();
    EXPECT_TRUE( same_rotation( tQ, tId ) );
}

TEST( Quaternion, FromRotationMatrixResultIsUnit )
{
    belfem::Vector< belfem::real > tAxis = { 1.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 1.5 );
    belfem::Matrix< belfem::real > tR( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR );
    auto tQ2 = belfem::quaternion_from_rotation_matrix( tR );
    EXPECT_NEAR( tQ2.norm(), 1.0, tEps );
}

TEST( Quaternion, MatrixToQuatToMatrixRoundTrip )
{
    // start with a known rotation matrix, convert to quaternion, convert back
    belfem::Vector< belfem::real > tAxis = { 0.0, 1.0, 0.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, M_PI / 3.0 );
    belfem::Matrix< belfem::real > tR1( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tR1 );

    auto tQ2 = belfem::quaternion_from_rotation_matrix( tR1 );
    belfem::Matrix< belfem::real > tR2( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ2, tR2 );

    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            EXPECT_NEAR( tR1( i, j ), tR2( i, j ), tTol );
        }
    }
}

TEST( Quaternion, AllThreeRotationPathsAgree )
{
    // member rotate(), free function, and matrix multiply all produce same result
    belfem::Vector< belfem::real > tAxis = { 1.0, 2.0, 3.0 };
    belfem::Quaternion< belfem::real > tQ( tAxis, 0.7 );
    belfem::Vector< belfem::real > tV = { 1.0, -2.0, 3.0 };

    // path 1: member rotate
    auto tR1 = tQ.rotate( tV );

    // path 2: optimized free function
    belfem::Vector< belfem::real > tR2( 3, 0.0 );
    belfem::quaternion_rotate_vector( tQ, tV, tR2 );

    // path 3: rotation matrix multiply
    belfem::Matrix< belfem::real > tMat( 3, 3 );
    belfem::quaternion_to_rotation_matrix( tQ, tMat );
    belfem::Vector< belfem::real > tR3( tMat * tV );

    for( size_t i = 0; i < 3; ++i )
    {
        EXPECT_NEAR( tR1( i ), tR2( i ), tTol );
        EXPECT_NEAR( tR1( i ), tR3( i ), tTol );
    }
}

// =============================================================================
// Debug tests  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( QuaternionDebug, FromVectorWrongLengthThrows )
{
    belfem::Vector< belfem::real > tV( 2, 1.0 );
    EXPECT_THROW( belfem::Quaternion< belfem::real > tQ( tV ), std::runtime_error );
}

TEST( QuaternionDebug, AxisAngleWrongLengthThrows )
{
    belfem::Vector< belfem::real > tAxis( 5, 1.0 );
    EXPECT_THROW(
        ( belfem::Quaternion< belfem::real >( tAxis, 1.0 ) ),
        std::runtime_error );
}

TEST( QuaternionDebug, AxisAngleZeroAxisThrows )
{
    belfem::Vector< belfem::real > tAxis = { 0.0, 0.0, 0.0 };
    EXPECT_THROW(
        ( belfem::Quaternion< belfem::real >( tAxis, 1.0 ) ),
        std::runtime_error );
}

TEST( QuaternionDebug, InitializerListWrongSizeThrows )
{
    belfem::Quaternion< belfem::real > tQ;
    EXPECT_THROW( ( tQ = { 1.0, 2.0, 3.0 } ), std::runtime_error );
}

TEST( QuaternionDebug, DivideByZeroThrows )
{
    // /= uses BELFEM_ERROR — always active, but only throws in debug
    // (idea from ChatGPT)
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    EXPECT_THROW( tQ /= 0.0, std::runtime_error );
}

TEST( QuaternionDebug, InvZeroThrows )
{
    // inv() uses BELFEM_ERROR (idea from ChatGPT)
    belfem::Quaternion< belfem::real > tQ;
    EXPECT_THROW( tQ.inv(), std::runtime_error );
}

TEST( QuaternionDebug, NormalizeZeroThrows )
{
    // normalize() uses BELFEM_ERROR (idea from ChatGPT)
    belfem::Quaternion< belfem::real > tQ;
    EXPECT_THROW( tQ.normalize(), std::runtime_error );
}

TEST( QuaternionDebug, RotateNonUnitThrows )
{
    belfem::Quaternion< belfem::real > tQ( 1.0, 2.0, 3.0, 4.0 );
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };
    EXPECT_THROW( tQ.rotate( tV ), std::runtime_error );
}

TEST( QuaternionDebug, RotateWrongVectorLengthThrows )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    belfem::Vector< belfem::real > tV( 2, 1.0 );
    EXPECT_THROW( tQ.rotate( tV ), std::runtime_error );
}

TEST( QuaternionDebug, OptimizedNonUnitThrows )
{
    belfem::Quaternion< belfem::real > tQ( 2.0, 0.0, 0.0, 0.0 );
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tOut( 3, 0.0 );
    EXPECT_THROW( belfem::quaternion_rotate_vector( tQ, tV, tOut ),
                  std::runtime_error );
}

TEST( QuaternionDebug, OptimizedInputLengthThrows )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    belfem::Vector< belfem::real > tV( 2, 1.0 );
    belfem::Vector< belfem::real > tOut( 3, 0.0 );
    EXPECT_THROW( belfem::quaternion_rotate_vector( tQ, tV, tOut ),
                  std::runtime_error );
}

TEST( QuaternionDebug, OptimizedOutputLengthThrows )
{
    // (idea from ChatGPT)
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tOut( 2, 0.0 );
    EXPECT_THROW( belfem::quaternion_rotate_vector( tQ, tV, tOut ),
                  std::runtime_error );
}

TEST( QuaternionDebug, ToMatrixNon3x3Throws )
{
    auto tQ = belfem::Quaternion< belfem::real >::identity();
    belfem::Matrix< belfem::real > tR( 2, 2, 0.0 );
    EXPECT_THROW( belfem::quaternion_to_rotation_matrix( tQ, tR ),
                  std::runtime_error );
}

TEST( QuaternionDebug, ToMatrixNonUnitQuatThrows )
{
    // (idea from ChatGPT)
    belfem::Quaternion< belfem::real > tQ( 2.0, 0.0, 0.0, 0.0 );
    belfem::Matrix< belfem::real > tR( 3, 3, 0.0 );
    EXPECT_THROW( belfem::quaternion_to_rotation_matrix( tQ, tR ),
                  std::runtime_error );
}

TEST( QuaternionDebug, FromMatrixNon3x3Throws )
{
    belfem::Matrix< belfem::real > tR( 2, 3, 0.0 );
    EXPECT_THROW( belfem::quaternion_from_rotation_matrix( tR ),
                  std::runtime_error );
}

TEST( QuaternionDebug, FromMatrixBadDetThrows )
{
    // matrix with det ≠ 1 (idea from ChatGPT)
    belfem::Matrix< belfem::real > tR = { { 2.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 1.0 } };
    EXPECT_THROW( belfem::quaternion_from_rotation_matrix( tR ),
                  std::runtime_error );
}

TEST( QuaternionDebug, SlerpNonUnitFirstThrows )
{
    belfem::Quaternion< belfem::real > tQ1( 2.0, 0.0, 0.0, 0.0 );
    auto tQ2 = belfem::Quaternion< belfem::real >::identity();
    EXPECT_THROW( belfem::slerp( tQ1, tQ2, 0.5 ), std::runtime_error );
}

TEST( QuaternionDebug, SlerpNonUnitSecondThrows )
{
    auto tQ1 = belfem::Quaternion< belfem::real >::identity();
    belfem::Quaternion< belfem::real > tQ2( 2.0, 0.0, 0.0, 0.0 );
    EXPECT_THROW( belfem::slerp( tQ1, tQ2, 0.5 ), std::runtime_error );
}

#endif // NDEBUG
