/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Solver facade, SolverParameters, and backend smoke tests.
 * See: tests_08_sparse.md §5–§7
 *
 * Backend tests are wrapped in #ifdef guards. Each uses the canonical
 * known-solution pattern: build tridiagonal K, compute b = K*x_known,
 * solve K*x = b, verify x ≈ x_known.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Solver.hpp"
#include <limits>
#include "cl_SolverParameters.hpp"
#include "en_SolverEnums.hpp"
#include "mumpstools.hpp"

namespace
{
    const belfem::real tEps = 1e-12;
    const belfem::real tTol = 1e-9;

    // Build an N×N tridiagonal Laplacian: 2 on diagonal, -1 off-diagonal
    belfem::Matrix< belfem::real > make_tridiag( belfem::uint aN )
    {
        belfem::Matrix< belfem::real > tK( aN, aN, 0.0 );
        for( belfem::uint k = 0; k < aN; ++k )
        {
            tK( k, k ) = 2.0;
            if( k > 0 )     tK( k, k - 1 ) = -1.0;
            if( k < aN - 1 ) tK( k, k + 1 ) = -1.0;
        }
        return tK;
    }

    // Canonical known-solution solve helper.
    // Builds tridiagonal K, computes b = K * x_known, solves, verifies.
    void solve_tridiag_and_verify(
        belfem::SolverType aType,
        belfem::uint aN = 6 )
    {
        belfem::Matrix< belfem::real > tK = make_tridiag( aN );
        belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

        belfem::Vector< belfem::real > tXknown( aN );
        for( belfem::uint k = 0; k < aN; ++k )
        {
            tXknown( k ) = static_cast< belfem::real >( k + 1 );
        }

        // b = K * x_known (dense multiply)
        belfem::Vector< belfem::real > tB( tK * tXknown );

        // solve
        belfem::Solver tSolver( aType );
        belfem::Vector< belfem::real > tX( aN, 0.0 );
        tSolver.solve( tM, tX, tB );

        // verify
        for( belfem::uint k = 0; k < aN; ++k )
        {
            EXPECT_NEAR( tX( k ), tXknown( k ), tTol );
        }
    }
}

// =============================================================================
// §5.1  Solver Construction  [semantic]
// =============================================================================

#ifdef BELFEM_SUITESPARSE

TEST( SolverConstruction, SolverTypeStored )
{
    belfem::Solver tSolver( belfem::SolverType::UMFPACK );
    EXPECT_EQ( tSolver.type(), belfem::SolverType::UMFPACK );
}

TEST( SolverConstruction, SolverWrapperNotNull )
{
    belfem::Solver tSolver( belfem::SolverType::UMFPACK );
    EXPECT_NE( tSolver.wrapper(), nullptr );
}

TEST( SolverConstruction, SolverFromParameters )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    belfem::Solver tSolver( tParams );
    EXPECT_EQ( tSolver.type(), belfem::SolverType::UMFPACK );
}

#endif // BELFEM_SUITESPARSE

// =============================================================================
// §5.2  Known-Solution Solve  [semantic]
// =============================================================================

#ifdef BELFEM_SUITESPARSE

TEST( SolverSolve, UMFPACKSolveTridiagonal )
{
    solve_tridiag_and_verify( belfem::SolverType::UMFPACK );
}

#endif

#ifdef BELFEM_MUMPS

TEST( SolverSolve, MUMPSSolveTridiagonal )
{
    solve_tridiag_and_verify( belfem::SolverType::MUMPS );
}

// The MUMPS symmetric-mode guard, always-active.
//
// Under SYM != 0 MUMPS wants exactly ONE representative of each symmetric
// coordinate -- either triangle will do, it does not insist on the lower one --
// and duplicate entries are SUMMED, so if both a(i,j) and a(j,i) are supplied
// they are added together ( MUMPS 5.5.1 user guide, section 5.2.2.1 ). BELFEM
// assembles and hands over the FULL matrix and nothing extracts a triangle, so
// a symmetric mode would double every off-diagonal and factorize a different
// operator -- WITHOUT FAILING. Measured 2026-08-29 on a matrix with an analytic
// spectrum: ||Ax-b||/||b|| = 7.4e16 and a converged eigenvalue of -2.5e-19
// against a true 2.46e-6, with the eigensolver reporting success throughout.
//
// The wrapper therefore refuses SYM != 0 outright. This test exists because
// NOTHING IN THE TREE ASKS FOR A SYMMETRIC MODE -- the IWG default is
// Unsymmetric and both direct derivers pass or inherit it -- so the guard has
// no production path that would exercise it, and "no path reaches it" is an
// argument rather than evidence. The test turns it into evidence, and it will
// also catch a future default flipping back.
TEST( SolverSymmetry, MUMPSRejectsSymmetricModes )
{
    // The matrix here is genuinely symmetric, which is the point: the guard
    // must refuse the MODE, not the matrix. A caller with a symmetric operator
    // is exactly who would reach for SYM != 0 and get a silently wrong answer
    const belfem::uint tN = 6 ;

    belfem::Matrix< belfem::real > tKdense = make_tridiag( tN );

    // both symmetric modes must be refused, and refused the SAME way -- a
    // downgrade to SYM = 0 would be a silent lie about what was solved
    {
        belfem::SpMatrix tK( tKdense, belfem::SpMatrixType::CSR );
        belfem::Vector< belfem::real > tX( tN, 0.0 );
        belfem::Vector< belfem::real > tB( tN, 1.0 );

        belfem::Solver tSolver( belfem::SolverType::MUMPS );
        tSolver.set_symmetry_mode( belfem::SymmetryMode::GeneralSymmetric );

        EXPECT_THROW( tSolver.solve( tK, tX, tB ), std::runtime_error );
    }

    {
        belfem::SpMatrix tK( tKdense, belfem::SpMatrixType::CSR );
        belfem::Vector< belfem::real > tX( tN, 0.0 );
        belfem::Vector< belfem::real > tB( tN, 1.0 );

        belfem::Solver tSolver( belfem::SolverType::MUMPS );
        tSolver.set_symmetry_mode( belfem::SymmetryMode::PositiveDefiniteSymmetric );

        EXPECT_THROW( tSolver.solve( tK, tX, tB ), std::runtime_error );
    }
}

// the companion: the mode BELFEM does support must still work, so the guard
// cannot be "passing" by refusing everything
TEST( SolverSymmetry, MUMPSAcceptsUnsymmetric )
{
    solve_tridiag_and_verify( belfem::SolverType::MUMPS );
}

// The mumpstools instance registry must RECYCLE freed slots. The old
// allocator handed out monotonically increasing IDs and only reset once
// EVERY instance was freed, so one long-lived solver plus a consumer that
// creates and frees per call -- the production-solver + per-timestep
// conditioning-diagnostic pattern -- exhausted the default pool of 8 after
// seven cycles ( the tapestack3d step-8 abort, 2026-08-29 ). The churn
// below drives that exact pattern past the pool size: it aborts on the
// monotonic allocator and passes on the first-free scan. This is the only
// executable gate for the allocator -- the production consumer now keeps
// its instance alive, so no solver run exercises recycling.
TEST( SolverLifecycle, MUMPSPoolRecyclesFreedSlots )
{
    const belfem::uint tN = 6 ;

    belfem::Matrix< belfem::real > tK = make_tridiag( tN );

    belfem::Vector< belfem::real > tXknown( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tXknown( k ) = static_cast< belfem::real >( k + 1 );
    }
    belfem::Vector< belfem::real > tB( tK * tXknown );

    // pin one instance for the whole test so the pool can never fully
    // drain -- without this, every free would reset the old allocator and
    // hide the defect
    belfem::SpMatrix tMpin( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tXpin( tN, 0.0 );
    belfem::Solver tPinned( belfem::SolverType::MUMPS );
    tPinned.solve( tMpin, tXpin, tB );

    // churn: init + free more times than the pool has slots ( 8 )
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Solver tChurn( belfem::SolverType::MUMPS );

    for( belfem::uint k = 0; k < 12; ++k )
    {
        belfem::Vector< belfem::real > tX( tN, 0.0 );

        tChurn.solve( tM, tX, tB );   // lazy init creates an instance

        for( belfem::uint i = 0; i < tN; ++i )
        {
            EXPECT_NEAR( tX( i ), tXknown( i ), tTol );
        }

        tChurn.free();                // ... and this releases it again
    }

    // The pinned instance must still be intact. Without this the test passes
    // even if a churn free() had destroyed the PIN's slot: first-free would
    // simply hand that slot back to tChurn, all twelve solves would succeed,
    // and the cross-wrapper destroy would go unnoticed. Solving through the
    // pin after the churn is what encodes that hazard rather than describing
    // it in a comment ( flagged by the code-audit round, 2026-08-29 )
    belfem::Vector< belfem::real > tXpin2( tN, 0.0 );
    tPinned.solve( tMpin, tXpin2, tB );

    for( belfem::uint i = 0; i < tN; ++i )
    {
        EXPECT_NEAR( tXpin2( i ), tXknown( i ), tTol );
    }
}

// An instance that is created and then freed WITHOUT EVER SOLVING must still
// give its slot back. Until 2026-08-29 it did not: the Fortran slot is
// reserved by the create, but MUMPS::mInitialized -- the only thing free()
// consults before issuing JOB = -2 -- was set by a SOLVE. So initialize +
// free leaked the slot, and so did an instance driven only through the
// matrix-RHS overload, which never set the flag on success at all ( two
// symptoms of one defect, closed by one fix ).
//
// Neither is reachable through Solver::solve, whose init is lazy and always
// followed by a solve -- but Solver::wrapper() is public and the shift-invert
// conditioning consumer already uses it, so the path is one public call away.
// That is also what makes this test possible without new API.
//
// Twelve cycles against a pool of eight: on the old code the first eight leak
// a slot each and the NINTH create finds the registry full, which aborts in
// MUMPS::initialize because soft-fail is off by default. Twelve, not eight --
// with no pinned occupant the eighth create still succeeds.
//
// This does NOT replace MUMPSPoolRecyclesFreedSlots: that one gates the
// first-free allocator on the solve+free path. Different gap, different path.
//
// STILL NOT COVERED, and deliberately said out loud: the matrix-RHS overload,
// which is the other half of the same defect ( it never set the flag on a
// SUCCESSFUL solve ). The fix is shared -- initialize() sets the flag for both
// overloads -- but no in-tree caller solves matrix-RHS-only, so there is
// nothing to drive it from here without inventing a consumer. Nor is any of
// this MPI: the exhaustion barrier and mixed-occupancy paths have no gate.
TEST( SolverLifecycle, MUMPSFreesAnInstanceThatNeverSolved )
{
    const belfem::uint tN = 6 ;

    belfem::Matrix< belfem::real > tK = make_tridiag( tN );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

    for( belfem::uint k = 0; k < 12; ++k )
    {
        belfem::Solver tSolver( belfem::SolverType::MUMPS );

        tSolver.wrapper()->initialize( tM, belfem::SymmetryMode::Unsymmetric, 1 );

        EXPECT_TRUE( tSolver.wrapper()->is_initialized() );

        // the release the old code skipped
        tSolver.free();

        EXPECT_FALSE( tSolver.wrapper()->is_initialized() );
    }

    // and the pool is still usable afterwards -- a test that only proved
    // "no abort" would also pass if the registry had been left in a state
    // where nothing can be created at all
    solve_tridiag_and_verify( belfem::SolverType::MUMPS );
}

#endif

#ifdef BELFEM_PARDISO

TEST( SolverSolve, PARDISOSolveTridiagonal )
{
    solve_tridiag_and_verify( belfem::SolverType::PARDISO );
}

#endif

// NOTE: PETSc and STRUMPACK may need MPI init. If they fail here,
// move these tests to a Tier 2 MPI test binary.

#ifdef BELFEM_PETSC

TEST( SolverSolve, PETScSolveTridiagonal )
{
    solve_tridiag_and_verify( belfem::SolverType::PETSc );
}

#endif

#ifdef BELFEM_STRUMPACK

TEST( SolverSolve, STRUMPACKSolveTridiagonal )
{
    solve_tridiag_and_verify( belfem::SolverType::STRUMPACK );
}

#endif

// =============================================================================
// §5.3  Re-Solve with Changed Values  [semantic]
// =============================================================================

#ifdef BELFEM_SUITESPARSE

TEST( SolverSolve, ReSolveReuseFactorization )
{
    belfem::uint tN = 6;
    belfem::Matrix< belfem::real > tK = make_tridiag( tN );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

    belfem::Vector< belfem::real > tXknown( tN );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tXknown( k ) = static_cast< belfem::real >( k + 1 );
    }

    belfem::Vector< belfem::real > tB( tK * tXknown );
    belfem::Vector< belfem::real > tX( tN, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::UMFPACK );
    tSolver.solve( tM, tX, tB );

    // verify first solve
    for( belfem::uint k = 0; k < tN; ++k )
    {
        EXPECT_NEAR( tX( k ), tXknown( k ), tTol );
    }

    // change values: scale diagonal by 2 (same structure)
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tM( k, k ) = 4.0;
    }

    // recompute RHS with new matrix values
    belfem::Matrix< belfem::real > tK2 = tK;
    for( belfem::uint k = 0; k < tN; ++k ) tK2( k, k ) = 4.0;

    belfem::Vector< belfem::real > tB2( tK2 * tXknown );
    belfem::Vector< belfem::real > tX2( tN, 0.0 );

    // re-solve (wrapper should detect already initialized, update values)
    tSolver.solve( tM, tX2, tB2 );

    for( belfem::uint k = 0; k < tN; ++k )
    {
        EXPECT_NEAR( tX2( k ), tXknown( k ), tTol );
    }
}

#endif // BELFEM_SUITESPARSE

// =============================================================================
// §5.4  Solver Lifecycle  [semantic]
// =============================================================================

#ifdef BELFEM_SUITESPARSE

TEST( SolverLifecycle, LazyInitialization )
{
    belfem::Solver tSolver( belfem::SolverType::UMFPACK );

    // before first solve, wrapper should not be initialized
    EXPECT_FALSE( tSolver.wrapper()->is_initialized() );

    // solve
    belfem::uint tN = 4;
    belfem::Matrix< belfem::real > tK = make_tridiag( tN );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tX( tN, 0.0 );
    belfem::Vector< belfem::real > tB( tN, 1.0 );

    tSolver.solve( tM, tX, tB );

    EXPECT_TRUE( tSolver.wrapper()->is_initialized() );
}

TEST( SolverLifecycle, FreeResetsState )
{
    belfem::Solver tSolver( belfem::SolverType::UMFPACK );

    belfem::uint tN = 4;
    belfem::Matrix< belfem::real > tK = make_tridiag( tN );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tX( tN, 0.0 );
    belfem::Vector< belfem::real > tB( tN, 1.0 );

    tSolver.solve( tM, tX, tB );
    ASSERT_TRUE( tSolver.wrapper()->is_initialized() );

    tSolver.free();

    EXPECT_FALSE( tSolver.wrapper()->is_initialized() );
}

#endif // BELFEM_SUITESPARSE

// =============================================================================
// §6  SolverParameters  [semantic]
// =============================================================================

TEST( SolverParameters, ReorderingMethodStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_reordering_method( belfem::ReorderingMethod::NATURAL );
    EXPECT_EQ( tParams.reordering_method(), belfem::ReorderingMethod::NATURAL );
}

TEST( SolverParameters, CompressionMethodStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_compression_method( belfem::CompressionMethod::BLR );
    EXPECT_EQ( tParams.compression_method(), belfem::CompressionMethod::BLR );
}

TEST( SolverParameters, InitialGuessFlagStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_use_initial_guess( true );
    EXPECT_TRUE( tParams.use_initial_guess() );
}

TEST( SolverParameters, RelativeToleranceStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_relative_tolerance( 1e-6 );
    EXPECT_NEAR( tParams.relative_tolerance(), 1e-6, tEps );
}

TEST( SolverParameters, DefaultParametersValid )
{
    // (idea from ChatGPT) — verify all defaults match source
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );

    EXPECT_EQ( tParams.type(), belfem::SolverType::UMFPACK );
    EXPECT_EQ( tParams.distributed_matrix_type(), belfem::DistributedMatrixType::CSR );
    EXPECT_EQ( tParams.reordering_method(), belfem::ReorderingMethod::AUTOMATIC );
    EXPECT_EQ( tParams.compression_method(), belfem::CompressionMethod::AUTOMATIC );
    EXPECT_EQ( tParams.krylov_method(), belfem::KrylovMethod::AUTO );

    // pins the source default in cl_SolverParameters.hpp — moved
    // 1e-6 -> 1e-8 -> 1e-10 over the solver campaigns; a deliberate
    // default change updates this line in the same commit
    EXPECT_NEAR( tParams.relative_tolerance(), 1e-10, tEps );

    EXPECT_FALSE( tParams.use_initial_guess() );

    // the cutoff default and its unset have-flag
    EXPECT_NEAR( tParams.compression_cutoff(), 1e-8, tEps );
    EXPECT_FALSE( tParams.have_compression_cutoff() );
}

TEST( SolverParameters, CompressionCutoffStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_compression_cutoff( 1e-12 );
    EXPECT_NEAR( tParams.compression_cutoff(), 1e-12, tEps );
    EXPECT_TRUE( tParams.have_compression_cutoff() );
}

TEST( SolverParameters, AbsoluteToleranceStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );

    // defaults: 1e-14 replaces STRUMPACK's 1e-10 library floor;
    // the have-flag stays false so PETSc keeps PETSC_DEFAULT
    EXPECT_NEAR( tParams.absolute_tolerance(), 1e-14, tEps );
    EXPECT_FALSE( tParams.have_absolute_tolerance() );

    tParams.set_absolute_tolerance( 1e-12 );
    EXPECT_NEAR( tParams.absolute_tolerance(), 1e-12, tEps );
    EXPECT_TRUE( tParams.have_absolute_tolerance() );
}

TEST( SolverParameters, AbsoluteToleranceRejectsInvalid )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    EXPECT_THROW( tParams.set_absolute_tolerance( 0.0 ), std::runtime_error );
    EXPECT_THROW( tParams.set_absolute_tolerance( -1e-14 ), std::runtime_error );
    EXPECT_THROW( tParams.set_absolute_tolerance(
        std::numeric_limits< belfem::real >::infinity() ), std::runtime_error );
    EXPECT_THROW( tParams.set_absolute_tolerance(
        std::numeric_limits< belfem::real >::quiet_NaN() ), std::runtime_error );

    // a rejected value must not have touched the state
    EXPECT_NEAR( tParams.absolute_tolerance(), 1e-14, tEps );
    EXPECT_FALSE( tParams.have_absolute_tolerance() );
}

TEST( SolverParameters, CompressionCutoffRejectsInvalid )
{
    // the setter is the single validation path ( the deck parser calls
    // it ): zero is MUMPS's lossless-BLR niche, deliberately unreachable;
    // negative and non-finite values must not reach rel_tol / CNTL(7)
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    EXPECT_THROW( tParams.set_compression_cutoff( 0.0 ), std::runtime_error );
    EXPECT_THROW( tParams.set_compression_cutoff( -1e-8 ), std::runtime_error );
    EXPECT_THROW( tParams.set_compression_cutoff(
        std::numeric_limits< belfem::real >::infinity() ), std::runtime_error );
    EXPECT_THROW( tParams.set_compression_cutoff(
        std::numeric_limits< belfem::real >::quiet_NaN() ), std::runtime_error );

    // a rejected value must not have touched the state
    EXPECT_NEAR( tParams.compression_cutoff(), 1e-8, tEps );
    EXPECT_FALSE( tParams.have_compression_cutoff() );
}

TEST( SolverParameters, MemoryBudgetStored )
{
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );

    // unset: 0 and no have-flag, which the MUMPS wrapper reads as
    // "measure the machine on the first workspace failure"
    EXPECT_EQ( tParams.memory_budget(), 0u );
    EXPECT_FALSE( tParams.have_memory_budget() );

    tParams.set_memory_budget( 4000 );
    EXPECT_EQ( tParams.memory_budget(), 4000u );
    EXPECT_TRUE( tParams.have_memory_budget() );
}

TEST( SolverParameters, MemoryBudgetRejectsZero )
{
    // 0 is the unset value; stating it would silently mean "no cap"
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    EXPECT_ANY_THROW( tParams.set_memory_budget( 0 ) );
}

TEST( SolverParameters, MemoryBudgetRejectsWhatTheShimCannotCarry )
{
    // the value reaches MUMPS as int_t; past its maximum the narrowing
    // turns negative and the guarded Fortran write would silently skip
    // ICNTL(23). The setter refuses the boundary; the largest value that
    // fits is accepted
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    const belfem::uint tMax = ( belfem::uint ) std::numeric_limits< belfem::int_t >::max();
    EXPECT_NO_THROW( tParams.set_memory_budget( tMax ) );
    EXPECT_EQ( tParams.memory_budget(), tMax );
#ifndef BELFEM_INT64
    EXPECT_ANY_THROW( tParams.set_memory_budget( tMax + 1u ) );
#endif
}

// the retry policy behind MUMPS::escalate_workspace(), pure so the table
// is checkable without a MUMPS instance. Column order:
// ( INFOG(1), ICNTL(14), ceiling, ICNTL(23) slot, budget MB, estimate MB )
TEST( SolverPolicy, MumpsWorkspaceActionTable )
{
    using mumps::WorkspaceAction;
    using mumps::next_workspace_action;

    // a code that is not a workspace failure is never retried
    EXPECT_EQ( next_workspace_action( -10, 30, 480, 0, 2000, 300 ), WorkspaceAction::GiveUp );
    EXPECT_EQ( next_workspace_action(   0, 30, 480, 0, 2000, 300 ), WorkspaceAction::GiveUp );

    // -19: the cap cannot be met, whatever else is true
    EXPECT_EQ( next_workspace_action( -19, 30, 480, 2000, 2000, 300 ), WorkspaceAction::GiveUp );
    EXPECT_EQ( next_workspace_action( -19, 30, 480, 0,    2000, 300 ), WorkspaceAction::GiveUp );

    // first -9 with a known budget above the estimate: cap
    EXPECT_EQ( next_workspace_action( -9, 30, 480, 0, 2000, 300 ), WorkspaceAction::Cap );
    // -8 takes the same route
    EXPECT_EQ( next_workspace_action( -8, 30, 480, 0, 2000, 300 ), WorkspaceAction::Cap );
    // an unknown estimate ( 0 ) does not block the cap
    EXPECT_EQ( next_workspace_action( -9, 30, 480, 0, 2000, 0 ), WorkspaceAction::Cap );
    // budget exactly at the estimate is still a cap
    EXPECT_EQ( next_workspace_action( -9, 30, 480, 0, 300, 300 ), WorkspaceAction::Cap );

    // known budget BELOW the estimate: the machine cannot hold it, and a
    // rung would ask for more -- give up, never ladder
    EXPECT_EQ( next_workspace_action( -9, 30, 480, 0, 200, 300 ), WorkspaceAction::GiveUp );

    // cap already set: the ladder, up to the ceiling
    EXPECT_EQ( next_workspace_action( -9, 30,  480, 2000, 0, 300 ), WorkspaceAction::Ladder );
    EXPECT_EQ( next_workspace_action( -9, 240, 480, 2000, 0, 300 ), WorkspaceAction::Ladder );
    EXPECT_EQ( next_workspace_action( -9, 480, 480, 2000, 0, 300 ), WorkspaceAction::GiveUp );

    // -17 / -20 ( MPI buffers ): the ladder only, cap or no cap -- a cap
    // does not widen the buffers, ICNTL(14) does
    EXPECT_EQ( next_workspace_action( -20, 30,  480, 0,    2000, 300 ), WorkspaceAction::Ladder );
    EXPECT_EQ( next_workspace_action( -17, 30,  480, 0,    2000, 300 ), WorkspaceAction::Ladder );
    EXPECT_EQ( next_workspace_action( -20, 60,  480, 2000, 0,    300 ), WorkspaceAction::Ladder );
    EXPECT_EQ( next_workspace_action( -20, 480, 480, 2000, 0,    300 ), WorkspaceAction::GiveUp );

    // budget unknown ( probe returned 0 ): the ladder as before
    EXPECT_EQ( next_workspace_action( -9, 30,  480, 0, 0, 300 ), WorkspaceAction::Ladder );
    EXPECT_EQ( next_workspace_action( -9, 480, 480, 0, 0, 300 ), WorkspaceAction::GiveUp );
}

TEST( SolverParameters, CopyPreservesTuningMembers )
{
    // the copy ctor silently dropped metis nodendp once ( found by
    // review, 2026-08-15 ) — Solver takes its parameters BY VALUE, so a
    // member missing there is reset to its default before the wrapper
    // sees it. This walks the members with public setters through a
    // copy, each set to a NON-default value so a dropped member cannot
    // hide behind its default. metis nodendp has no setter, so its
    // check is default-only and would NOT catch a re-dropped member —
    // weak by necessity, noted honestly
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_compression_method( belfem::CompressionMethod::BLR );
    tParams.set_compression_cutoff( 1e-13 );
    tParams.set_absolute_tolerance( 1e-12 );
    tParams.set_relative_tolerance( 1e-6 );
    tParams.set_use_initial_guess( true );
    tParams.set_matrix_matching( false );
    tParams.set_reordering_method( belfem::ReorderingMethod::NATURAL );
    tParams.set_preconditioner( belfem::Preconditioner::GAMG );
    tParams.set_krylov_method( belfem::KrylovMethod::GMRES );
    tParams.set_memory_budget( 4000 );

    belfem::SolverParameters tCopy( tParams );
    EXPECT_EQ( tCopy.compression_method(), belfem::CompressionMethod::BLR );
    EXPECT_NEAR( tCopy.compression_cutoff(), 1e-13, tEps );
    EXPECT_TRUE( tCopy.have_compression_cutoff() );
    EXPECT_NEAR( tCopy.absolute_tolerance(), 1e-12, tEps );
    EXPECT_TRUE( tCopy.have_absolute_tolerance() );
    EXPECT_NEAR( tCopy.relative_tolerance(), 1e-6, tEps );
    EXPECT_TRUE( tCopy.have_relative_tolerance() );
    EXPECT_TRUE( tCopy.use_initial_guess() );
    EXPECT_FALSE( tCopy.use_matrix_matching() );
    EXPECT_EQ( tCopy.reordering_method(), belfem::ReorderingMethod::NATURAL );
    EXPECT_EQ( tCopy.preconditioner(), belfem::Preconditioner::GAMG );
    EXPECT_EQ( tCopy.krylov_method(), belfem::KrylovMethod::GMRES );
    EXPECT_TRUE( tCopy.use_metis_nodendp() );
    EXPECT_EQ( tCopy.memory_budget(), 4000u );
    EXPECT_TRUE( tCopy.have_memory_budget() );
}

TEST( SolverParameters, SetPreconditioner )
{
    // (idea from ChatGPT)
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_preconditioner( belfem::Preconditioner::ILU );
    EXPECT_EQ( tParams.preconditioner(), belfem::Preconditioner::ILU );
}

TEST( SolverParameters, SetKrylovMethod )
{
    // (idea from ChatGPT)
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_krylov_method( belfem::KrylovMethod::GMRES );
    EXPECT_EQ( tParams.krylov_method(), belfem::KrylovMethod::GMRES );
}

TEST( SolverParameters, SetDistributedMatrixType )
{
    // (idea from ChatGPT)
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    tParams.set_distributed_matrix_type( belfem::DistributedMatrixType::CSC );
    EXPECT_EQ( tParams.distributed_matrix_type(), belfem::DistributedMatrixType::CSC );
}

// =============================================================================
// Unlinked Backend Tests  [semantic]
// (idea from ChatGPT — BELFEM_ERROR is always active, no #ifndef NDEBUG)
// =============================================================================

#ifndef BELFEM_SUITESPARSE
TEST( SolverUnlinked, UnlinkedUMFPACKThrows )
{
    EXPECT_THROW( belfem::Solver( belfem::SolverType::UMFPACK ), std::runtime_error );
}
#endif

#ifndef BELFEM_MUMPS
TEST( SolverUnlinked, UnlinkedMUMPSThrows )
{
    EXPECT_THROW( belfem::Solver( belfem::SolverType::MUMPS ), std::runtime_error );
}
#endif

#ifndef BELFEM_PARDISO
TEST( SolverUnlinked, UnlinkedPARDISOThrows )
{
    EXPECT_THROW( belfem::Solver( belfem::SolverType::PARDISO ), std::runtime_error );
}
#endif

#ifndef BELFEM_PETSC
TEST( SolverUnlinked, UnlinkedPETScThrows )
{
    EXPECT_THROW( belfem::Solver( belfem::SolverType::PETSc ), std::runtime_error );
}
#endif

#ifndef BELFEM_STRUMPACK
TEST( SolverUnlinked, UnlinkedSTRUMPACKThrows )
{
    EXPECT_THROW( belfem::Solver( belfem::SolverType::STRUMPACK ), std::runtime_error );
}
#endif
