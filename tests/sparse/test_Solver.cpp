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
#include "arpacktools.hpp"
#include "cl_Cell.hpp"

#include <algorithm>
#include <string>

#ifdef BELFEM_MUMPS
#include "cl_SolverMUMPS.hpp"
#endif

#ifdef BELFEM_STRUMPACK
#include "strumpacktools.hpp"
#endif

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

// MKL PARDISO returns no determinant: the wrapper must say so rather than hand back a number
TEST( SolverSolve, PARDISODeterminantIsUnsupported )
{
    belfem::Matrix< belfem::real > tK = make_tridiag( 4 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tB( 4, 1.0 );
    belfem::Vector< belfem::real > tX( 4, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::PARDISO );
    tSolver.solve( tM, tX, tB );

    EXPECT_THROW( tSolver.wrapper()->get_determinant(), std::runtime_error );
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

TEST( SolverParameters, RelativeToleranceRejectsInvalid )
{
    // the setter is the single validation path ( the deck parser calls it ).
    // Mutation: drop the positivity conjunct and the first two lines go red
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );
    EXPECT_THROW( tParams.set_relative_tolerance( 0.0 ), std::runtime_error );
    EXPECT_THROW( tParams.set_relative_tolerance( -1e-6 ), std::runtime_error );
    EXPECT_THROW( tParams.set_relative_tolerance(
        std::numeric_limits< belfem::real >::infinity() ), std::runtime_error );
    EXPECT_THROW( tParams.set_relative_tolerance(
        std::numeric_limits< belfem::real >::quiet_NaN() ), std::runtime_error );

    // a rejected value must not have touched the state
    EXPECT_NEAR( tParams.relative_tolerance(), 1e-10, tEps );
    EXPECT_FALSE( tParams.have_relative_tolerance() );
}

TEST( SolverParameters, RelativeToleranceAcceptsValid )
{
    // the pair to the refusals: a guard that refused everything would pass them
    belfem::SolverParameters tParams( belfem::SolverType::UMFPACK );

    tParams.set_relative_tolerance( 1e-14 );
    EXPECT_EQ( tParams.relative_tolerance(), 1e-14 );
    EXPECT_TRUE( tParams.have_relative_tolerance() );

    tParams.set_relative_tolerance( 1e-2 );
    EXPECT_EQ( tParams.relative_tolerance(), 1e-2 );
}

TEST( SolverParameters, MaxIterationsSetter )
{
    // unset means the library default; a count must be positive, and a
    // negative one is refused rather than wrapped into a huge unsigned budget
    belfem::SolverParameters tParams( belfem::SolverType::PETSc );
    EXPECT_EQ( tParams.max_iterations(), 0u );

    EXPECT_THROW( tParams.set_max_iterations( 0 ), std::runtime_error );
    EXPECT_THROW( tParams.set_max_iterations( -1 ), std::runtime_error );
    EXPECT_EQ( tParams.max_iterations(), 0u );

    tParams.set_max_iterations( 7 );
    EXPECT_EQ( tParams.max_iterations(), 7u );

    belfem::SolverParameters tCopy( tParams );
    EXPECT_EQ( tCopy.max_iterations(), 7u );
}

#ifdef BELFEM_PETSC
// an exhausted iteration budget is a divergence, and with soft failure off
// it is a hard error naming the reason. Unpreconditioned CG is used because
// the default ASM preconditioner degenerates to ILU on one block, which is
// the exact LU of a tridiagonal and would converge in one iteration.
// Mutation: drop the negative-reason branch after KSPSolve and this goes red
TEST( PetscSolve, OneIterationCapIsAHardDivergence )
{
    const belfem::uint tN = 6;
    belfem::Matrix< belfem::real > tK = make_tridiag( tN );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tB( tN, 1.0 );
    belfem::Vector< belfem::real > tX( tN, 0.0 );

    // the solver copies its parameters at construction, so they are set first
    belfem::SolverParameters tParams( belfem::SolverType::PETSc );
    tParams.set_krylov_method( belfem::KrylovMethod::CG );
    tParams.set_preconditioner( belfem::Preconditioner::NONE );
    tParams.set_relative_tolerance( 1e-14 );
    tParams.set_max_iterations( 1 );

    belfem::Solver tSolver( tParams );

    bool tThrown = false;
    try
    {
        tSolver.solve( tM, tX, tB );
    }
    catch( const std::runtime_error & aError )
    {
        tThrown = true;
        const belfem::string tWhat( aError.what() );
        EXPECT_NE( tWhat.find( "did not converge" ), belfem::string::npos ) << tWhat;
        EXPECT_NE( tWhat.find( "DIVERGED_ITS" ), belfem::string::npos ) << tWhat;
    }
    EXPECT_TRUE( tThrown );
}
#endif

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
    tParams.set_distributed_matrix_type( belfem::DistributedMatrixType::CSC );

    belfem::SolverParameters tCopy( tParams );
    EXPECT_EQ( tCopy.type(), belfem::SolverType::UMFPACK );
    EXPECT_EQ( tCopy.distributed_matrix_type(), belfem::DistributedMatrixType::CSC );
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

// =============================================================================
// §7  Enum string round trip  [semantic]
// =============================================================================

namespace
{
    // every enumerator below UNDEFINED has its own name, and the parser maps
    // that name back to it. aAlias names the one enumerator that is allowed
    // to share a spelling with an earlier one; it is excluded from the
    // distinctness and round-trip checks and pinned by the caller instead
    template< typename E >
    void
    expect_round_trip(
            E ( *aParse )( const belfem::string & ),
            const E aAlias = E::UNDEFINED )
    {
        belfem::Cell< belfem::string > tNames;

        for( belfem::uint k = 0; k < static_cast< belfem::uint >( E::UNDEFINED ); ++k )
        {
            const E tValue = static_cast< E >( k );
            const belfem::string tName = belfem::to_string( tValue );

            EXPECT_NE( belfem::string_to_lower( tName ), "unknown" ) << "enumerator " << k;

            if( tValue == aAlias )
            {
                continue;
            }

            EXPECT_EQ( aParse( tName ), tValue ) << "enumerator " << k << " '" << tName << "'";
            tNames.push( belfem::string_to_lower( tName ) );
        }

        std::sort( tNames.begin(), tNames.end() );
        EXPECT_TRUE( std::adjacent_find( tNames.begin(), tNames.end() ) == tNames.end() )
            << "two enumerators share a name";

        EXPECT_THROW( aParse( "no-such-enumerator" ), std::runtime_error );
    }
}

TEST( SolverEnums, RoundTripsThroughStrings )
{
    expect_round_trip< belfem::SolverType >( belfem::solver_type );
    expect_round_trip< belfem::DistributedMatrixType >( belfem::distributed_matrix_type );
    expect_round_trip< belfem::KrylovMethod >( belfem::krylov_method );
    expect_round_trip< belfem::ReorderingMethod >( belfem::reordering_method );
    expect_round_trip< belfem::CompressionMethod >( belfem::compression_method );

    // GASM is PCASM with the restricted type, so it stringifies to PETSc's
    // "asm" and the parser returns ASM for that name: GASM is not selectable
    // by name, which is the documented input contract
    expect_round_trip< belfem::Preconditioner >( belfem::preconditioner,
                                                 belfem::Preconditioner::GASM );
    EXPECT_EQ( belfem::to_string( belfem::Preconditioner::GASM ), "asm" );
    EXPECT_EQ( belfem::preconditioner( belfem::to_string( belfem::Preconditioner::GASM ) ),
               belfem::Preconditioner::ASM );
}

TEST( SolverEnums, ParsingIsCaseInsensitive )
{
    EXPECT_EQ( belfem::solver_type( "MuMpS" ), belfem::SolverType::MUMPS );
    EXPECT_EQ( belfem::solver_type( "petsc" ), belfem::SolverType::PETSc );
    EXPECT_EQ( belfem::krylov_method( "GMRES" ), belfem::KrylovMethod::GMRES );
    EXPECT_EQ( belfem::preconditioner( "ILU" ), belfem::Preconditioner::ILU );
    EXPECT_EQ( belfem::reordering_method( "ParMETIS" ), belfem::ReorderingMethod::PARMETIS );
    EXPECT_EQ( belfem::compression_method( "BLR" ), belfem::CompressionMethod::BLR );
    EXPECT_EQ( belfem::distributed_matrix_type( "CSR" ), belfem::DistributedMatrixType::CSR );
}

// =============================================================================
// §8  ARPACK flag decoding  [semantic]
// =============================================================================

#ifdef BELFEM_ARPACK

// 0, 1, 3 and -9999 leave the verdict to the converged count; every other
// ARPACK flag and every driver code >= 100 aborts with the flag in the message
TEST( ArpackMessages, CheckNaupdClassifiesFlags )
{
    EXPECT_NO_THROW( belfem::arpack::check_naupd( 0 ) );
    EXPECT_NO_THROW( belfem::arpack::check_naupd( 1 ) );
    EXPECT_NO_THROW( belfem::arpack::check_naupd( 3 ) );
    EXPECT_NO_THROW( belfem::arpack::check_naupd( -9999 ) );

    const belfem::int_t tFatal[] = { -1, -3, -9, -12, 100, 101, 102, 103, 104, 105 };

    for( const belfem::int_t tFlag : tFatal )
    {
        bool tThrown = false;
        try
        {
            belfem::arpack::check_naupd( tFlag );
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true;
            const belfem::string tWhat( aError.what() );
            EXPECT_NE( tWhat.find( "error " + std::to_string( ( long ) tFlag ) + " " ),
                       belfem::string::npos ) << tWhat;
            EXPECT_EQ( tWhat.find( '%' ), belfem::string::npos ) << tWhat;
        }
        EXPECT_TRUE( tThrown ) << "flag " << tFlag;
    }
}

#endif // BELFEM_ARPACK

// =============================================================================
// §9  STRUMPACK option mapping  [semantic]
// =============================================================================

#ifdef BELFEM_STRUMPACK

namespace
{
    // an options object whose every field the wrapper owns is set away from
    // what the wrapper writes, so a missing assignment cannot hide behind a
    // library default that happens to agree
    void
    dirty_strumpack_options( belfem::sparse::StrumpackOptions & aOpts )
    {
        aOpts.set_rel_tol( 0.5 );
        aOpts.set_abs_tol( 0.5 );
        aOpts.set_maxit( 7 );
        aOpts.set_compression( strumpack::CompressionType::BLR );
        aOpts.set_Krylov_solver( strumpack::KrylovSolver::REFINE );
        aOpts.enable_METIS_NodeND();
    }
}

TEST( StrumpackOptions, DefaultsReachTheWrapper )
{
    belfem::SolverParameters tParams( belfem::SolverType::STRUMPACK );

    belfem::sparse::StrumpackOptions tOpts;
    dirty_strumpack_options( tOpts );
    belfem::sparse::set_strumpack_options( tParams, tOpts, 1000 );

    EXPECT_EQ( tOpts.rel_tol(), 1e-10 );
    EXPECT_EQ( tOpts.abs_tol(), 1e-14 );
    EXPECT_EQ( tOpts.maxit(), 50 );
    EXPECT_EQ( tOpts.compression(), strumpack::CompressionType::NONE );
    EXPECT_EQ( tOpts.Krylov_solver(), strumpack::KrylovSolver::PREC_GMRES );
    EXPECT_NE( tOpts.matching(), strumpack::MatchingJob::NONE );
    EXPECT_TRUE( tOpts.use_METIS_NodeNDP() );

    // OFF compresses no more than the AUTOMATIC default does
    tParams.set_compression_method( belfem::CompressionMethod::OFF );
    belfem::sparse::StrumpackOptions tOff;
    dirty_strumpack_options( tOff );
    belfem::sparse::set_strumpack_options( tParams, tOff, 1000 );
    EXPECT_EQ( tOff.compression(), strumpack::CompressionType::NONE );

    // the stated tolerances win over the class defaults
    tParams.set_relative_tolerance( 1e-7 );
    tParams.set_absolute_tolerance( 1e-12 );
    belfem::sparse::StrumpackOptions tStated;
    dirty_strumpack_options( tStated );
    belfem::sparse::set_strumpack_options( tParams, tStated, 1000 );
    EXPECT_EQ( tStated.rel_tol(), 1e-7 );
    EXPECT_EQ( tStated.abs_tol(), 1e-12 );
}

TEST( StrumpackOptions, BlrCutoffAndLeafSize )
{
    belfem::SolverParameters tParams( belfem::SolverType::STRUMPACK );
    tParams.set_compression_method( belfem::CompressionMethod::BLR );
    tParams.set_compression_cutoff( 1e-13 );

    // leaf size 128 below 100000 rows, 256 below 250000, 512 from there on
    const belfem::index_t tRows[]  = { 50000, 99999, 100000, 150000, 249999, 250000, 400000 };
    const int             tLeafs[] = { 128,   128,   256,    256,    256,    512,    512 };

    for( belfem::uint k = 0; k < 7; ++k )
    {
        belfem::sparse::StrumpackOptions tOpts;
        dirty_strumpack_options( tOpts );
        belfem::sparse::set_strumpack_options( tParams, tOpts, tRows[ k ] );

        EXPECT_EQ( tOpts.compression(), strumpack::CompressionType::BLR );
        EXPECT_EQ( tOpts.BLR_options().rel_tol(), 1e-13 );
        EXPECT_EQ( tOpts.BLR_options().leaf_size(), tLeafs[ k ] ) << tRows[ k ] << " rows";
    }
}

TEST( StrumpackOptions, KrylovMapping )
{
    // CG, CGS, IBCGS and TFQMR have no STRUMPACK counterpart: refused, not
    // silently replaced by the library's own AUTO
    const belfem::KrylovMethod tUnmapped[] = {
        belfem::KrylovMethod::CG, belfem::KrylovMethod::CGS,
        belfem::KrylovMethod::IBCGS, belfem::KrylovMethod::TFQMR };

    for( const belfem::KrylovMethod tMethod : tUnmapped )
    {
        belfem::SolverParameters tParams( belfem::SolverType::STRUMPACK );
        tParams.set_krylov_method( tMethod );
        belfem::sparse::StrumpackOptions tOpts;
        EXPECT_THROW( belfem::sparse::set_strumpack_options( tParams, tOpts, 1000 ),
                      std::runtime_error ) << belfem::to_string( tMethod );
    }

    struct Case
    {
        belfem::KrylovMethod   mMethod;
        belfem::Preconditioner mPreconditioner;
        strumpack::KrylovSolver mExpected;
    };

    const Case tCases[] = {
        { belfem::KrylovMethod::AUTO,    belfem::Preconditioner::ASM,  strumpack::KrylovSolver::PREC_GMRES },
        { belfem::KrylovMethod::AUTO,    belfem::Preconditioner::NONE, strumpack::KrylovSolver::PREC_GMRES },
        { belfem::KrylovMethod::PREONLY, belfem::Preconditioner::ASM,  strumpack::KrylovSolver::DIRECT },
        { belfem::KrylovMethod::GMRES,   belfem::Preconditioner::ASM,  strumpack::KrylovSolver::PREC_GMRES },
        { belfem::KrylovMethod::GMRES,   belfem::Preconditioner::NONE, strumpack::KrylovSolver::GMRES },
        { belfem::KrylovMethod::BCGS,    belfem::Preconditioner::ASM,  strumpack::KrylovSolver::PREC_BICGSTAB },
        { belfem::KrylovMethod::BCGS,    belfem::Preconditioner::NONE, strumpack::KrylovSolver::BICGSTAB } };

    for( const Case & tCase : tCases )
    {
        belfem::SolverParameters tParams( belfem::SolverType::STRUMPACK );
        tParams.set_krylov_method( tCase.mMethod );
        tParams.set_preconditioner( tCase.mPreconditioner );

        belfem::sparse::StrumpackOptions tOpts;
        dirty_strumpack_options( tOpts );
        EXPECT_NO_THROW( belfem::sparse::set_strumpack_options( tParams, tOpts, 1000 ) );
        EXPECT_EQ( tOpts.Krylov_solver(), tCase.mExpected )
            << belfem::to_string( tCase.mMethod ) << " / " << belfem::to_string( tCase.mPreconditioner );
    }
}

// NodeNDP and MC64 matching are both on by default; matching has a setter
// for its opt-out. The NodeNDP opt-out is reachable only from a deck
TEST( StrumpackOptions, MetisNodeNdpAndMatchingDefaults )
{
    belfem::SolverParameters tParams( belfem::SolverType::STRUMPACK );

    belfem::sparse::StrumpackOptions tOpts;
    dirty_strumpack_options( tOpts );
    belfem::sparse::set_strumpack_options( tParams, tOpts, 1000 );
    EXPECT_TRUE( tOpts.use_METIS_NodeNDP() );
    EXPECT_NE( tOpts.matching(), strumpack::MatchingJob::NONE );

    tParams.set_matrix_matching( false );
    belfem::sparse::StrumpackOptions tOff;
    belfem::sparse::set_strumpack_options( tParams, tOff, 1000 );
    EXPECT_EQ( tOff.matching(), strumpack::MatchingJob::NONE );
}

#endif // BELFEM_STRUMPACK

// =============================================================================
// §10  MUMPS warnings, lifecycle and diagnostics  [semantic]
// =============================================================================

#ifdef BELFEM_MUMPS

namespace
{
    // warning_message() is protected; the decoder needs no MUMPS instance,
    // only the object, so a derived class exposes it
    class MumpsWarningProbe : public belfem::solver::MUMPS
    {
    public:
        explicit MumpsWarningProbe( const belfem::SolverParameters * aParams ) :
            belfem::solver::MUMPS( aParams )
        {
        }

        using belfem::solver::MUMPS::warning_message;
    };

    void
    decode_mumps_warning(
            const belfem::int_t aInfo1,
            const belfem::int_t aInfo2,
            belfem::Cell< belfem::string > & aWarnings )
    {
        belfem::SolverParameters tParams( belfem::SolverType::MUMPS );
        MumpsWarningProbe tProbe( &tParams );
        const belfem::int_t tInfo[ 2 ] = { aInfo1, aInfo2 };
        aWarnings.clear();
        tProbe.warning_message( tInfo, aWarnings );
    }

    bool
    contains( const belfem::string & aText, const belfem::string & aPattern )
    {
        return aText.find( aPattern ) != belfem::string::npos;
    }

    belfem::int_t
    mumps_pool_occupancy()
    {
        belfem::int_t tCount = -1;
        mumpstools_num_solvers( tCount );
        return tCount;
    }
}

TEST( MumpsWarnings, SingleBitsDecode )
{
    belfem::Cell< belfem::string > tWarnings;

    for( const belfem::int_t tBit : { 1, 2, 4, 8, 16 } )
    {
        decode_mumps_warning( tBit, 4242, tWarnings );
        ASSERT_EQ( tWarnings.size(), 1u ) << "INFO(1) = " << tBit;
        EXPECT_FALSE( contains( tWarnings( 0 ), "%" ) ) << tWarnings( 0 );
    }

    // bits 1 and 16 alone own INFO(2) and report it
    decode_mumps_warning( 1, 4242, tWarnings );
    EXPECT_TRUE( contains( tWarnings( 0 ), "4242" ) ) << tWarnings( 0 );
    decode_mumps_warning( 16, 4242, tWarnings );
    EXPECT_TRUE( contains( tWarnings( 0 ), "4242" ) ) << tWarnings( 0 );

    // the other bits carry no detail, so INFO(2) must not leak into them
    for( const belfem::int_t tBit : { 2, 4, 8 } )
    {
        decode_mumps_warning( tBit, 4242, tWarnings );
        EXPECT_FALSE( contains( tWarnings( 0 ), "4242" ) ) << tWarnings( 0 );
    }
}

TEST( MumpsWarnings, CombinedBitsDecodeIndependently )
{
    belfem::Cell< belfem::string > tWarnings;

    decode_mumps_warning( 6, 4242, tWarnings );
    EXPECT_EQ( tWarnings.size(), 2u );

    decode_mumps_warning( 31, 4242, tWarnings );
    EXPECT_EQ( tWarnings.size(), 5u );

    // with bits 1 and 16 both raised, INFO(2) belongs to whichever MUMPS set
    // last and cannot be attributed: neither line may carry it
    decode_mumps_warning( 17, 4242, tWarnings );
    ASSERT_EQ( tWarnings.size(), 2u );
    EXPECT_FALSE( contains( tWarnings( 0 ), "4242" ) ) << tWarnings( 0 );
    EXPECT_FALSE( contains( tWarnings( 1 ), "4242" ) ) << tWarnings( 1 );
}

TEST( MumpsWarnings, UnknownBitStillReports )
{
    belfem::Cell< belfem::string > tWarnings;

    decode_mumps_warning( 64, 4242, tWarnings );
    ASSERT_EQ( tWarnings.size(), 1u );
    EXPECT_TRUE( contains( tWarnings( 0 ), "64" ) ) << tWarnings( 0 );
    EXPECT_FALSE( contains( tWarnings( 0 ), "%" ) ) << tWarnings( 0 );

    // an unknown bit next to a known one is reported too, not swallowed by it.
    // Mutation: report unknown bits only when no known bit matched, and this
    // goes red
    decode_mumps_warning( 65, 4242, tWarnings );
    ASSERT_EQ( tWarnings.size(), 2u );
    EXPECT_TRUE( contains( tWarnings( 1 ), "bits 64" ) ) << tWarnings( 1 );
}

// a second initialize() on a live wrapper is refused BEFORE a solver slot is
// reserved, so the refusal leaves the pool occupancy unchanged
TEST( SolverLifecycle, MumpsRefusesReInitialization )
{
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

    const belfem::int_t tBase = mumps_pool_occupancy();

    belfem::Solver tSolver( belfem::SolverType::MUMPS );
    tSolver.wrapper()->initialize( tM, belfem::SymmetryMode::Unsymmetric, 1 );
    ASSERT_EQ( mumps_pool_occupancy(), tBase + 1 );

    EXPECT_THROW( tSolver.wrapper()->initialize( tM, belfem::SymmetryMode::Unsymmetric, 1 ),
                  std::runtime_error );
    EXPECT_EQ( mumps_pool_occupancy(), tBase + 1 );
    EXPECT_TRUE( tSolver.wrapper()->is_initialized() );

    tSolver.free();
    EXPECT_EQ( mumps_pool_occupancy(), tBase );
}

// the guard keys on the initialized state, not on having been used once
TEST( SolverLifecycle, MumpsReInitializesAfterFree )
{
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

    const belfem::int_t tBase = mumps_pool_occupancy();

    belfem::Solver tSolver( belfem::SolverType::MUMPS );
    tSolver.wrapper()->initialize( tM, belfem::SymmetryMode::Unsymmetric, 1 );
    tSolver.free();
    ASSERT_EQ( mumps_pool_occupancy(), tBase );

    EXPECT_NO_THROW( tSolver.wrapper()->initialize( tM, belfem::SymmetryMode::Unsymmetric, 1 ) );
    EXPECT_EQ( mumps_pool_occupancy(), tBase + 1 );

    tSolver.free();
    EXPECT_EQ( mumps_pool_occupancy(), tBase );
}

// an error analysis armed before the first solve survives the lazy
// initialize() inside it, so the capture after the solve is a number
TEST( MumpsDiagnostics, ArmSurvivesInitialize )
{
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tB( 6, 1.0 );
    belfem::Vector< belfem::real > tX( 6, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::MUMPS );
    tSolver.set_mumps_error_analysis( belfem::MumpsErrorAnalysis::Full );
    tSolver.solve( tM, tX, tB );

    const belfem::real tCond1 = tSolver.wrapper()->get_cond1();
    EXPECT_TRUE( std::isfinite( tCond1 ) ) << tCond1;
    EXPECT_GE( tCond1, 1.0 );
}

// unarmed, the capture reads NaN rather than whatever RINFOG holds
TEST( MumpsDiagnostics, UnarmedCaptureIsNaN )
{
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tB( 6, 1.0 );
    belfem::Vector< belfem::real > tX( 6, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::MUMPS );
    tSolver.set_mumps_error_analysis( belfem::MumpsErrorAnalysis::None );
    tSolver.solve( tM, tX, tB );

    EXPECT_TRUE( std::isnan( tSolver.wrapper()->get_cond1() ) );
    EXPECT_TRUE( std::isnan( tSolver.wrapper()->get_cond2() ) );
}

#endif // BELFEM_MUMPS

#ifdef BELFEM_STRUMPACK

// the controller arms the MUMPS error analysis unconditionally, so on any
// other solver type the call must do nothing
TEST( SolverDiagnostics, ErrorAnalysisIsANoOpOnNonMumps )
{
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

    belfem::Vector< belfem::real > tXknown( 6 );
    for( belfem::uint k = 0; k < 6; ++k )
    {
        tXknown( k ) = static_cast< belfem::real >( k + 1 );
    }
    belfem::Vector< belfem::real > tB( tK * tXknown );
    belfem::Vector< belfem::real > tX( 6, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::STRUMPACK );
    EXPECT_NO_THROW( tSolver.set_mumps_error_analysis( belfem::MumpsErrorAnalysis::Full ) );
    tSolver.solve( tM, tX, tB );

    for( belfem::uint k = 0; k < 6; ++k )
    {
        EXPECT_NEAR( tX( k ), tXknown( k ), tTol );
    }
}

#endif // BELFEM_STRUMPACK

// =============================================================================
// §11  SuperLU  [semantic]
// =============================================================================

namespace
{
    // nonsymmetric and diagonally dominant: 4 on the diagonal, -1 below,
    // -2 above, and one far corner entry. A symmetric matrix would hide a
    // transpose: A and A^T give the same answer only when A = A^T
    [[maybe_unused]] belfem::Matrix< belfem::real >
    make_asymmetric( const belfem::uint aN, const belfem::real aDiagonal = 4.0 )
    {
        belfem::Matrix< belfem::real > tK( aN, aN, 0.0 );
        for( belfem::uint k = 0; k < aN; ++k )
        {
            tK( k, k ) = aDiagonal;
            if( k > 0 )      tK( k, k - 1 ) = -1.0;
            if( k < aN - 1 ) tK( k, k + 1 ) = -2.0;
        }
        tK( aN - 1, 0 ) = 0.5;
        return tK;
    }

    [[maybe_unused]] belfem::Vector< belfem::real >
    make_known_solution( const belfem::uint aN )
    {
        belfem::Vector< belfem::real > tX( aN );
        for( belfem::uint k = 0; k < aN; ++k )
        {
            tX( k ) = static_cast< belfem::real >( k + 1 );
        }
        return tX;
    }

    [[maybe_unused]] void
    expect_solution(
            const belfem::Vector< belfem::real > & aX,
            const belfem::Vector< belfem::real > & aXknown )
    {
        ASSERT_EQ( aX.length(), aXknown.length() );
        for( belfem::uint k = 0; k < aXknown.length(); ++k )
        {
            EXPECT_NEAR( aX( k ), aXknown( k ), tTol ) << "row " << k;
        }
    }
}

#ifdef BELFEM_SUPERLU

// CSR arrays are the CSC arrays of A^T, so the wrapper hands both formats to
// the same SLU_NC store and undoes the transpose with Trans = TRANS for CSR
TEST( SuperLUSolve, CsrAndCscAgree )
{
    const belfem::uint tN = 6;
    belfem::Matrix< belfem::real > tK = make_asymmetric( tN );
    belfem::Vector< belfem::real > tXknown = make_known_solution( tN );
    belfem::Vector< belfem::real > tB( tK * tXknown );

    for( const belfem::SpMatrixType tType : { belfem::SpMatrixType::CSR, belfem::SpMatrixType::CSC } )
    {
        belfem::SpMatrix tM( tK, tType );
        belfem::Vector< belfem::real > tX( tN, 0.0 );

        belfem::Solver tSolver( belfem::SolverType::SUPERLU );
        tSolver.solve( tM, tX, tB );

        SCOPED_TRACE( tType == belfem::SpMatrixType::CSR ? "CSR" : "CSC" );
        expect_solution( tX, tXknown );
    }
}

// a value change on the same matrix object is picked up by the next solve
TEST( SuperLULifecycle, ReSolveFollowsNewValues )
{
    const belfem::uint tN = 6;
    belfem::Matrix< belfem::real > tK = make_asymmetric( tN );
    belfem::Vector< belfem::real > tXknown = make_known_solution( tN );
    belfem::Vector< belfem::real > tB( tK * tXknown );

    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tX( tN, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::SUPERLU );
    tSolver.solve( tM, tX, tB );
    expect_solution( tX, tXknown );

    // same pattern, new diagonal
    belfem::Matrix< belfem::real > tK2 = make_asymmetric( tN, 7.0 );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    for( belfem::uint k = 0; k < tN; ++k )
    {
        tM( k, k ) = 7.0;
    }

    belfem::Vector< belfem::real > tB2( tK2 * tXknown );
    belfem::Vector< belfem::real > tX2( tN, 0.0 );
    tSolver.solve( tM, tX2, tB2 );
    expect_solution( tX2, tXknown );
}

// a different matrix object on a live wrapper replaces the factorization,
// with or without an intervening free()
TEST( SuperLULifecycle, MatrixReplacementRefactorizes )
{
    const belfem::uint tN1 = 6;
    const belfem::uint tN2 = 9;

    belfem::Matrix< belfem::real > tK1 = make_asymmetric( tN1 );
    belfem::Matrix< belfem::real > tK2 = make_asymmetric( tN2, 5.0 );

    belfem::Vector< belfem::real > tX1known = make_known_solution( tN1 );
    belfem::Vector< belfem::real > tX2known = make_known_solution( tN2 );
    belfem::Vector< belfem::real > tB1( tK1 * tX1known );
    belfem::Vector< belfem::real > tB2( tK2 * tX2known );

    belfem::SpMatrix tM1( tK1, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM2( tK2, belfem::SpMatrixType::CSR );

    belfem::Solver tSolver( belfem::SolverType::SUPERLU );

    belfem::Vector< belfem::real > tX( tN1, 0.0 );
    tSolver.solve( tM1, tX, tB1 );
    expect_solution( tX, tX1known );

    // no free(): the symbolic re-entry path
    belfem::Vector< belfem::real > tY( tN2, 0.0 );
    tSolver.solve( tM2, tY, tB2 );
    expect_solution( tY, tX2known );

    tSolver.free();

    belfem::Vector< belfem::real > tZ( tN1, 0.0 );
    tSolver.solve( tM1, tZ, tB1 );
    expect_solution( tZ, tX1known );
}

TEST( SuperLULifecycle, FreeAndReuse )
{
    const belfem::uint tN = 6;
    belfem::Matrix< belfem::real > tK = make_asymmetric( tN );
    belfem::Vector< belfem::real > tXknown = make_known_solution( tN );
    belfem::Vector< belfem::real > tB( tK * tXknown );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );

    belfem::Solver tSolver( belfem::SolverType::SUPERLU );
    EXPECT_FALSE( tSolver.wrapper()->is_initialized() );

    belfem::Vector< belfem::real > tX( tN, 0.0 );
    tSolver.solve( tM, tX, tB );
    EXPECT_TRUE( tSolver.wrapper()->is_initialized() );

    tSolver.free();
    EXPECT_FALSE( tSolver.wrapper()->is_initialized() );

    belfem::Vector< belfem::real > tX2( tN, 0.0 );
    tSolver.solve( tM, tX2, tB );
    EXPECT_TRUE( tSolver.wrapper()->is_initialized() );
    expect_solution( tX2, tXknown );
}

// each column of a multi-column solve is its own single-RHS solve
TEST( SuperLUSolve, MultiRhsScatter )
{
    const belfem::uint tN = 6;
    const belfem::uint tM = 3;

    belfem::Matrix< belfem::real > tK = make_asymmetric( tN );
    belfem::SpMatrix tA( tK, belfem::SpMatrixType::CSR );

    belfem::Matrix< belfem::real > tXknown( tN, tM );
    for( belfem::uint j = 0; j < tM; ++j )
    {
        for( belfem::uint i = 0; i < tN; ++i )
        {
            tXknown( i, j ) = static_cast< belfem::real >( ( i + 1 ) * ( j + 1 ) ) - 2.0 * j;
        }
    }
    belfem::Matrix< belfem::real > tB( tK * tXknown );
    belfem::Matrix< belfem::real > tX( tN, tM, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::SUPERLU );
    tSolver.solve( tA, tX, tB );

    ASSERT_EQ( tX.n_rows(), tN );
    ASSERT_EQ( tX.n_cols(), tM );
    for( belfem::uint j = 0; j < tM; ++j )
    {
        for( belfem::uint i = 0; i < tN; ++i )
        {
            EXPECT_NEAR( tX( i, j ), tXknown( i, j ), tTol ) << "( " << i << ", " << j << " )";
        }
    }
}

// the two columns of [ 1 1 ; 1 1 ] are equal, so whatever the pivoting, the
// second pivot is exactly zero and the message names it
TEST( SuperLUFailure, SingularMatrixNamesThePivot )
{
    belfem::Matrix< belfem::real > tK( 2, 2, 1.0 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSC );
    belfem::Vector< belfem::real > tB( 2, 1.0 );
    belfem::Vector< belfem::real > tX( 2, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::SUPERLU );

    bool tThrown = false;
    try
    {
        tSolver.solve( tM, tX, tB );
    }
    catch( const std::runtime_error & aError )
    {
        tThrown = true;
        const belfem::string tWhat( aError.what() );
        EXPECT_NE( tWhat.find( "singular" ), belfem::string::npos ) << tWhat;
        EXPECT_NE( tWhat.find( "U( 2, 2 )" ), belfem::string::npos ) << tWhat;
    }
    EXPECT_TRUE( tThrown );
}

TEST( SuperLUSolve, RejectsNonSquare )
{
    belfem::Matrix< belfem::real > tK( 2, 3, 0.0 );
    tK( 0, 0 ) = 1.0;
    tK( 1, 1 ) = 1.0;
    tK( 0, 2 ) = 1.0;
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    belfem::Vector< belfem::real > tB( 2, 1.0 );
    belfem::Vector< belfem::real > tX( 2, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::SUPERLU );
    EXPECT_THROW( tSolver.solve( tM, tX, tB ), std::runtime_error );
}

#endif // BELFEM_SUPERLU

#ifdef BELFEM_ARPACK

// the n x n 1-D Laplacian has the closed-form spectrum
// 2 - 2 cos( k pi / ( n + 1 ) ), k = 1 ... n
TEST( ArpackSubspace, FiveEigenvaluesOnAKnownSpectrum )
{
    const belfem::uint tN = 30;
    belfem::Matrix< belfem::real > tK = make_tridiag( tN );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );

    const belfem::int_t tNumRows  = tN;
    const belfem::int_t tNNZ      = tM.number_of_nonzeros();
    const belfem::int_t tJob      = 1;
    const belfem::int_t tNev      = 5;
    const belfem::int_t tNcvMin   = 0;
    const belfem::real  tTolerance = 1e-12;
    const belfem::int_t tMaxIter  = 1000;
    const belfem::real  tSigma    = 0.0;

    belfem::Vector< belfem::real > tLambdaReal( tNev + 1, std::numeric_limits< belfem::real >::quiet_NaN() );
    belfem::Vector< belfem::real > tLambdaImag( tNev + 1, std::numeric_limits< belfem::real >::quiet_NaN() );
    belfem::Vector< belfem::int_t > tInfo( belfem::arpack::gNumInfoEntries, 0 );

    belfem::arpack::arpack_symmetric_eigen(
            &tNumRows, &tNNZ, tM.data(), tM.cols(), tM.pointers(),
            &tJob, &tNev, &tNcvMin, &tTolerance, &tMaxIter, &tSigma,
            tLambdaReal.data(), tLambdaImag.data(), tInfo.data() );

    EXPECT_NO_THROW( belfem::arpack::check_saupd( tInfo( belfem::arpack::gInfoNaupd ) ) );
    EXPECT_NO_THROW( belfem::arpack::check_seupd( tInfo( belfem::arpack::gInfoNeupd ) ) );
    EXPECT_GE( tInfo( belfem::arpack::gInfoNumConverged ), tNev );

    belfem::Cell< belfem::real > tComputed;
    for( belfem::int_t k = 0; k < tNev; ++k )
    {
        tComputed.push( tLambdaReal( k ) );
        EXPECT_EQ( tLambdaImag( k ), 0.0 );
    }
    std::sort( tComputed.begin(), tComputed.end() );

    // the five largest, ascending: k = n - 4 ... n
    for( belfem::int_t k = 0; k < tNev; ++k )
    {
        const belfem::real tExact = 2.0 - 2.0 * std::cos(
                static_cast< belfem::real >( tN - tNev + 1 + k ) * std::acos( -1.0 )
                / static_cast< belfem::real >( tN + 1 ) );
        EXPECT_NEAR( tComputed( k ), tExact, 1e-10 ) << "eigenvalue " << k;
    }
}

#endif // BELFEM_ARPACK

// =============================================================================
// §12  PARDISO restores the caller's indexing base  [semantic]
// =============================================================================

// SpMatrix::operator() reads through the base, so a base left at Fortran by
// the solver silently shifts every later element access
TEST( PardisoBase, SolveLeavesTheMatrixInTheCallersBase )
{
#ifdef BELFEM_PARDISO
    for( const belfem::SpMatrixIndexingBase tBase :
         { belfem::SpMatrixIndexingBase::Cpp, belfem::SpMatrixIndexingBase::Fortran } )
    {
        belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
        belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
        tM.set_indexing_base( tBase );

        belfem::Vector< belfem::real > tXknown = make_known_solution( 6 );
        belfem::Vector< belfem::real > tB( tK * tXknown );
        belfem::Vector< belfem::real > tX( 6, 0.0 );

        belfem::Solver tSolver( belfem::SolverType::PARDISO );
        tSolver.solve( tM, tX, tB );

        SCOPED_TRACE( tBase == belfem::SpMatrixIndexingBase::Cpp ? "Cpp" : "Fortran" );
        EXPECT_EQ( tM.indexing_base(), static_cast< belfem::int_t >( tBase ) );
        expect_solution( tX, tXknown );

        // a repeated solve takes the already-initialized path
        tSolver.solve( tM, tX, tB );
        EXPECT_EQ( tM.indexing_base(), static_cast< belfem::int_t >( tBase ) );

        if( tBase == belfem::SpMatrixIndexingBase::Cpp )
        {
            const belfem::SpMatrix & tConst = tM;
            EXPECT_EQ( tConst( 0, 0 ), 2.0 );
            EXPECT_EQ( tConst( 1, 0 ), -1.0 );
        }
    }
#else
    GTEST_SKIP() << "built without MKL PARDISO";
#endif
}

TEST( PardisoBase, MultiRhsSolveRestoresToo )
{
#ifdef BELFEM_PARDISO
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Matrix< belfem::real > tXknown( 6, 2 );
    for( belfem::uint j = 0; j < 2; ++j )
    {
        for( belfem::uint i = 0; i < 6; ++i )
        {
            tXknown( i, j ) = static_cast< belfem::real >( i + 1 + 10 * j );
        }
    }
    belfem::Matrix< belfem::real > tB( tK * tXknown );
    belfem::Matrix< belfem::real > tX( 6, 2, 0.0 );

    belfem::Solver tSolver( belfem::SolverType::PARDISO );
    tSolver.solve( tM, tX, tB );

    EXPECT_EQ( tM.indexing_base(), 0 );
    for( belfem::uint j = 0; j < 2; ++j )
    {
        for( belfem::uint i = 0; i < 6; ++i )
        {
            EXPECT_NEAR( tX( i, j ), tXknown( i, j ), tTol ) << "( " << i << ", " << j << " )";
        }
    }
#else
    GTEST_SKIP() << "built without MKL PARDISO";
#endif
}

TEST( PardisoBase, InitializeAndFreeRestore )
{
#ifdef BELFEM_PARDISO
    belfem::Matrix< belfem::real > tK = make_tridiag( 6 );
    belfem::SpMatrix tM( tK, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Solver tSolver( belfem::SolverType::PARDISO );
    tSolver.wrapper()->initialize( tM, belfem::SymmetryMode::Unsymmetric, 1 );
    EXPECT_EQ( tM.indexing_base(), 0 );

    tSolver.free();
    EXPECT_EQ( tM.indexing_base(), 0 );
    EXPECT_FALSE( tSolver.wrapper()->is_initialized() );
#else
    GTEST_SKIP() << "built without MKL PARDISO";
#endif
}
