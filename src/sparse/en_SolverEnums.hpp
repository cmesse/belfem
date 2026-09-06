/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_EN_SOLVERENUMS_HPP
#define BELFEM_EN_SOLVERENUMS_HPP

#include "typedefs.hpp"
#include "fn_to_enum.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    enum class SolverType
    {
        UMFPACK,
        SUPERLU,
        MUMPS,
        STRUMPACK,
        PARDISO,
        PETSc,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class SymmetryMode
    {
        Unsymmetric               = 0,
        PositiveDefiniteSymmetric = 1,
        GeneralSymmetric          = 2,
        UNDEFINED
    };

//------------------------------------------------------------------------------


    enum class EulerMethod
    {
        Static,            //
        ForwardExplicit,   // theta = 0.0
        CrankNicolson,     // theta = 0.5
        Galerkin,          // theta = 2/3
        BackwardDifference1,              // theta = 1.0
        BackwardDifference2,
        BackwardDifference3,
        BackwardDifference4,
        BackwardDifference5,
        StiffnessOnly, //
        MassOnly, //
        Derivative,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class DistributedMatrixType
    {
        AIJ = 0,
        CSR = 1,
        CSC = 2,
        UNDEFINED = 3
    };

//------------------------------------------------------------------------------

    /**
     * PETSC only
     */
    enum class Preconditioner
    {
        NONE,    // No preconditioner
        ASM,     // Additive Schwarz
        GASM,    // Restricted additive Schwarz ( PCASM with PC_ASM_RESTRICT, NOT PETSc's PCGASM )
        GAMG,    // Geometric algebraic multigrid
        JACOBI,  // i.e. diagonal scaling preconditioning)
        BJACOBI, // Block Jacobi
        LU,      // direct solver, based on LU factorization, as a preconditioner
        ICC,     // Incomplete Cholesky factorization
        ILU,     // Incomplete factorization
        HMG,     // Hybrid of PETSc preconditioners
        SPAI,    // Use the Sparse Approximate Inverse method of Grote and Barnard as a preconditioner
        UNDEFINED
    };


//------------------------------------------------------------------------------

    /**
     * PETSC only
     */
    enum class KrylovMethod
    {
        PREONLY, // No Krylov Method
        CG,      // Conjugate Gradient
        CGS,     // Conjugate Gradient Squared
        IBCGS,   // Improved Stabilized version of BiConjugate Gradient
        GMRES,   // Generalized Minimal Residual
        TFQMR,   // transpose free QMR
        BCGS,    // Stabilized version of BiConjugate Gradient
        AUTO,    // PETSc: PREONLY if the preconditioner is LU, else GMRES; STRUMPACK: always factorization-preconditioned GMRES
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class ReorderingMethod
    {
        NATURAL   = 0,
        METIS     = 1,
        SCOTCH    = 2,
        AUTOMATIC = 3,
        PARMETIS  = 4,   // explicit parallel ND: PETSc path uses BELFEM's collective
        PTSCOTCH  = 5,   // wrapper; MUMPS / STRUMPACK map it as metis / scotch do in parallel
        UNDEFINED = 6    // appended, never inserted: the value travels as a uint
    };

//------------------------------------------------------------------------------

    enum class CompressionMethod
    {
        OFF       = 0,
        BLR       = 1,
        AUTOMATIC = 2,
        UNDEFINED = 3
    };

//------------------------------------------------------------------------------

    /**
     * MUMPS only
     */
    enum class MumpsSerialReodrdering
    {
        AMD       = 0,   //Approximate Minimum Degree
        USERPIVOT = 1,   // User Pivot ( currently not supported )
        AMF       = 2,   // Approximate Minimum Fill
        SCOTCH    = 3,   // SCOTCH
        PORD      = 4,   // PORD ordering ( shipped with MUMPS )
        METIS     = 5,   // METIS nested dissection
        QAMD      = 6,   // Approximate Minimum Degree with automatic quasi-dense row detection
        AUTOMATIC = 7,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    /**
     * MUMPS only
     */
    enum class MumpsParallelReodrdering
    {
        AUTOMATIC = 0,
        PTSCOTCH  = 1,
        PARMETIS  = 2,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    /**
     * MUMPS only
     */
     enum class MumpsBlockLowRanking
     {
         Off                      = 0,
         Automatic                = 1,
         FactorizationAndSolution = 2,
         FactorizationOnly        = 3,
         UNDEFINED                = 4
     };

//------------------------------------------------------------------------------

    /**
     * MUMPS only
     */
    enum class MumpsErrorAnalysis
    {
        None                     = 0,
        Full                     = 1,
        Main                     = 2,
        UNDEFINED                = 3
    };


//------------------------------------------------------------------------------

// set the default solver, PETSc is deliberately not set by default
#ifdef BELFEM_STRUMPACK
    const SolverType gDefaultSolver = SolverType::STRUMPACK;
#elif BELFEM_MUMPS
    const SolverType gDefaultSolver = SolverType::MUMPS;
#elif  BELFEM_PARDISO
    const SolverType gDefaultSolver = SolverType::PARDISO;
#elif BELFEM_SUITESPARSE
    const SolverType gDefaultSolver = SolverType::UMFPACK;
#elif BELFEM_SUPERLU
    const SolverType gDefaultSolver = SolverType::SUPERLU;
#else
    const SolverType gDefaultSolver = SolverType::UNDEFINED;
#endif


//------------------------------------------------------------------------------
// conversion tools
//------------------------------------------------------------------------------

    string
    to_string( const SolverType aSolverType );

    SolverType
    solver_type( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const DistributedMatrixType aDistributedMatrixType );

    DistributedMatrixType
    distributed_matrix_type( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const EulerMethod aEulerMethod );

    EulerMethod
    euler_method( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const Preconditioner aPreconditioner ) ;

    Preconditioner
    preconditioner( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const KrylovMethod aKrylovMethod ) ;

    KrylovMethod
    krylov_method( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const ReorderingMethod aReorderingMethod ) ;

    ReorderingMethod
    reordering_method( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const CompressionMethod aCompressionMethod ) ;

    CompressionMethod
    compression_method( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const MumpsSerialReodrdering aMumpsSerialReodrdering ) ;

    MumpsSerialReodrdering
    serial_reordering( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const MumpsParallelReodrdering aMumpsParallelReodrdering ) ;

    MumpsParallelReodrdering
    parallel_reordering( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const MumpsBlockLowRanking aBLR );

    MumpsBlockLowRanking
    block_low_ranking( const string & aString );

//------------------------------------------------------------------------------
}

#endif //BELFEM_EN_SOLVERENUMS_HPP
