//
// Created by Christian Messe on 10.07.20.
//

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
        MUMPS,
        STRUMPACK,
        PARDISO,
        PETSC,
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

    /**
     * For TATCAD IWG
     */
    enum class EulerMethod
    {
        ForwardExplicit,   // theta = 0.0
        CrankNicolson,     // theta = 0.5
        Galerkin,          // theta = 2/3
        BackwardImplicit,  // theta = 1.0
        UNDEFINED
    };

//------------------------------------------------------------------------------

    /**
     * PETSC only
     */
    enum class Preconditioner
    {
        NONE,    // No preconditioner
        ASM,     // Additive Schwarz
        GASM,    // Restricted additive Schwarz
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
        AUTO,    // for PETSc: GMRES, for STRUMPACK: GMRES if iterative, otherwise none
        UNDEFINED
    };

//------------------------------------------------------------------------------

    /**
     * MUMPS only
     */
    enum class SerialReodrdering
    {
        AMD       = 0,   //Approximate Minimum Degree
        USERPIVOT = 1,   // User Pivot ( currently not supported )
        AMF       = 2,   // Approximate Minimum Fill
        SCOTCH    = 3,   // SCOTCH
        PORD      = 4,   // Generalized Minimal Residual
        METIS     = 5,   // transpose free QMR
        QAMD      = 6,   // Approximate Minimum Degree with automatic quasi-dense row detection
        AUTOMATIC = 7,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    /**
     * MUMPS only
     */
    enum class ParallelReodrdering
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
     enum class BlockLowRanking
     {
         Off                      = 0,
         Automatic                = 1,
         FactorizationAndSolution = 2,
         FactorizationOnly        = 3,
         UNDEFINED                = 4
     };

//------------------------------------------------------------------------------

// set the default solver, PETSc is deliberately not set by default
#ifdef  BELFEM_PARDISO
    const SolverType gDefaultSolver = SolverType::PARDISO;
#elif BELFEM_SUITESPARSE
    const SolverType gDefaultSolver = SolverType::UMFPACK;
#elif BELFEM_MUMPS
    const SolverType gDefaultSolver = SolverType::MUMPS;
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
    to_string( const EulerMethod aEulerMethod );

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
    to_string( const SerialReodrdering aSerialReodrdering ) ;

    SerialReodrdering
    serial_reordering( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const ParallelReodrdering aParallelReodrdering ) ;

    ParallelReodrdering
    parallel_reordering( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const BlockLowRanking aBLR );

    BlockLowRanking
    block_low_ranking( const string & aString );

//------------------------------------------------------------------------------
}

#endif //BELFEM_EN_SOLVERENUMS_HPP
