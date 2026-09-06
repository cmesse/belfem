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

#ifndef BELFEM_MUMPSTOOLS_HPP
#define BELFEM_MUMPSTOOLS_HPP

#include "typedefs.hpp"
#ifdef BELFEM_MUMPS
#ifdef __cplusplus
extern"C" {
#endif

//------------------------------------------------------------------------------

    /**
     * Creates a new MUMPS solver instance
     * @param[out] aSolverID - unique ID of the created solver (1-based)
     * @param[out] aInfo     - info status (0 = success, -1000 = max solvers reached)
     * @param[in]  aHostIsWorking - MUMPS PAR flag (1 = working host, 0 = host not working)
     * @param[in]  aSymmetryMode  - MUMPS SYM flag (0 = unsymmetric, 1 = SPD, 2 = general symmetric)
     */
    void
    mumpstools_create_solver(
            belfem::int_t    & aSolverID,
            belfem::int_t    & aInfo,
            const belfem::int_t & aHostIsWorking,
            const belfem::int_t & aSymmetryMode );

//------------------------------------------------------------------------------

    /**
     * Frees a MUMPS solver instance and releases its resources
     * @param[in]  aSolverID - ID of the solver to free
     * @param[out] aInfo     - MUMPS INFO array (80 elements)
     */
    void
    mumpstools_free_solver(
            belfem::int_t    & aSolverID,
            belfem::int_t    * aInfo );

//------------------------------------------------------------------------------

    /**
     * Number of live solver instances in the shim's registry ( occupancy ).
     * Rank-uniform because create and free are collective.
     * @param[out] aCount - number of occupied slots
     */
    void
    mumpstools_num_solvers( belfem::int_t & aCount );

//------------------------------------------------------------------------------

    /**
     * Solves a linear system using MUMPS
     * @param[in]  aIParameters - array of integer parameters (see Parameter enum)
     * @param[in]  aRParameters - array of real parameters (compression tolerance, etc.)
     * @param[in]  aN           - matrix dimension
     * @param[in]  aNNZ         - number of non-zeros
     * @param[in]  aNHRS        - number of right-hand side columns
     * @param[in]  aRowIndices  - COO row indices (1-based)
     * @param[in]  aColIndices  - COO column indices (1-based)
     * @param[in]  aValues      - matrix values in COO format
     * @param[out] aX           - solution vector (initially copied from aY)
     * @param[in]  aY           - right-hand side vector
     * @param[out] aInfo        - MUMPS INFO array (80 elements, rank-local)
     * @param[out] aInfoG       - MUMPS INFOG array (80 elements, rank-uniform)
     * @param[out] aRInfoG      - MUMPS RINFOG array (20 elements)
     */
    void
    mumpstools_solve(
        const belfem::int_t  * aIParameters,
        const double * aRParameters,
        const belfem::int_t    & aN,
        const belfem::int_t    & aNNZ,
        const belfem::int_t    & aNHRS,
        const belfem::int_t    * aRowIndices,
        const belfem::int_t    * aColIndices,
        const double * aValues,
        double       * aX,
        const double * aY,
        belfem::int_t * aInfo,
        belfem::int_t * aInfoG,
        double        * aRInfoG ) ;

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
#endif

namespace mumps
{
    enum class Parameter
    {
        SolverID,
        Job,
        MasterRank,
        InfoLevel,
        ErrorAnalysis,
        WorkingHost,
        SymmetryMode,
        SerialReordering,
        ParallelReordering,
        CompressionMode,
        NumRefinementSteps,
        MemoryRelaxation,
        ComputeDeterminant,
        MemoryBudget        // ICNTL(23), MB per process; 0 = not set
    };

//------------------------------------------------------------------------------

    /**
     * what the wrapper does after a workspace failure of the factorization
     */
    enum class WorkspaceAction
    {
        Ladder,   // double ICNTL(14) and factorize again
        Cap,      // hand MUMPS the machine budget as ICNTL(23) and factorize again
        GiveUp    // hand the failure to the soft-fail arms
    };

//------------------------------------------------------------------------------

    /**
     * The retry policy after a factorization returned INFOG(1) = -9, -8,
     * -17, -20 or -19, as a pure function so it can be unit-tested without
     * a MUMPS instance. The MUMPS 5.9 user guide ( sec. 2.11 and the error
     * list ): -9 and -8 ( workarrays S / IS ) ask for a larger ICNTL(14)
     * and "may still occur" once ICNTL(23) is set; -17 and -20 ( the MPI
     * send / reception buffers, sized from ICNTL(14) before anything else
     * once a cap is in place ) ask for a larger ICNTL(14) only -- a cap
     * does not help them; -19 means ICNTL(23) itself cannot be met.
     *
     * @param aInfoG1        rank-uniform INFOG(1)
     * @param aRelax         current ICNTL(14) in percent
     * @param aRelaxCeiling  ladder ceiling for ICNTL(14)
     * @param aBudgetSlot    current ICNTL(23) slot ( 0 = not set )
     * @param aBudgetMB      per-process budget the machine affords ( 0 = unknown )
     * @param aBoundMB       MUMPS's own estimate for the factorization
     *                       ( INFOG(16), or INFOG(36) under BLR; 0 = unknown )
     *
     * Rules, in order:
     *   any other code                              -> GiveUp ( nothing to retry )
     *   -19                                         -> GiveUp
     *   -17 / -20                                   -> Ladder below the ceiling, else GiveUp
     *   -9 / -8, slot 0, budget known and >= bound  -> Cap
     *   slot 0, budget known and <  bound           -> GiveUp ( the machine cannot
     *                                                  hold the estimate; a rung
     *                                                  would ask for more )
     *   otherwise ( slot set, or budget unknown ),
     *   relaxation below the ceiling                -> Ladder
     *   relaxation at the ceiling                   -> GiveUp
     */
    WorkspaceAction
    next_workspace_action(
            const belfem::int_t aInfoG1,
            const belfem::int_t aRelax,
            const belfem::int_t aRelaxCeiling,
            const belfem::int_t aBudgetSlot,
            const belfem::int_t aBudgetMB,
            const belfem::int_t aBoundMB );

}

#endif //BELFEM_MUMPSTOOLS_HPP
