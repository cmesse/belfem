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

#ifndef PARPACKTOOLS_HPP
#define PARPACKTOOLS_HPP

#include "typedefs.hpp"

// the two drivers share the info layout and therefore the decoders. USE_PARPACK
// hard-errors without USE_ARPACK ( config/linalg/config_parpack.cmake ), so
// arpacktools is always available here
#include "arpacktools.hpp"

namespace belfem
{
    namespace parpack
    {
        // the info array layout is identical to the serial driver's -- use
        // arpack::gInfoNaupd, arpack::gInfoNumConverged, ... to index it, and
        // arpack::check_naupd() / arpack::check_neupd() to decode it

#ifdef __cplusplus
        extern"C" {
#endif

            /**
             * distributed reverse communication driver around pdnaupd /
             * pdneupd. Every rank owns a contiguous block of ROWS, and PARPACK
             * distributes the Arnoldi basis with it, so the n x ncv basis is
             * split across ranks rather than held whole.
             *
             * COLLECTIVE over MPI_COMM_WORLD. Every rank must call this.
             *
             * Row distribution contract, owned by the caller:
             *   - pointers are LOCAL and one-based, so pointers( 1 ) == 1 and
             *     pointers( nloc + 1 ) == nnz + 1
             *   - indices are GLOBAL and one-based, NOT rebased: a local row
             *     may reference any column of the global vector
             *   - row blocks are contiguous and ordered by rank, so rank p
             *     owns the rows following rank p-1
             *   - every rank owns at least one row ( else info = 102 )
             *
             * Cross-rank uniformity, also owned by the caller and NOT enforced:
             * nglobal, nev, job, tol, maxit and sigma must be identical on every
             * rank -- sigma changes the operator, so a divergent value would put the
             * ranks on different problems.
             * PARPACK runs one algorithm whose control flow is executed
             * redundantly on all ranks, and they stay in lockstep only because
             * they take the same branches. A divergent nev or job changes the
             * trip count inside pdnaupd and HANGS the job rather than failing.
             * ( nglobal is the exception -- it falls out of the row map check
             * and surfaces as info = 101. )
             *
             * info must have arpack::gNumInfoEntries entries; see
             * arpack::check_naupd() for the driver codes 100-105.
             *
             * lambdareal and lambdaimag must have nev + 1 entries and come back
             * identical on every rank. Eigenvalues only -- the driver does not
             * compute Ritz vectors.
             */
            void
            parpack_standard_eigen(
                const int_t * nloc,      // rows owned by this rank
                const int_t * nglobal,   // rows of the whole matrix
                const int_t * nnz,       // nonzeros in the local block
                const real  * values,    // local values
                const int_t * indices,   // GLOBAL column indices
                const int_t * pointers,  // LOCAL row pointers
                const int_t * job,
                const int_t * nev,
                const int_t * ncvmin,   // floor for the Krylov subspace
                const real  * tol,
                const int_t * maxit,
                const real  * sigma,    // 0 for OP = A, else OP = sigma*I - A
                real * lambdareal,
                real * lambdaimag,
                int_t * info ) ;

            /**
             * symmetric counterpart of parpack_standard_eigen. Identical
             * contract, identical row-distribution rules, identical info
             * layout ; pdsaupd / pdseupd instead of pdnaupd / pdneupd.
             * The same cross-rank uniformity requirements apply, sigma
             * included.
             *
             * Flags are decoded by arpack::check_saupd / check_seupd.
             */
            void
            parpack_symmetric_eigen(
                const int_t * nloc,      // rows owned by this rank
                const int_t * nglobal,   // rows of the whole matrix
                const int_t * nnz,       // nonzeros in the local block
                const real  * values,    // local values
                const int_t * indices,   // GLOBAL column indices
                const int_t * pointers,  // LOCAL row pointers
                const int_t * job,
                const int_t * nev,
                const int_t * ncvmin,   // floor for the Krylov subspace
                const real  * tol,
                const int_t * maxit,
                const real  * sigma,    // 0 for OP = A, else OP = sigma*I - A
                real * lambdareal,
                real * lambdaimag,
                int_t * info ) ;

#ifdef __cplusplus
        }
#endif

    }

}
#endif //PARPACKTOOLS_HPP
