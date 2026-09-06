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

#ifndef ARPACKTOOLS_HPP
#define ARPACKTOOLS_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace arpack
    {
        //! number of entries the info array of the eigen drivers expects.
        //! arpacktools.f90 and parpacktools.f90 share this layout on purpose,
        //! so check_naupd() / check_neupd() decode both
        constexpr index_t gNumInfoEntries = 7 ;

        //! index of the dnaupd flag within that array
        constexpr index_t gInfoNaupd = 0 ;

        //! index of the dneupd flag within that array
        constexpr index_t gInfoNeupd = 1 ;

        //! index of the number of converged Ritz values
        constexpr index_t gInfoNumConverged = 2 ;

        //! index of the number of restart iterations taken
        constexpr index_t gInfoNumIterations = 3 ;

        //! index of the number of OP*x operations
        constexpr index_t gInfoNumOperations = 4 ;

        //! index of the number of re-orthogonalizations
        constexpr index_t gInfoNumReorthogonalizations = 5 ;

        //! index of the Krylov subspace size actually used
        constexpr index_t gInfoSubspaceSize = 6 ;

#ifdef __cplusplus
        extern"C" {
#endif

            /**
             * serial reverse communication driver around dnaupd / dneupd.
             * The whole matrix lives on the calling rank.
             *
             * info must have gNumInfoEntries entries and is filled with
             *   0 : dnaupd flag, or a driver code >= 100 ( see check_naupd )
             *   1 : dneupd flag
             *   2 : NCONV,  number of converged Ritz values
             *   3 : MXITER, restart iterations taken
             *   4 : NUMOP,  number of OP*x operations
             *   5 : NUMREO, number of re-orthogonalizations
             *   6 : NCV,    Krylov subspace size actually used
             *
             * lambdareal and lambdaimag must have nev + 1 entries.
             * Eigenvalues only -- the driver does not compute Ritz vectors.
             */
            void
            arpack_standard_eigen(
                const int_t * n,
                const int_t * nnz,
                const real * values,
                const int_t * indices,
                const int_t * pointers,
                const int_t * job,
                const int_t * nev,
                const int_t * ncvmin,
                const real  * tol,
                const int_t * maxit,
                const real  * sigma,    // 0 for OP = A, else OP = sigma*I - A
                real * lambdareal,
                real * lambdaimag,
                int_t * info ) ;

            /**
             * symmetric counterpart of arpack_standard_eigen: dsaupd / dseupd
             * instead of dnaupd / dneupd. Same argument list, same info
             * layout, same fold argument -- and the SAME BUFFER SIZES, so a
             * caller can switch between the two on the strength of the
             * matrix alone.
             *
             * lambdaimag comes back as exact zeros: the spectrum of a
             * symmetric matrix is real.
             *
             * The flags in info are decoded by check_saupd / check_seupd, NOT
             * by the nonsymmetric pair -- the error tables differ.
             */
            void
            arpack_symmetric_eigen(
                const int_t * n,
                const int_t * nnz,
                const real * values,
                const int_t * indices,
                const int_t * pointers,
                const int_t * job,
                const int_t * nev,
                const int_t * ncvmin,
                const real  * tol,
                const int_t * maxit,
                const real  * sigma,    // 0 for OP = A, else OP = sigma*I - A
                real * lambdareal,
                real * lambdaimag,
                int_t * info ) ;

            /**
             * ONE step of a shift-invert ( mode 3 ) reverse-communication
             * loop. The loop itself lives in C++, because OP = A^-1 needs the
             * sparse solver -- see run_shift_invert().
             *
             * Only numeric arguments cross this boundary: BMAT = 'I' and
             * WHICH = 'LM' stay Fortran-side, so no character descriptor is
             * passed. Every array belongs to the caller and MUST NOT be
             * touched between the last step call and the extract call.
             *
             * Service the returned ido as:
             *   -1 : y = A^-1 * x, x at workd[ ipntr[0] - 1 ]
             *    1 : y = A^-1 * x, x at workd[ ipntr[2] - 1 ]   <-- ipntr(3)!
             *    2 : y = x                       ( legal, not expected here )
             *   99 : done
             * result always to workd[ ipntr[1] - 1 ]. The ipntr values are
             * ONE-based Fortran indices
             */
            void
            arpack_si_step(
                int_t * ido,
                const int_t * n,
                const int_t * nev,
                const int_t * ncv,
                real  * tol,           // dsaupd may overwrite it
                real  * resid,         // [ n ]
                real  * vectors,       // [ ldv * ncv ]
                const int_t * ldv,
                int_t * iparam,        // [ 11 ]
                int_t * ipntr,         // [ 11 ]
                real  * workd,         // [ 3 * n ]
                real  * workl,         // [ lworkl ]
                const int_t * lworkl,
                int_t * info ) ;

            /**
             * extract the eigenvalues after arpack_si_step returned ido = 99.
             *
             * The values come back ALREADY transformed to eigenvalues of A --
             * dseupd applies lambda = 1/theta + sigma itself. Do NOT take a
             * reciprocal afterwards; that would invert them a second time.
             *
             * lambdareal must hold nev entries. Eigenvalues only: no Ritz
             * vectors are computed
             */
            void
            arpack_si_extract(
                const real  * sigma,
                const int_t * n,
                const int_t * nev,
                const int_t * ncv,
                real  * tol,
                real  * resid,
                real  * vectors,
                const int_t * ldv,
                int_t * iparam,
                int_t * ipntr,
                real  * workd,
                real  * workl,
                const int_t * lworkl,
                real  * lambdareal,    // [ nev ], eigenvalues of A
                int_t * info ) ;

#ifdef __cplusplus
        }
#endif

        /**
         * decode the dnaupd / pdnaupd flag, info( gInfoNaupd ).
         * throws on a negative flag except -9999 ( no Arnoldi factorization
         * could be built ), which is reported like the non-converged exits
         * and left to the converged count; warns on a non-converged exit.
         *
         * Flags >= 100 are raised by the BELFEM drivers rather than by
         * ARPACK, and mean nothing was computed:
         *   100 unexpected reverse communication request
         *   101 the local row blocks do not sum to the global row count,
         *       or the global row count differs between ranks
         *   102 at least one rank owns no rows
         *   103 the local nonzero count contradicts the row pointers
         *   104 the job flag is neither 0 nor 1
         *   105 an MPI call failed while the row map was built
         * 101-105 can only come from the distributed driver.
         */
        void
        check_naupd( const int_t aInfo );

        /**
         * decode the dneupd flag, info( gInfoNeupd ).
         * dnaupd and dneupd do NOT share an error table, so they must not
         * share a decoder either
         */
        void
        check_neupd( const int_t aInfo );

        /**
         * decode the dsaupd / pdsaupd flag, info( gInfoNaupd ).
         *
         * dsaupd does NOT share dnaupd's table, which is why this exists:
         * -5 admits 'LA', 'SA' and 'BE' that the nonsymmetric routine
         * rejects, -13 means something different. info = 3 ( no shifts
         * applied ) is documented by both and decoded as a warning here as
         * well. Decoding a symmetric run with check_naupd would misreport.
         *
         * The driver codes >= 100 have the same meaning as in check_naupd.
         */
        void
        check_saupd( const int_t aInfo );

        /**
         * decode the dseupd / pdseupd flag, info( gInfoNeupd ).
         * Separate from check_neupd for the same reason check_saupd is
         * separate from check_naupd
         */
        void
        check_seupd( const int_t aInfo );

    }

}
#endif //ARPACKTOOLS_HPP
