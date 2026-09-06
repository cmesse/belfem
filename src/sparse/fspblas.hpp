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

#ifndef BELFEM_FSPBLAS_HPP
#define BELFEM_FSPBLAS_HPP
#include "typedefs.hpp"

// sparse matrix-vector kernels ( splinalg.f90 ); base-agnostic - the
// caller passes the indexing base ( pointers[ 0 ], 0 or 1 ) and the
// kernel shifts internally, so a matvec never rewrites the matrix.
// They overwrite y, so scaling with alpha and beta is handled by the caller

#ifdef __cplusplus
extern"C" {
#endif

        void
        matvec_csr (
            const belfem::int_t * n,
            const belfem::int_t * m,
            const belfem::int_t * nnz,
            const belfem::real  * values,
            const belfem::int_t * indices,
            const belfem::int_t * pointers,
            const belfem::real  * x,
                  belfem::real  * y,
            const belfem::int_t * base
        );

        void
        matvec_csc (
            const belfem::int_t * n,
            const belfem::int_t * m,
            const belfem::int_t * nnz,
            const belfem::real  * values,
            const belfem::int_t * indices,
            const belfem::int_t * pointers,
            const belfem::real  * x,
                  belfem::real  * y,
            const belfem::int_t * base
        );

#ifdef __cplusplus
}
#endif

// for the extended product y = alpha * op(A) * x + beta * y, MKL's
// Inspector-Executor sparse BLAS is used when available ( the classic
// NIST-style mkl_dcscmm / mkl_dcsrmm was removed in oneMKL 2026 )
#ifdef BELFEM_MKL
#include <mkl_spblas.h>
#endif

#endif //BELFEM_FSPBLAS_HPP
