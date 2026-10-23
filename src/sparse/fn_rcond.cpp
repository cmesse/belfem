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

#ifdef BELFEM_SUITESPARSE
#include <cholmod.h>
#endif
#include "fn_rcond.hpp"
#include "assert.hpp"

namespace belfem
{
    real
    rcond( SpMatrix & aMatrix )
    {
#ifdef _OPENMP
#ifdef BELFEM_SUITESPARSE
        // a matrix of undefined type exposes no compressed index array
        BELFEM_ERROR( aMatrix.type() == SpMatrixType::CSC || aMatrix.type() == SpMatrixType::CSR,
            "rcond needs a CSC or CSR matrix" );

        BELFEM_ERROR( aMatrix.n_rows() == aMatrix.n_cols(),
            "rcond needs a square matrix, got %lu x %lu",
            ( long unsigned int ) aMatrix.n_rows(), ( long unsigned int ) aMatrix.n_cols() );

        // the view below hands int_t arrays to CHOLMOD as CHOLMOD_INT
        BELFEM_ERROR( sizeof( int_t ) == sizeof( int ),
            "rcond calls the 32-bit CHOLMOD interface, but int_t is %lu bytes",
            ( long unsigned int ) sizeof( int_t ) );

        // CHOLMOD reads zero-based indices; the caller's base is restored below
        const SpMatrixIndexingBase tBase = aMatrix.indexing_base() == 0
                ? SpMatrixIndexingBase::Cpp : SpMatrixIndexingBase::Fortran ;
        aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

        // stype = 1: symmetric, CHOLMOD reads the upper triangle of a CSC
        // matrix. A CSR matrix seen as CSC is its transpose, so CHOLMOD then
        // reads its lower triangle: CSC needs the upper triangle, CSR needs
        // the lower one, and a matrix with both stored works either way
        cholmod_sparse tMatrix ;
        tMatrix.nrow = aMatrix.n_rows();
        tMatrix.ncol = aMatrix.n_cols();
        tMatrix.nzmax = aMatrix.number_of_nonzeros();
        tMatrix.p = aMatrix.pointers();
        tMatrix.i = aMatrix.indices();
        tMatrix.x = aMatrix.data();
        tMatrix.nz = nullptr ;
        tMatrix.z = nullptr ;
        tMatrix.stype = 1 ;
        tMatrix.itype = CHOLMOD_INT ;
        tMatrix.xtype = CHOLMOD_REAL ;
        tMatrix.dtype = CHOLMOD_DOUBLE ;
        tMatrix.sorted = 1;
        tMatrix.packed = 1;

        cholmod_common tCommon ;
        cholmod_start (&tCommon) ;

        // Cholesky factorization: negative status is an error, NOT_POSDEF
        // means A is not SPD, DSMALL is only a warning about a tiny pivot
        cholmod_factor * tFactor = cholmod_analyze( &tMatrix, &tCommon) ;

        bool tOk = tFactor != nullptr && tCommon.status >= CHOLMOD_OK ;

        if ( tOk )
        {
            cholmod_factorize( &tMatrix, tFactor, &tCommon );
            tOk = tCommon.status >= CHOLMOD_OK && tCommon.status != CHOLMOD_NOT_POSDEF ;
        }

        const int tStatus = tCommon.status ;

        // CHOLMOD's rough estimate from the diagonal of the factor, not a
        // 1-norm condition estimate
        const real aCond = tOk ? cholmod_rcond( tFactor, &tCommon ) : BELFEM_QUIET_NAN ;

        if ( tFactor != nullptr )
        {
            cholmod_free_factor ( &tFactor, &tCommon) ;
        }
        cholmod_finish (&tCommon);

        aMatrix.set_indexing_base( tBase );

        BELFEM_ERROR( tOk,
            "rcond: CHOLMOD could not factorize the matrix ( status %i: %s )",
            tStatus,
            tStatus == CHOLMOD_NOT_POSDEF ? "not positive definite"
            : tStatus == CHOLMOD_OUT_OF_MEMORY ? "out of memory"
            : tStatus == CHOLMOD_INVALID ? "invalid input"
            : "CHOLMOD error" );

        return aCond ;
#else
        BELFEM_ERROR( false, "We are not linked against CHOLMOD." );
        return BELFEM_QUIET_NAN ;
#endif
#else
        BELFEM_ERROR( false, "We are not linked against OPENMP." );
        return BELFEM_QUIET_NAN ;
#endif
    }
}