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

/**
 * @file
 * @brief Eigenvalues of a general square matrix.
 * @ingroup grp_linalg
 */

/**
 * @fn int_t belfem::eigen( const Matrix<real> & aMatrix, Vector<real> & aValues, const bool aAbortOnComplex )
 * @brief Eigenvalues of a general (not necessarily symmetric) square matrix.
 * @ingroup grp_linalg
 *
 * A general matrix may have complex eigenvalues, which a real Vector cannot hold. By
 * default that is treated as a caller error and aborts with a BELFEM_ERROR naming the
 * offending index and its imaginary part -- an always-active check, not an assert, so it
 * survives release builds.
 *
 * Pass @p aAbortOnComplex as false to handle the case instead: the complex entries are
 * then written as BELFEM_QUIET_NAN and the return value counts them. This mirrors the
 * `AbortOnError` argument of belfem::gesv and belfem::posv.
 *
 * If the matrix is symmetric, prefer belfem::eigen_sym(), whose eigenvalues are real by
 * construction, so the question cannot arise, and whose solver is faster.
 *
 * @param aMatrix          square matrix; not modified
 * @param aValues          resized to the matrix order and filled with the eigenvalues,
 *                         in the order the backend returns them, which is not sorted
 * @param aAbortOnComplex  abort on the first complex eigenvalue (default); pass false to
 *                         receive a count instead
 * @return the number of eigenvalues found to be complex; 0 when all are real
 *
 * @note An eigenvalue counts as complex when the magnitude of its imaginary part exceeds
 * BELFEM_EPSILON. Both backends use that same threshold.
 */

/**
 * @fn void belfem::eigen_sym( const Matrix<real> & aMatrix, Vector<real> & aValues )
 * @brief Eigenvalues of a symmetric matrix.
 * @ingroup grp_linalg
 *
 * A real symmetric matrix has real eigenvalues, so there is no complex case to handle and
 * no NaN to guard against. Prefer this over belfem::eigen whenever the matrix is known to
 * be symmetric -- which covers most of what a finite-element code produces.
 *
 * The symmetry of @p aMatrix is **not** checked. Only the **upper** triangle is read, on
 * both backends, so passing a non-symmetric matrix yields the eigenvalues of the
 * symmetric matrix implied by that triangle rather than an error.
 *
 * @param aMatrix square, symmetric matrix; not modified
 * @param aValues resized to the matrix order and filled with the eigenvalues in ascending
 *                order
 */

#ifndef BELFEM_FN_EIGEN_HPP
#define BELFEM_FN_EIGEN_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_eigen.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_eigen.hpp"
#endif

#endif //BELFEM_FN_EIGEN_HPP
