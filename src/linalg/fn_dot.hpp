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
 * @brief Inner products of vectors and matrices.
 * @ingroup grp_linalg
 *
 * The header selects the Armadillo or Blaze implementation at compile time; both
 * present the same interface, documented below.
 *
 * Some overloads are templated on an expression type `ET`. These accept unevaluated
 * backend views and expressions -- a row, a column, or the result of `a + b` -- without
 * first storing them in a temporary `Vector`. They compute the same result as the
 * corresponding concrete overload.
 */

/**
 * @fn template<typename T> auto belfem::dot( const Vector<T> & aA, const Vector<T> & aB )
 * @brief Scalar product of two vectors.
 * @ingroup grp_linalg
 * @param aA left input vector; must have the same length as @p aB
 * @param aB right input vector; must have the same length as @p aA
 * @return scalar sum of the elementwise products
 */

/**
 * @fn template<typename T> auto belfem::dot( const Matrix<T> & aA, const Vector<T> & aB )
 * @brief Scalar product of a matrix and a vector, both read as flat sequences.
 * @ingroup grp_linalg
 *
 * This is **not** a matrix-vector product: it forwards to the backend's own `dot`, which
 * walks both operands as flat element sequences and returns a single number. The matrix
 * must therefore hold exactly as many elements as the vector.
 *
 * @param aA matrix, read element by element
 * @param aB vector with as many entries as @p aA has elements
 * @return scalar sum of the elementwise products
 */

/**
 * @fn template<typename T> auto belfem::dot( const Vector<T> & aA, const Matrix<T> & aB )
 * @brief Scalar product of a vector and a matrix, both read as flat sequences.
 * @ingroup grp_linalg
 *
 * The mirror of the overload above, and equally not a vector-matrix product.
 *
 * @param aA vector with as many entries as @p aB has elements
 * @param aB matrix, read element by element
 * @return scalar sum of the elementwise products
 */

#ifndef BELFEM_FN_DOT_HPP
#define BELFEM_FN_DOT_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_dot.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_dot.hpp"
#endif

#endif //BELFEM_FN_DOT_HPP
