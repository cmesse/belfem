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

#ifndef BELFEM_ARMADILLO_HPP
#define BELFEM_ARMADILLO_HPP

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK
#define ARMA_USE_ARPACK
//#define ARMA_USE_SUPERLU
#define ARMA_USE_WRAPPER

// Bounds and size-conformance checks follow the build mode, as Blaze's BLAZE_USER_ASSERTION does
// in blaze_config.hpp, on the predicate that defines BELFEM_ASSERT in assert.hpp: on in debug,
// off in release. ARMA_NO_DEBUG would not do this; it only drops ARMA_DEBUG. ARMA_WARN_LEVEL 1
// silences the per-translation-unit note Armadillo emits for the disabled checks and drops its
// level-2 runtime warnings in release.
#if !defined( NDEBUG ) || defined( DEBUG )
#else
#  ifndef ARMA_DONT_CHECK_CONFORMANCE
#    define ARMA_DONT_CHECK_CONFORMANCE
#  endif
#  ifndef ARMA_WARN_LEVEL
#    define ARMA_WARN_LEVEL 1
#  endif
#endif

#ifdef OMP
#define ARMA_USE_OPENMP
#endif

#ifdef BELFEM_MKL
#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
//#define ARMA_USE_MKL_ALLOC
#endif

#include <armadillo>

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif

#endif //BELFEM_ARMADILLO_HPP
