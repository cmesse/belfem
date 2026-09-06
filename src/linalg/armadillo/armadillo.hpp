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
