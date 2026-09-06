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


#ifndef BELFEM_BLAZE_CONFIG_HPP
#define BELFEM_BLAZE_CONFIG_HPP

//------------------------------------------------------------------------------
// DEBUG SETTINGS
//------------------------------------------------------------------------------

#if !defined( NDEBUG ) || defined( DEBUG )
#ifndef BLAZE_USER_ASSERTION
#define BLAZE_USER_ASSERTION 1
#endif
#ifndef BLAZE_USE_DEBUG_MODE
#define BLAZE_USE_DEBUG_MODE 1
#endif
#ifndef BLAZE_USE_FUNCTION_TRACES
#define BLAZE_USE_FUNCTION_TRACES 0
#endif
#else
#ifndef BLAZE_USER_ASSERTION
#define BLAZE_USER_ASSERTION 0
#endif
#ifndef BLAZE_USE_DEBUG_MODE
#define BLAZE_USE_DEBUG_MODE 0
#endif
#ifndef BLAZE_USE_FUNCTION_TRACES
#define BLAZE_USE_FUNCTION_TRACES 0
#endif
#endif

//------------------------------------------------------------------------------
// PARALLEL MODE
//------------------------------------------------------------------------------

#ifdef BELFEM_MPI
#define BLAZE_MPI_PARALLEL_MODE 1
#else
#define BLAZE_MPI_PARALLEL_MODE 0
#endif

//------------------------------------------------------------------------------
// STORAGE ORDER
//------------------------------------------------------------------------------
#define BLAZE_DEFAULT_STORAGE_ORDER blaze::columnMajor

//------------------------------------------------------------------------------
// BLAS
//------------------------------------------------------------------------------

#ifdef BELFEM_NETLIB
#define BLAZE_BLAS_MODE 0
#elif BELFEM_ACCELLERATE
#define BLAZE_BLAS_MODE 1
#define BLAZE_BLAS_INCLUDE_FILE <vecLib/cblas.h>
#elif BELFEM_MKL
#define BLAZE_BLAS_MODE 1
#define BLAZE_BLAS_INCLUDE_FILE <mkl_cblas.h>
#endif

#if !defined( BLAZE_BLAS_MODE )
#define BLAZE_BLAS_MODE 1
#endif

// keep Blaze's BLAS integer in sync with belfem::int_t ( SCLS builds all
// third party libraries with one consistent integer width ); without this,
// blaze::blas_int_t stays 32 bit under an ILP64 build and the static_assert
// in lapacktools.hpp fires
#ifdef BELFEM_INT64
#define BLAZE_BLAS_IS_64BIT 1
#endif

//------------------------------------------------------------------------------
// COMPILER
//------------------------------------------------------------------------------

#ifdef BELFEM_CLANG
#define _BLAZE_SYSTEM_COMPILER_H_
#define BLAZE_GNU_COMPILER 0
#define BLAZE_CLANG_COMPILER 1
#define BLAZE_MSC_COMPILER 0
#define BLAZE_INTEL_COMPILER 0
#elif BELFEM_GCC
#define _BLAZE_SYSTEM_COMPILER_H_
#define BLAZE_GNU_COMPILER 1
#define BLAZE_CLANG_COMPILER 0
#define BLAZE_MSC_COMPILER 0
#define BLAZE_INTEL_COMPILER 0
#elif BELFEM_INTEL
#define _BLAZE_SYSTEM_COMPILER_H_
#define BLAZE_GNU_COMPILER 0
#define BLAZE_CLANG_COMPILER 0
#define BLAZE_MSC_COMPILER 0
#define BLAZE_INTEL_COMPILER 1
#endif
#include <blaze/util/typetraits/HasSize.h>
//------------------------------------------------------------------------------

#endif //BELFEM_BLAZE_CONFIG_HPP
