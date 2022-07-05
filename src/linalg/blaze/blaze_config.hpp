//
// Created by Christian Messe on 2018-12-23.
//

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
#else
#define BLAZE_BLAS_MODE 1
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
