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

#ifndef BELFEM_TYPEDEFS_HPP
#define BELFEM_TYPEDEFS_HPP

#include <string>
#include <limits>
#include <complex>
#include <array>
#include <cstdint>

namespace belfem
{
//------------------------------------------------------------------------------

    typedef size_t              size_t;
    typedef std::string         string;

    typedef int                      sint;
    typedef long int                 lsint;
    typedef unsigned int             uint;
    typedef unsigned char            uchar;
    typedef short unsigned int       suint;
    typedef long unsigned int        luint;
    typedef long long unsigned int   lluint;

    typedef double                   real;
    typedef std::complex<real>       cplx;

//------------------------------------------------------------------------------

    typedef unsigned int             id_t;
    typedef int                      proc_t;
    typedef long long unsigned int   key_t;
    typedef __uint128_t              key128_t;

//------------------------------------------------------------------------------
#ifdef BELFEM_INT64
    typedef int64_t             int_t ;
    typedef uint64_t            index_t;
#else
    typedef int32_t             int_t ;
    typedef uint32_t            index_t;
#endif

//------------------------------------------------------------------------------

    constexpr index_t gNoIndex = std::numeric_limits<index_t>::max();
    constexpr id_t    gNoID    = std::numeric_limits<id_t>::max();
    constexpr proc_t  gNoOwner = std::numeric_limits<proc_t>::max();

    // todo: probably better to move to constants
    constexpr real    gTfreeze = 273.15 ;
    constexpr real    gTref    = 288.15 ;
    constexpr real    gTroom   = 293.15 ;
    constexpr real    gTmin    = 1.0 ;

    constexpr real    gFinDiffDeltaT      = 0.001 ;
    constexpr real    gFinDiffDeltaB      = 0.001 ;
    constexpr real    gFinDiffDeltaAngle  = 1.74532925199433e-3 ;

    // L, M, T, I, theta, N, J
    typedef std::array< real, 7 >   unit ;

    typedef std::pair< real, unit > value ;

//------------------------------------------------------------------------------
#define BELFEM_UCHAR_MAX      std::numeric_limits<uchar>::max()
#define BELFEM_SINT_MAX      std::numeric_limits<sint>::max()
#define BELFEM_SUINT_MAX     std::numeric_limits<suint>::max()
#define BELFEM_UINT_MAX      std::numeric_limits<uint>::max()
#define BELFEM_REAL_MAX      std::numeric_limits<real>::max()
#define BELFEM_REAL_MIN      std::numeric_limits<real>::min()
#define BELFEM_INT_MAX       std::numeric_limits<int>::max()
#define BELFEM_LUINT_MAX     std::numeric_limits<luint>::max()
#define BELFEM_KEY_MAX       std::numeric_limits<key_t>::max()
#define BELFEM_SIGNALING_NAN std::numeric_limits<real>::signaling_NaN()
#define BELFEM_QUIET_NAN     std::numeric_limits<real>::quiet_NaN()
#define BELFEM_INFINITY      std::numeric_limits<real>::infinity()

    constexpr   real BELFEM_EPSILON = 10 * std::numeric_limits<real>::epsilon();
    constexpr   real BELFEM_EPS     =      std::numeric_limits<real>::epsilon();
    constexpr   real BELFEM_MESH_EPSILON = 1e-9 ;

//------------------------------------------------------------------------------
}
#endif //BELFEM_TYPEDEFS_HPP
