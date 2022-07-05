//
// Created by Christian Messe on 17.11.18.
//

#ifndef BELFEM_TYPEDEFS_HPP
#define BELFEM_TYPEDEFS_HPP

#include <string>
#include <limits>
#include <complex>
#include <array>

namespace belfem
{
//------------------------------------------------------------------------------

    typedef size_t              size_t;
    typedef std::string         string;

    typedef int                 sint;
    typedef long int            lsint;
    typedef long unsigned int   luint;

    typedef unsigned int        uint;
    typedef double              real;
    typedef std::complex<real>  cplx;

//------------------------------------------------------------------------------

    typedef unsigned int        id_t;
    typedef unsigned int        index_t;
    typedef int                 proc_t;

//------------------------------------------------------------------------------

    const index_t               gNoIndex = std::numeric_limits<index_t>::max();
    const id_t                  gNoID    = std::numeric_limits<id_t>::max();
    const proc_t                gNoOwner = std::numeric_limits<proc_t>::max();

    // L, M, T, I, theta, N, J
    typedef std::array< real, 7 >   unit ;

    typedef std::pair< real, unit > value ;

//------------------------------------------------------------------------------

#define BELFEM_SINT_MAX      std::numeric_limits<sint>::max()
#define BELFEM_UINT_MAX      std::numeric_limits<uint>::max()
#define BELFEM_REAL_MAX      std::numeric_limits<real>::max()
#define BELFEM_REAL_MIN      std::numeric_limits<real>::min()
#define BELFEM_INT_MAX       std::numeric_limits<int>::max()
#define BELFEM_LUINT_MAX     std::numeric_limits<luint>::max()

#define BELFEM_SIGNALING_NAN std::numeric_limits<real>::signaling_NaN()
#define BELFEM_QUIET_NAN     std::numeric_limits<real>::quiet_NaN()
    const   real BELFEM_EPSILON = 10 * std::numeric_limits<real>::epsilon();
    const   real BELFEM_EPS     =      std::numeric_limits<real>::epsilon();
    const   real BELFEM_MESH_EPSILON = 1e-9 ;

//------------------------------------------------------------------------------
}
#endif //BELFEM_TYPEDEFS_HPP