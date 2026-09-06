/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_GRAPH_TYPEDEFS_HPP
#define BELFEM_GRAPH_TYPEDEFS_HPP


// compiler settings
#ifdef BELFEM_GCC
#ifndef __GNUC__
#define __GNUC__
#endif
#elif BELFEM_INTEL
#ifndef __ICC
#define __ICC
#endif
#endif

#ifdef BELFEM_METIS
#include <metis.h>
#endif


#ifdef BELFEM_PARMETIS
#include <parmetis.h>
#endif

#ifdef BELFEM_SCOTCH
#include <iostream>
#include <cstdint>
#include <scotch.h>
#endif

// bugfix to avoid error with BLAZE
#ifdef BELFEM_BLAZE
#undef abs
#undef iabs
#endif

namespace belfem
{
#ifdef BELFEM_METIS
    typedef ::idx_t metis_t;
#else
    typedef int   metis_t;
#endif
#ifdef BELFEM_SCOTCH
    typedef SCOTCH_Num scotch_t;
#else
    typedef int scotch_t;
#endif
}
#endif //BELFEM_GRAPH_TYPEDEFS_HPP