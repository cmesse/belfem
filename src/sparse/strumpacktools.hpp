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

#include "typedefs.hpp"

#include "cl_SolverDistMatrix.hpp"

#ifdef BELFEM_STRUMPACK
#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#elif BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#endif

#include <StrumpackSparseSolver.hpp>
#include <StrumpackSparseSolverMPIDist.hpp>

namespace belfem
{
    namespace sparse
    {
        typedef strumpack::SPOptions< real >   StrumpackOptions;
        typedef strumpack::ReturnCode          StrumpackReturnCode;

        typedef DistMatrixCSR< int_t > StrumpackCSR ;
        typedef DistMatrixAIJ< int_t > StrumpackAIJ ;
    }
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_CLANG
#pragma clang diagnostic pop
#endif
#else
namespace belfem
{
    namespace sparse
    {
        typedef int  StrumpackOptions ;
        typedef int  StrumpackReturnCode ;
    }
}
#endif

#include "cl_SolverParameters.hpp"

namespace belfem
{
    namespace sparse
    {
        string
	    strumpack_message( const StrumpackReturnCode aCode ) ;

        void
        set_strumpack_options( const SolverParameters & aParams,
                               StrumpackOptions & aOpts,
                               const index_t aNumRows );
    }
}
