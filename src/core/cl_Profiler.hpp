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

#ifndef BELFEM_CL_PROFILER_HPP
#define BELFEM_CL_PROFILER_HPP



#include "typedefs.hpp"

namespace belfem
{
    /**
     * @brief gperftools CPU profiling with Callgrind-format output.
     *
     * @ingroup grp_core
     * @see @ref core_core_usage_guide
     */
    class Profiler
    {
        string mLogFile ;
        string mCallgrindFile ;
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Profiler( const string aLogFile="profiler.log");

//------------------------------------------------------------------------------

        ~Profiler() = default ;

//------------------------------------------------------------------------------

        void
        start();

//------------------------------------------------------------------------------

        void
        stop();

//------------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_PROFILER_HPP
