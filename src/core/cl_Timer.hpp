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

#ifndef BELFEM_CL_TIMER_HPP
#define BELFEM_CL_TIMER_HPP

#include <chrono>

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief High-resolution wall-clock timing.
     *
     * @ingroup grp_core
     * @see @ref core_core_usage_guide
     */
    class Timer
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> mStart;
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        inline Timer() : mStart( std::chrono::high_resolution_clock::now() ) {}

//------------------------------------------------------------------------------

        ~Timer() = default;

//------------------------------------------------------------------------------

        inline uint64_t
        stop()
        {
            return ( unsigned int )
                   ( std::chrono::duration_cast<std::chrono::milliseconds>
                   ( std::chrono::high_resolution_clock::now() - mStart ).count() );
        }

//------------------------------------------------------------------------------

        inline uint64_t
        next()
        {
            unsigned int aTime = this->stop();
            mStart = std::chrono::high_resolution_clock::now();
            return aTime ;
        }

//------------------------------------------------------------------------------

        inline void
        reset()
        {
            mStart = std::chrono::high_resolution_clock::now();
        }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_TIMER_HPP
