//
// Created by Christian Messe on 2019-01-04.
//

#ifndef BELFEM_CL_TIMER_HPP
#define BELFEM_CL_TIMER_HPP

#include <chrono>

namespace belfem
{
//------------------------------------------------------------------------------

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

        inline unsigned int
        stop()
        {
            return ( unsigned int )
                   ( std::chrono::duration_cast<std::chrono::milliseconds>
                   ( std::chrono::high_resolution_clock::now() - mStart ).count() );
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
