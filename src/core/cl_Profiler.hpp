//
// Created by christian on 7/22/21.
//

#ifndef BELFEM_CL_PROFILER_HPP
#define BELFEM_CL_PROFILER_HPP



#include "typedefs.hpp"

namespace belfem
{
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
