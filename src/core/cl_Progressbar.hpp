//
// Created by Christian Messe on 03.11.19.
//

#ifndef BELFEM_CL_PROGRESSBAR_HPP
#define BELFEM_CL_PROGRESSBAR_HPP

#include "typedefs.hpp"

namespace belfem
{
    class Progressbar
    {
        // Width of this bar
        const index_t mNumSteps;

        FILE * mFile;

        // actual progress as part of mNumSteps
        uint mProgress = 0;

        // progess in percent
        uint mStep = 0;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------


        Progressbar( const uint aNumSteps=100, FILE * aFile = stdout );

//------------------------------------------------------------------------------

        ~Progressbar() = default;

//------------------------------------------------------------------------------

        void
        reset();

//------------------------------------------------------------------------------

        /**
         * prescribe a step
         * @param aProgress
         */
        void
        step( const uint & aProgress );

//------------------------------------------------------------------------------

        /**
         * just go to next setp
         */
         void
         step();

//------------------------------------------------------------------------------

        void
        finish();

//------------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_PROGRESSBAR_HPP
