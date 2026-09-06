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


#ifndef BELFEM_CL_PROGRESSBAR_HPP
#define BELFEM_CL_PROGRESSBAR_HPP

#include "typedefs.hpp"

namespace belfem
{
    /**
     * @brief Progress display for long-running loops.
     *
     * @ingroup grp_core
     * @see @ref core_core_usage_guide
     */
    class Progressbar
    {
        // width of the bar in characters
        const index_t mWidth = 65 ;

        // number of steps the caller will report
        const index_t mNumSteps;

        FILE * mFile;

        // actual progress as part of mNumSteps
        uint mProgress = 0;

        // bar width already drawn, in characters. the bar only redraws when
        // this grows, so the number of frames is bounded by mWidth
        uint mStep = 0;

        // true once the cursor has been handed back to the terminal
        bool mCursorRestored = false;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------


        Progressbar( const uint aNumSteps=100, FILE * aFile = stdout );

//------------------------------------------------------------------------------

        ~Progressbar();

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
    private:
//------------------------------------------------------------------------------

        /**
         * write one frame of the bar and push it out immediately
         */
        void
        draw( const uint aStep, const uint aProgress );

//------------------------------------------------------------------------------

        /**
         * flush the stream, so the frame reaches the user before the
         * next newline ( required under mpirun, where stdout is a pipe )
         */
        void
        flush();

//------------------------------------------------------------------------------

        /**
         * hand the cursor back to the terminal ( idempotent )
         */
        void
        show_cursor();

//------------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_PROGRESSBAR_HPP
