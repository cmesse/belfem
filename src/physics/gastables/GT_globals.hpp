//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_GT_GLOBALS_HPP
#define BELFEM_GT_GLOBALS_HPP

#include "typedefs.hpp"
#include "constants.hpp"

namespace belfem
{
    namespace gastables
    {
        const real gTref = BELFEM_TREF;
        const real gPref = BELFEM_PREF;
        const real gTmax = 13000.0; // << -- increasing this might fail tests
                                    // default value is 6000

        const real gTmin = 100.0; // << -- eg. for Prandtl-Meyer expansion
        const real gPmin = 1e-4 ; // << -- eg. for Prandtl-Meyer expansion

        // default value is 6000
        const real gDeltaT = 5.0;

        // calculate number of samples
        const uint gNumberOfSplinePoints = ( uint ) gTmax / gDeltaT + 1;

    }
}

#endif //BELFEM_GT_GLOBALS_HPP
