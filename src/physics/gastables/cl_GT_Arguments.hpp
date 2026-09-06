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

#ifndef BELFEM_CL_GT_ARGUMENTS_HPP
#define BELFEM_CL_GT_ARGUMENTS_HPP

#include "typedefs.hpp"
#include "cl_Arguments.hpp"
#include "GT_globals.hpp"

namespace belfem
{
    namespace gastables
    {
        enum class State
        {
            PrintHelp    = 0,
            PrintTable   = 1,
            PrintUsage   = 2,
            PrintBanner  = 3,
            Undefined    = 4
        };

//------------------------------------------------------------------------------

        /**
         * Command line front end of the gastable executable.
         *
         * Collects the species, the temperature sweep and the molar or mass
         * specific output flag, and reports which of them the user actually
         * set. Defaults are a nitrogen table from 200 to 2000 K in 50 K steps.
         */
        class Arguments : public belfem::Arguments
        {
            State  mState      = State::Undefined;
            real   mTmin       = 200.0;
            real   mTmax       = 2000.0;
            real   mDeltaT     = 50.0;
            string mGasName    = "N2";
            bool   mMolarFlag  = false;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Arguments( int & argc, char * argv[] );

//------------------------------------------------------------------------------

            ~Arguments() = default;

//------------------------------------------------------------------------------

            const State &
            state() const;

//------------------------------------------------------------------------------

            const real &
            T_min() const;
//------------------------------------------------------------------------------

            const real &
            T_max() const;

//------------------------------------------------------------------------------

            const real &
            delta_T() const;

//------------------------------------------------------------------------------

            const string &
            gasname() const;

//------------------------------------------------------------------------------

            bool
            molar() const;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            check_arguments();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GT_ARGUMENTS_HPP
