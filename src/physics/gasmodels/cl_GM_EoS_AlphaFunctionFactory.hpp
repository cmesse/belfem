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

#ifndef BELFEM_CL_GM_ALPHAFUNCTIONFACTORY_HPP
#define BELFEM_CL_GM_ALPHAFUNCTIONFACTORY_HPP

#include "cl_GM_EoS_AlphaFunction.hpp"

namespace belfem
{
    namespace gastables
    {
        class GasData;
    }

    namespace gasmodels
    {

//----------------------------------------------------------------------------
        /**
         * Creates the alpha function of a cubic equation of state for one
         * species.
         *
         * The alpha function carries the whole temperature dependence of the
         * attraction term, so which one is picked matters more than the choice
         * between SRK and PR. The selection is by data availability, best
         * first: a fitted Mahmoodi and Sedigh form when cubicalpha.inp has a
         * record for the CAS number, otherwise the Coquelet, Chapoy and Richon
         * form - from a table of 22 hard coded species where one exists, and
         * from the generalized acentric factor correlation where it does not.
         *
         * The plain Soave and Peng Robinson 1978 factories are kept for
         * benchmarking against the classical forms and are not wired into
         * EoS_Cubic.
         *
         * Every function returns a heap object that the caller owns.
         */
        class AlphaFunctionFactory
        {
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

             AlphaFunctionFactory() = default;

            ~AlphaFunctionFactory() = default;

//----------------------------------------------------------------------------

            AlphaFunction *
            create_empty();

//----------------------------------------------------------------------------

            AlphaFunction *
            create_srk( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_pr78( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_ccr_pr( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_ccr_mc_srk( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_pm_srk(  const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_pm_pr(  const gastables::GasData * aData );

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_ALPHAFUNCTIONFACTORY_HPP
