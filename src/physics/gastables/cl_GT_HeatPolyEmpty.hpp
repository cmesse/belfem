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

#ifndef BELFEM_CL_GT_HEATPOLYEMPTY_HPP
#define BELFEM_CL_GT_HEATPOLYEMPTY_HPP

#include "GT_globals.hpp"
#include "cl_GT_HeatPoly.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
          * An empty polynomial that returns only zeros.
          * Creted if no data has been found in thermo.inp
          */
        class HeatPolyEmpty : public HeatPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPolyEmpty() :
                HeatPoly( 0.0, gTmax, 0.0, 0.0, Vector<real> ( 1, 0.0 ) )
            {

            }

            ~HeatPolyEmpty() = default;

//------------------------------------------------------------------------------

            real
            Cp( const real T ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            H( const real T ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            S( const real T ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            dSdT( const real T ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            dCpdT( const real T ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            d2CpdT2( const real T ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_GT_HEATPOLYEMPTY_HPP
