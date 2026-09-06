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

#ifndef BELFEM_CL_GT_INPUTALPHA_HPP
#define BELFEM_CL_GT_INPUTALPHA_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_GT_GasData.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        // reads the alpha coefficients for a cubic gas
        /**
         * Reader for cubicalpha.inp, the fitted alpha function coefficients of
         * Mahmoodi and Sedigh ( doi 10.1016/j.fluid.2016.12.015 ).
         *
         * Records are keyed on CAS number rather than on species label, which
         * is what lets the join work across the different nomenclatures of the
         * source tables. read_data applies the structural factors of the
         * paper's Eq. ( 14 ) - c1 = 2 C1, c2 = -C2^2, c3 = (2/3) C3^3 - so the
         * evaluator downstream sees a plain polynomial and the transform lives
         * in exactly one place.
         *
         * It also overwrites the critical temperature and pressure of the
         * GasData with the values of this file, so that the alpha function and
         * the critical point it was fitted against stay consistent.
         *
         * A species with no record here simply keeps has_cubic() false and
         * falls back to the Coquelet correlation, which needs only the acentric
         * factor. That is a designed fallback, not a failure.
         */
        class InputAlpha : public Ascii
        {
            // map to line in buffer
            Map <string, uint> mMap;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            InputAlpha( const string & aPath );

//----------------------------------------------------------------------------

            ~InputAlpha() = default;

//----------------------------------------------------------------------------

            bool
            entry_exists( GasData * aData );

//----------------------------------------------------------------------------

            void
            read_data( GasData * aData );

//----------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_CL_GT_INPUTALPHA_HPP
