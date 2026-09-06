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

#ifndef BELFEM_CL_GT_HEATPOLYCUSTOM_HPP
#define BELFEM_CL_GT_HEATPOLYCUSTOM_HPP

#include "cl_Vector.hpp"
#include "cl_GT_HeatPoly.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        /**
         * A caloric interval whose exponents are not the NASA-9 ones ( -2 … 4 );
         * used for the fitted cryogenic and hot extrapolations.
         */
        class HeatPolyCustom : public HeatPoly
        {
            // exponents
            const Vector<real> mExponents;

            const uint mNumberOfCoeffs;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPolyCustom(
                    const real           aTmin,
                    const real           aTmax,
                    const real           aEnthalpyConstant,
                    const real           aEntropyConstant,
                    const Vector<real> & aCoefficients,
                    const Vector<real> & aExponents );

//------------------------------------------------------------------------------

            ~HeatPolyCustom() = default;

//------------------------------------------------------------------------------

            real
            Cp( const real T ) const;

//------------------------------------------------------------------------------

            real
            H( const real T ) const;

//------------------------------------------------------------------------------

            real
            S( const real T ) const;

//------------------------------------------------------------------------------

            real
            dSdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            dCpdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            d2CpdT2( const real T ) const;

//------------------------------------------------------------------------------
        };
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_HEATPOLYCUSTOM_HPP
