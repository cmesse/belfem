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

#ifndef BELFEM_CL_GT_TRANSPORTPOLYCUSTOM_HPP
#define BELFEM_CL_GT_TRANSPORTPOLYCUSTOM_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        /**
         * A transport interval whose exponents are not the CEA ones.
         *
         * Used where the tabulated correlation does not reach, in particular
         * for the viscosity synthesized from the Lucas correlation for species
         * that have a critical point but no trans.inp record. Carries its own
         * exponent vector instead of the fixed ln T, 1/T, 1/T^2, 1 of the base.
         */
        class TransportPolyCustom: public TransportPoly
        {
            const Vector<real> mExponents;
            const uint mNumberOfCoeffs;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPolyCustom(
                    const enum TransportPolyType aType,
                    const real aTmin,
                    const real aTmax,
                    const Vector<real> & aCoefficients,
                    const Vector<real> & aExponents );

//------------------------------------------------------------------------------

            ~TransportPolyCustom() = default;

//------------------------------------------------------------------------------

            real
            rawpoly( const real T ) const;

//------------------------------------------------------------------------------

            real
            drawpoly( const real T ) const;

//------------------------------------------------------------------------------

            real
            ddrawpoly( const real T ) const;

//------------------------------------------------------------------------------
        };
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_TRANSPORTPOLYCUSTOM_HPP
