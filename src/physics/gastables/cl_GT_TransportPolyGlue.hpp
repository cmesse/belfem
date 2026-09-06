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

#ifndef BELFEM_CL_GT_TRANSPORTPOLYGLUE_HPP
#define BELFEM_CL_GT_TRANSPORTPOLYGLUE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        /**
         * A short interval bridging two tabulated transport intervals.
         *
         * The CEA intervals meet with a kink, and a solver that differentiates
         * viscosity or conductivity will see it. This class carries the
         * polynomial that RefGas fits across the junction so that the property
         * and its first derivative stay continuous.
         */
        class TransportPolyGlue: public TransportPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPolyGlue(
                    const enum TransportPolyType aType,
                    const real aTmin,
                    const real aTmax,
                    const Vector<real> & aCoefficients );

//------------------------------------------------------------------------------

            ~TransportPolyGlue() = default;

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
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */


#endif //BELFEM_CL_GT_TRANSPORTPOLYGLUE_HPP
