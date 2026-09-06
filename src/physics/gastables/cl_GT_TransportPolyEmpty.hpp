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

#ifndef BELFEM_CL_GT_TRANSPORTPOLYEMPTY_HPP
#define BELFEM_CL_GT_TRANSPORTPOLYEMPTY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        /**
         * Placeholder for a species with no transport data at all.
         *
         * eval(), deval() and ddeval() return zero rather than raising, so that
         * a gas which is only ever used for caloric properties can still be
         * constructed. rawpoly(), drawpoly() and ddrawpoly() have no meaningful
         * value and abort. Every virtual
         * of the base must be overridden here - a signature that does not match
         * silently leaves the base version in place, which then reads
         * coefficients that were never set.
         */
        class TransportPolyEmpty: public TransportPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPolyEmpty( const enum TransportPolyType aType  );

//------------------------------------------------------------------------------

            ~TransportPolyEmpty() = default;

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

            real
            eval( const real T ) const;

//------------------------------------------------------------------------------

            real
            deval( const real T ) const;

//------------------------------------------------------------------------------

            real
            ddeval( const real T ) const;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_CL_GT_TRANSPORTPOLYEMPTY_HPP
