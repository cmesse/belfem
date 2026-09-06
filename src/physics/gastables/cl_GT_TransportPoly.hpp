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

#ifndef BELFEM_CL_GT_TRANSPORTPOLY_HPP
#define BELFEM_CL_GT_TRANSPORTPOLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        enum class TransportPolyType
        {
            VISCOSITY,
            CONDUCTIVITY,
            UNDEFINED
        };

//------------------------------------------------------------------------------

        enum class TransportPolyKind
        {
            DEFAULT,
            CUSTOM,
            GLUE,
            EMPTY,
            UNDEFINED
        };
//------------------------------------------------------------------------------

        /**
         * One temperature interval of the transport description of a species.
         *
         * Evaluates the CEA correlation of NASA RP-1311 Eq. ( 5.1 ),
         *
         *     ln( X )  =  A ln T + B/T + C/T^2 + D
         *
         * for either viscosity or thermal conductivity, selected by mType.
         *
         * mScale carries the unit conversion and is the reason the two kinds
         * cannot share one object: the tabulated coefficients give viscosity in
         * micropoise and conductivity in microwatt per centimetre kelvin, so
         * the factors to SI are 1e-7 and 1e-4 respectively. Everything this
         * class returns is already SI - Pa s and W/(m K).
         *
         * @ingroup grp_physics_gastables
         * @see @ref physics_gastables_gastables_usage_guide
         */
        class TransportPoly
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            const TransportPolyType mType;
            const TransportPolyKind mKind ;

            // unit scale
            const real mScale;

            // minimum temperature
            real mTmin;

            // maximum temperature
            real mTmax;

            // coefficients
            const Vector <real> mCoefficients;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPoly(
                    const enum TransportPolyType aType,
                    const real aTmin,
                    const real aTmax,
                    const Vector <real> & aCoefficients,
                    const TransportPolyKind aKind = TransportPolyKind::DEFAULT );

            virtual ~TransportPoly() = default;

//------------------------------------------------------------------------------

            virtual real
            rawpoly( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            drawpoly( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            ddrawpoly( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            eval( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            deval( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            ddeval( const real T ) const;

//------------------------------------------------------------------------------

            inline const real &
            T_min() const;

//------------------------------------------------------------------------------

            inline const real &
            T_max() const;

//------------------------------------------------------------------------------

            TransportPolyType
            type() const;

//------------------------------------------------------------------------------

            void
            set_T_min( const real aTmin );

//------------------------------------------------------------------------------

            void
            set_T_max( const real aTmax );

//------------------------------------------------------------------------------

            TransportPolyKind
            kind() const ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        const real &
        TransportPoly::T_min() const
        {
            return mTmin;
        }

//------------------------------------------------------------------------------

        const real &
        TransportPoly::T_max() const
        {
            return mTmax;
        }
//------------------------------------------------------------------------------
        inline TransportPolyType
        TransportPoly::type() const
        {
            return mType;
        }

//------------------------------------------------------------------------------

        inline TransportPolyKind
        TransportPoly::kind() const
        {
            return mKind ;
        }


//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_TRANSPORTPOLY_HPP
