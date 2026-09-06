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

#ifndef BELFEM_CL_GT_THERMOPOLY_HPP
#define BELFEM_CL_GT_THERMOPOLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
         * One temperature interval of the caloric description of a species.
         *
         * Evaluates the NASA-9 form of NASA RP-1311 Eqs. ( 4.9 ) to ( 4.11 ),
         *
         *     cp/R  =  a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
         *
         * with H and S following by integration, which is where the two
         * constants come in: mEnthalpyConstant is b1 of the record and
         * mEntropyConstant is b2. They are not properties of the polynomial but
         * of the reference state, and RefGas rewrites both when it stitches the
         * intervals together, so an interval taken out of its parent is not
         * self contained.
         *
         * Returns molar quantities: Cp and S in J/(mol K), H in J/mol.
         *
         * The subclasses cover what the tables do not: HeatPolyGlue smooths a
         * junction between two tabulated intervals, HeatPolyCustom carries a
         * fitted interval such as the cryogenic extrapolation, and
         * HeatPolyEmpty stands in for a species with no caloric data at all.
         *
         * @ingroup grp_physics_gastables
         * @see @ref physics_gastables_gastables_usage_guide
         */
        class HeatPoly
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            // minimum temperature of this poly
            real mTmin;

            // maximum temperature
            real mTmax;

            // offset for enthalpy
            real mEnthalpyConstant;

            // offset for entropy
            real mEntropyConstant;

            // coefficients
            const Vector <real> mCoefficients;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPoly(
                    const real aTmin,
                    const real aTmax,
                    const real aEnthalpyConstant,
                    const real aEntropyConstant,
                    const Vector <real> & aCoefficients );

//------------------------------------------------------------------------------

            virtual ~HeatPoly() = default;

//------------------------------------------------------------------------------

            virtual real
            Cp( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            H( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            S( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            dSdT( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            dCpdT( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            d2CpdT2( const real T ) const;

//------------------------------------------------------------------------------

            inline const real &
            T_min() const;

//------------------------------------------------------------------------------

            inline const real &
            T_max() const;

//------------------------------------------------------------------------------

            void
            set_T_min( const real aTmin );

//------------------------------------------------------------------------------

            void
            set_T_max( const real aTmax );

//------------------------------------------------------------------------------

            void
            set_enthalpy_constant( const real & aEnthalpyConstant );

//------------------------------------------------------------------------------

            void
            set_entropy_constant( const real & aEntropyConstant );

//------------------------------------------------------------------------------

            real
            enthalpy_constant() const;
//------------------------------------------------------------------------------

            real
            entropy_constant() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        const real &
        HeatPoly::T_min() const
        {
            return mTmin;
        }

//------------------------------------------------------------------------------

        const real &
        HeatPoly::T_max() const
        {
            return mTmax;
        }

//------------------------------------------------------------------------------


    }
}
#endif //BELFEM_CL_GT_THERMOPOLY_HPP
