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

#ifndef BELFEM_CL_GM_EOS_CUBIC_HPP
#define BELFEM_CL_GM_EOS_CUBIC_HPP

#include "typedefs.hpp"
#include "cl_Spline.hpp"
#include "cl_GM_EoS.hpp"
#include "en_GM_GasModel.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace gasmodels
    {
        class AlphaFunction;

//----------------------------------------------------------------------------

        /**
         * cubic equation of state
         *
         * \f$ p = \frac{R \, T }{ v - b} - \frac{a \, \alpha}{ \left( v-b \, r_1\right) \, \left( v-b \, r_2\right)} \f$
         *
         * @ingroup grp_physics_gasmodels
         * @see @ref physics_gasmodels_gasmodels_usage_guide
         */
        class EoS_Cubic : public EoS
        {
//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------
            // equation constants
            real mR1;
            real mR2;
            real mOmegaA;
            real mOmegaB;

            Vector< real > mAc;
            Vector< real > mBc;

            // values of alpha, recomputed by eval_a for the current T
            mutable Vector< real > mA;
            mutable Vector< real > mdAdT;
            mutable Vector< real > md2AdT2;

            // b, b*r1, b*r2 and b*( r2 - r1 ), mixed in remix()
            Vector< real > mB;

            // alpha functions
            Cell< AlphaFunction * > mAlpha;

            //! scratch: cache tags and values of eval_a()
            mutable Vector< real > mCubicTemperatures;
            mutable Vector< real > mCubicStatevals;

            //! scratch for cardano
            mutable Vector< real > mWorkA;
            mutable Vector< real > mWorkZ;

            //! scratch: cache of the chi() departure term
            mutable real mDepartureV;
            mutable real mDepartureValue; // log( ( v - b*r2)/(v-b*r1))

            // splines for departure function at p = gPref
            Cell< Matrix<real> > mDepartureCoefficients;

            Spline mDepartureSpline;

            // component wise properties
            //! scratch: cache of update_component_parameters()
            mutable index_t mComponentCol; // column for departure spline
            mutable real mComponentT = BELFEM_QUIET_NAN;
            mutable real mComponentP = BELFEM_QUIET_NAN;
            mutable Vector< real > mComponentHDEP;
            mutable Vector< real > mComponentCPDEP;
            mutable Vector< real > mComponentV;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            EoS_Cubic( Gas & aParent, const GasModel aGasModel );

//----------------------------------------------------------------------------

            ~EoS_Cubic();

//----------------------------------------------------------------------------

            void
            remix();

//----------------------------------------------------------------------------
// Thermodynamic State
//----------------------------------------------------------------------------

            real
            p( const real T, const real v ) const;

            real
            v( const real T, const real p ) const;

            real
            T( const real p, const real v ) const;

//----------------------------------------------------------------------------
// State Derivatives
//----------------------------------------------------------------------------
           /**
             *
             * \f$ \frac{\partial p}{\partial T } = \frac{R}{ v - b} - \frac{a \,\frac{\partial \alpha}{\partial T }}{ \left( v-b \, r_1\right) \, \left( v-b \, r_2\right)} \f$
             */
            real
            dpdT( const real T, const real v ) const;

//----------------------------------------------------------------------------

            /**
              *
              * \f$ \frac{\partial^2 p}{\partial T^2} =  - \frac{a \, \frac{\partial^2 \alpha}{\partial T^2}}{ \left( v-b \, r_1\right) \, \left( v-b \, r_2\right)}
              * \f$
              */
            real
            d2pdT2( const real T, const real v ) const;

//----------------------------------------------------------------------------

            /**
              *
              * \f$ \frac{\partial p}{\partial v} = -\frac{R \, T }{\left( v - b\right)^2} - \frac{a \, \alpha \left[ b \, \left( r_1 + r_2 \right) - 2 \, v \right] }{ \left[\left( v-b \, r_1\right) \, \left( v-b \,r_2\right) \right]^2}
              * \f$
              */
            real
            dpdv( const real T, const real v ) const;

//----------------------------------------------------------------------------
            /**
              *
              * \f$ \frac{\partial^2 p}{\partial v^2} = 2 \, \left\{ \frac{R \, T }{\left( v - b\right)^3} - \frac{a \, \alpha \left[ b^2 \, \left( r_1^2 + r_1\,r_2 +  r_2^2 \right) - 3 \, b \, \left(r_1 + r_2 \right) \, v + 3 \, v^2 \right] }{ \left[\left( v-b \, r_1\right) \, \left( v-b \,r_2\right) \right]^3} \right\}
              * \f$
              */
            real
            d2pdv2( const real T, const real v ) const;

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            /**
             * thermal expansion coefficient
             *
             * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
             *
             * @param T temperature in K
             * @param p pressure in Pa
             *
             */
            real
            alpha( const real T, const real p ) const;

//------------------------------------------------------------------------------

            /**
             * isochoric stress coefficient
             *
             * \f$ \beta = \frac{1}{p} \left( \frac{\partial p}{\partial T}\right)_v \f$
             *
             * @param T temperature in K
             * @param p pressure in Pa
             *
             */
            real
            beta( const real T, const real p ) const;

//------------------------------------------------------------------------------

            /**
             * isothermal compressibility coefficient
             *
             * \f$ \kappa = -\frac{1}{v} \left( \frac{\partial v}{\partial p}\right)_T \f$
             *
             * @param T temperature in K
             * @param p pressure in Pa
             *
             */
            real
            kappa( const real T, const real p ) const;

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

            /**
             * enthalpy departure
            * \f$  h - h^{\circ} = p\, v - R\,T + \frac{T \, \frac{\partial a}{\partial T} - a}{b \, \left( r_2 - r_1\right)} \, \ln \left[ \frac{v - b \, r_1}{v - b \, r_2}\right]  \f$
            */
            real
            hdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            cpdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            sdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // pressure derivative of enthalpy departure
            real
            dhdepdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            real
            dsdepdT( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // pressure derivative of entropy departure
            real
            dsdepdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            hdep0( const real T ) const;

//------------------------------------------------------------------------------

            real
            cpdep0( const real T ) const;

//------------------------------------------------------------------------------

            real
            sdep0( const real T ) const;

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            real
            dsdepdT0( const real T ) const;

//------------------------------------------------------------------------------

            void
            eval_critical_point( real & T, real & p, real & v ) const;

//------------------------------------------------------------------------------
// Component Wise Specific Volume and Departure Functions
//------------------------------------------------------------------------------

            real
            v( const uint aIndex, const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            hdep( const uint aIndex, const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            cpdep( const uint aIndex, const real T, const real p ) const;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
//      INITIALIZATION
//------------------------------------------------------------------------------

            void
            init_srk();

//------------------------------------------------------------------------------

            void
            init_pr();

//------------------------------------------------------------------------------

            void
            init_common();

//------------------------------------------------------------------------------

            void
            init_departure_splines();

//------------------------------------------------------------------------------
//      A-FACTOR
//------------------------------------------------------------------------------

            /**
             * calculates the a-function, only to be used by
             * a, dadT and d2adT2
             * @param T       temperature in K
             * @param aDeriv   derivative 0, 1 or 2
             */
            void
            eval_a( const real T, const int aDeriv ) const;

//------------------------------------------------------------------------------

            /**
             * the "a" function for cubic gas
             */
            real
            a( const real T ) const;

//------------------------------------------------------------------------------

            /**
             * first derivative of the "a" function for cubic gas
             */
            real
            dadT( const real T ) const;

//------------------------------------------------------------------------------

            /**
             * second derivative of the "a" function for cubic gas
             */
            real
            d2adT2( const real T ) const;

//------------------------------------------------------------------------------

            /**
             * for departure function, returns
             * \f$ \ln \left( \frac{v - b \, r_1}{v - b \, r_2}\right)
             *     / \left( b \, ( r_2 - r_1 ) \right) \f$
             */
             real
             chi( const real v ) const;


             real
             dchidT( const real T, const real p, const real v ) const;

             real
             dchidp( const real T, const real p, const real v ) const;

//------------------------------------------------------------------------------
// component wise density etc
//------------------------------------------------------------------------------

            void
            component_wise_parameters(
                    const uint aIndex,
                    const real T,
                    const real p,
                    real & v,
                    real & aHDEP,
                    real & aCPDEP ) const;

//------------------------------------------------------------------------------

            void
            update_component_parameters(
                    const real T,
                    const real p ) const;

        };

//----------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */
#endif //BELFEM_CL_GM_EOS_CUBIC_HPP
