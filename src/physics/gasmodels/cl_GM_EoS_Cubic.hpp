//
// Created by Christian Messe on 13.09.19.
//

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

            // values of alpha
            Vector< real > mA;
            Vector< real > mdAdT;
            Vector< real > md2AdT2;

            // contains b, r1 and r2
            Vector< real > mB;

            // alpha functions
            Cell< AlphaFunction * > mAlpha;

            Vector< real > mCubicTemperatures;
            Vector< real > mCubicStatevals;

            // for cardano
            Vector< real > mWorkA;
            Vector< real > mWorkZ;

            // for departure
            real mDepartureV;
            real mDepartureValue; // log( ( v - b*r2)/(v-b*r1))

            // splines for departure function at p = gPref
            Cell< Matrix<real> > mDepartureCoefficients;

            Spline mDepartureSpline;

            // Vector< real > mVM;

            // special flag
            bool mUseAlwaysCardano = true;

            // component wise properties
            index_t mComponentCol; // column for departure spline
            real mComponentT = BELFEM_QUIET_NAN;
            real mComponentP = BELFEM_QUIET_NAN;
            Vector< real > mComponentHDEP;
            Vector< real > mComponentCPDEP;
            Vector< real > mComponentV;

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
            p( const real aT, const real aV );

            real
            v( const real aT, const real aP );

            real
            T( const real aP, const real aV );

//----------------------------------------------------------------------------
// State Derivatives
//----------------------------------------------------------------------------
           /**
             *
             * \f$ \frac{\partial p}{\partial T } = \frac{R}{ v - b} - \frac{a \,\frac{\partial \alpha}{\partial T }}{ \left( v-b \, r_1\right) \, \left( v-b \, r_2\right)} \f$
             */
            real
            dpdT( const real aT, const real aV );

//----------------------------------------------------------------------------

            /**
              *
              * \f$ \frac{\partial^2 p}{\partial T^2} =  - \frac{a \, \frac{\partial^2 \alpha}{\partial T^2}}{ \left( v-b \, r_1\right) \, \left( v-b \, r_2\right)}
\f$
              */
            real
            d2pdT2( const real aT, const real aV );

//----------------------------------------------------------------------------

            /**
              *
              * \f$ \frac{\partial p}{\partial v} = -\frac{R \, T }{\left( v - b\right)^2} - \frac{a \, \alpha \left[ b \, \left( r_1 + r_2 \right) - 2 \, v \right] }{ \left[\left( v-b \, r_1\right) \, \left( v-b \,r_2\right) \right]^2}
 \f$
              */
            real
            dpdv( const real aT, const real aV );

//----------------------------------------------------------------------------
            /**
              *
              * \f$ \frac{\partial^2 p}{\partial v^2} = 2 \, \left\{ \frac{R \, T }{\left( v - b\right)^3} - \frac{a \, \alpha \left[ b^2 \, \left( r_1^2 + r_1\,r_2 +  r_2^2 \right) - 3 \, b \, \left(r_1 + r_2 \right) \, v + 3 \, v^2 \right] }{ \left[\left( v-b \, r_1\right) \, \left( v-b \,r_2\right) \right]^3} \right\}
 \f$
              */
            real
            d2pdv2( const real aT, const real aV );

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            /**
             * thermal expansion coefficient
             *
             * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
             *
             * @param aT temperature in K
             * @param aP pressure in Pa
             *
             */
            real
            alpha( const real aT, const real aP );

//------------------------------------------------------------------------------

            /**
             * isochoric stress coefficient
             *
             * \f$ \beta = \frac{1}{p} \left( \frac{\partial p}{\partial T}\right)_v \f$
             *
             * @param aT temperature in K
             * @param aP pressure in Pa
             *
             */
            real
            beta( const real aT, const real aP );

//------------------------------------------------------------------------------

            /**
             * isothermal compressibility coefficient
             *
             * \f$ \kappa = -\frac{1}{v} \left( \frac{\partial v}{\partial p}\right)_T \f$
             *
             * @param aT temperature in K
             * @param aP pressure in Pa
             *
             */
            real
            kappa( const real aT, const real aP );

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

            /**
             * enthalpy departure
            * \f$  h - h^{\circ} = p\, v - R\,T + \frac{T \, \frac{\partial a}{\partial T} - a}{b \, \left( r_2 - r_1\right)} \, \ln \left[ \frac{v - b \, r_1}{v - b \, r_2}\right]  \f$
            */
            real
            hdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            cpdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            sdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            // pressure derivative of enthalpy departure
            real
            dhdepdp( const real aT, const real aP );

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            real
            dsdepdT( const real aT, const real aP );

//------------------------------------------------------------------------------

            // pressure derivative of entropy departure
            real
            dsdepdp( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            hdep0( const real aT );

//------------------------------------------------------------------------------

            real
            cpdep0( const real aT );

//------------------------------------------------------------------------------

            real
            sdep0( const real aT );

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            real
            dsdepdT0( const real aT );

//------------------------------------------------------------------------------

            void
            eval_critical_point( real & aT, real & aP, real & aV );

//------------------------------------------------------------------------------
// Component Wise Specofic Volume and Departure Functions
//------------------------------------------------------------------------------

            real
            v( const uint aIndex, const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            hdep( const uint aIndex, const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            cpdep( const uint aIndex, const real aT, const real aP );

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
             * @param aT       temperature in K
             * @param aDeriv   derivative 0, 1 or 2
             */
            void
            eval_a( const real aT, const int aDeriv );

//------------------------------------------------------------------------------

            /**
             * the "a" function for cubic gas
             */
            real
            a( const real aT );

//------------------------------------------------------------------------------

            /**
             * first derivative of the "a" function for cubic gas
             */
            real
            dadT( const real aT );

//------------------------------------------------------------------------------

            /**
             * second derivative of the "a" function for cubic gas
             */
            real
            d2adT2( const real aT );

//------------------------------------------------------------------------------

            /**
             * for departure function, returns
             * \f$ \ln \left( \frac{v - b \, r_2}{v - b \, r_1}\right) \f$
             */
             real
             chi( const real aV );


             real
             dchidT( const real aT, const real aP, const real aV );

             real
             dchidp( const real aT, const real aP, const real aV );

//------------------------------------------------------------------------------
// component wise density etc
//------------------------------------------------------------------------------

            void
            component_wise_parameters(
                    const uint aIndex,
                    const real aT,
                    const real aP,
                    real & aV,
                    real & aHDEP,
                    real & aCPDEP );

//------------------------------------------------------------------------------

            void
            update_component_parameters(
                    const real aT,
                    const real aP );

        };

//----------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */
#endif //BELFEM_CL_GM_EOS_CUBIC_HPP
