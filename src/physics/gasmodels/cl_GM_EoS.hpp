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

#ifndef BELFEM_CL_GM_EOS_HPP
#define BELFEM_CL_GM_EOS_HPP

#include "typedefs.hpp"

namespace belfem
{
    // forward declaration for parent
    class Gas;

    namespace gasmodels
    {
        // forward declaration for statevals
        class Statevals;

        /**
         * virtual class for equation of state
         *
         * @ingroup grp_physics_gasmodels
         * @see @ref physics_gasmodels_gasmodels_usage_guide
         */
        class EoS
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * gas object that owns this class
             */
            Gas & mParent;

            /**
             * statevals object
             */
            Statevals & mStatevals;

             /**
              * specific gas constant
              */
            const real & mR;

            /**
             * molar mass in kg/Mol
             */
            const real & mM;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            EoS( Gas & aParent );

//------------------------------------------------------------------------------

            virtual ~EoS() = default;

//------------------------------------------------------------------------------

            /**
             * call the remix function of the equation of state
             */
            virtual void
            remix() ;

//------------------------------------------------------------------------------
// Thermodynamic States
//------------------------------------------------------------------------------

            /**
             * pressure in Pa
             */
            virtual real
            p( const real T, const real v ) const = 0;

//------------------------------------------------------------------------------

            /**
             * specific volume in m^3/kg
             */
            virtual real
            v( const real T, const real p ) const = 0;

//------------------------------------------------------------------------------

            /**
             * temperature in K
             */
            virtual real
            T( const real p, const real v ) const = 0;

//------------------------------------------------------------------------------
// State Derivatives
//------------------------------------------------------------------------------

            virtual real
            dpdT( const real T, const real v ) const = 0;

//------------------------------------------------------------------------------

            virtual real
            d2pdT2( const real T, const real v ) const;

//------------------------------------------------------------------------------

            virtual real
            dpdv( const real T, const real v ) const = 0;

//------------------------------------------------------------------------------

            virtual real
            d2pdv2( const real T, const real v ) const;

//------------------------------------------------------------------------------

            virtual real
            dvdT( const real T, const real v ) const;

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            /**
             * thermal expansion coefficient
             *
             * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
             */
            virtual real
            alpha( const real T, const real p ) const;

//------------------------------------------------------------------------------


            /**
             * isochoric stress coefficient
             *
             * \f$ \beta = \frac{1}{p} \left( \frac{\partial p}{\partial T}\right)_v \f$
             */
            virtual real
            beta( const real T, const real p ) const;

//------------------------------------------------------------------------------

            /**
             * isothermal compressibility coefficient
             *
             * \f$ \kappa = -\frac{1}{v} \left( \frac{\partial v}{\partial p}\right)_T \f$
             */
            virtual real
            kappa( const real T, const real p ) const;

//------------------------------------------------------------------------------
// Caloric Functions ( only for special fluid models )
//------------------------------------------------------------------------------

            virtual real
            h( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            hvap( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            cv( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            cp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            s( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            dsdT( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            dsdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            /**
             * speed of sound in m/s ( helmholtz only )
             * @param T
             * @param p
             * @return
             */
            virtual real
            w( const real T, const real p ) const;

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

            virtual real
            hdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            cpdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            sdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            dhdepdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            virtual real
            dsdepdT( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // pressure derivative of entropy departure
            virtual real
            dsdepdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            virtual real
            hdep0( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            cpdep0( const real T ) const;

//------------------------------------------------------------------------------

            virtual real
            sdep0( const real T ) const;

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            virtual real
            dsdepdT0( const real T ) const;

//------------------------------------------------------------------------------

            virtual void
            eval_critical_point( real & T, real & p, real & v ) const = 0;

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

            virtual real
            v( const uint aIndex, const real T, const real p ) const ;

            virtual real
            hdep( const uint aIndex, const real T, const real p ) const ;

            virtual real
            cpdep( const uint aIndex, const real T, const real p ) const ;

//------------------------------------------------------------------------------

            // special function for tablegas
            virtual real
            pi( const real T, const real p ) const;

            // special function for tablegas
            virtual real
            dpidp( const real T, const real p ) const;

//------------------------------------------------------------------------------
// vaporization curve
//----------------------------------------------------------------------------

            /**
             * @param T  vapor temperature in K
             * @return    vapor pressure in Pa
             */
            virtual real
            p_vap( const real T ) const;

//----------------------------------------------------------------------------

            /**
             * @param p  vapor pressure in Pa
             * @return    vapor temperature in K
             */
            virtual real
            T_vap( const real p ) const;

//------------------------------------------------------------------------------

            /**
             * expose the parent of this EoS
             */
             Gas *
             parent();

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------

        inline Gas *
        EoS::parent()
        {
            return & mParent ;
        }

//----------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */
#endif //BELFEM_CL_GM_EOS_HPP
