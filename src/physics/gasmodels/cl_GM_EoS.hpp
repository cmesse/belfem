//
// Created by Christian Messe on 09.09.19.
//

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
            p( const real aT, const real aV ) = 0;

//------------------------------------------------------------------------------

            /**
             * specific volume in m^3/kg
             */
            virtual real
            v( const real aT, const real aP ) = 0;

//------------------------------------------------------------------------------

            /**
             * temperature in K
             */
            virtual real
            T( const real aP, const real aV ) = 0;

//------------------------------------------------------------------------------
// State Derivatives
//------------------------------------------------------------------------------

            virtual real
            dpdT( const real aT, const real aV ) = 0;

//------------------------------------------------------------------------------

            virtual real
            d2pdT2( const real aT, const real aV );

//------------------------------------------------------------------------------

            virtual real
            dpdv( const real aT, const real aV ) = 0;

//------------------------------------------------------------------------------

            virtual real
            d2pdv2( const real aT, const real aV );

//------------------------------------------------------------------------------

            virtual real
            dvdT( const real aT, const real aV );

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            /**
             * thermal expansion coefficient
             *
             * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
             */
            virtual real
            alpha( const real aT, const real aP );

//------------------------------------------------------------------------------


            /**
             * isochoric stress coefficient
             *
             * \f$ \beta = \frac{1}{p} \left( \frac{\partial p}{\partial T}\right)_v \f$
             */
            virtual real
            beta( const real aT, const real aP );

//------------------------------------------------------------------------------

            /**
             * isothermal compressibility coefficient
             *
             * \f$ \kappa = -\frac{1}{v} \left( \frac{\partial v}{\partial p}\right)_T \f$
             */
            virtual real
            kappa( const real aT, const real aP );

//------------------------------------------------------------------------------
// Caloric Functions ( only for special fluid models )
//------------------------------------------------------------------------------

            virtual real
            h( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            hvap( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            cv( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            cp( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            s( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            dsdT( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            dsdp( const real aT, const real aP );

//------------------------------------------------------------------------------

            /**
             * speed of sound in m/s ( helmholtz only )
             * @param aT
             * @param aP
             * @return
             */
            virtual real
            w( const real aT, const real aP );

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

            virtual real
            hdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            cpdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            sdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            dhdepdp( const real aT, const real aP );

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            virtual real
            dsdepdT( const real aT, const real aP );

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            virtual real
            dsdepdp( const real aT, const real aP );

//------------------------------------------------------------------------------

            virtual real
            hdep0( const real aT );

//------------------------------------------------------------------------------

            virtual real
            cpdep0( const real aT );

//------------------------------------------------------------------------------

            virtual real
            sdep0( const real aT );

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            virtual real
            dsdepdT0( const real aT );

//------------------------------------------------------------------------------

            virtual void
            eval_critical_point( real & aT, real & aP, real & aV ) = 0;

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

            virtual real
            v( const uint aIndex, const real aT, const real aP ) ;

            virtual real
            hdep( const uint aIndex, const real aT, const real aP ) ;

            virtual real
            cpdep( const uint aIndex, const real aT, const real aP ) ;

//------------------------------------------------------------------------------

            // special function for tablegas
            virtual real
            pi( const real aT, const real aP );

            // special function for tablegas
            virtual real
            dpidp( const real aT, const real aP );

//------------------------------------------------------------------------------
// vaporization curve
//----------------------------------------------------------------------------

            /**
             * @param aT  vapor temperature in K
             * @return    vapor pressure in Pa
             */
            virtual real
            p_vap( const real aT );

//----------------------------------------------------------------------------

            /**
             * @param aP  vapor pressure in Pa
             * @return    vapor temperature in K
             */
            virtual real
            T_vap( const real aP );

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
