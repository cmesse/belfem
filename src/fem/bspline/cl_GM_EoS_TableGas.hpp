//
// Created by Christian Messe on 01.05.20.
//

#ifndef BELFEM_CL_GM_EOS_TABLEGAS_HPP
#define BELFEM_CL_GM_EOS_TABLEGAS_HPP

#include "typedefs.hpp"
#include "cl_GM_EoS.hpp"
#include "cl_TableGas.hpp"
#include "cl_BS_LookupTable.hpp"

namespace belfem
{
    namespace gasmodels
    {
        /**
         * this is the ideal gas equaiton of state for a frozen table gas
         *
         * \f$ p \, v = R \, T \f$
         */
        class EoS_TableGas : public EoS
        {
            bspline::LookupTable & mTable ;

            const uint mIndexM ;

            // mimumum pressure on table
            const real mPimin ;

            // maximum pressure on table
            const real mPimax ;

            // mimumum pressure on table
            const real mPmin ;

            // maximum pressure on table
            const real mPmax ;

            // constant for scaling derivative to p
            const real mdpscale = 1000.0 / std::log( 10.0 );


//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            EoS_TableGas( TableGas & aParent );

//----------------------------------------------------------------------------

            ~EoS_TableGas() = default;

//----------------------------------------------------------------------------

            void
            remix();

//----------------------------------------------------------------------------
// Special
//----------------------------------------------------------------------------

            real
            pi( const real aT, const real aP );

            real
            dpidp(const real &aT, const real &aP);

            const real &
            M( const real aT, const real aP );

            real
            dMdT( const real aT, const real aP );

            real
            dMdp( const real aT, const real aP );

//----------------------------------------------------------------------------
// Thermodynamic State
//----------------------------------------------------------------------------

            real
            p( const real aT, const real aV );

            real
            v( const real aT, const real aP );

            real
            T( const real aP, const real aV );

//------------------------------------------------------------------------------
// State Derivatives
//------------------------------------------------------------------------------

            real
            dpdT( const real aT, const real aV );

//------------------------------------------------------------------------------

            real
            d2pdT2( const real aT, const real aV );

//------------------------------------------------------------------------------

            real
            dpdv( const real aT, const real aV );

//------------------------------------------------------------------------------

            real
            d2pdv2( const real aT, const real aV );

//------------------------------------------------------------------------------

            real
            dvdT( const real aT, const real aV );

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            /**
             * thermal expansion coefficient
             *
             * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
             */
            real
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
// Departure Functions ( all zero for this class )
//------------------------------------------------------------------------------

            inline real hdep( const real aT, const real aP ) { return 0.0; };
            inline real cpdep( const real aT, const real aP ) { return 0.0; };
            inline real sdep( const real aT, const real aP )  { return 0.0; };
            inline real dhdepdp( const real aT, const real aP ) { return 0.0; };
            inline real dsdepdT( const real aT, const real aP ) { return 0.0; };
            inline real dsdepdp( const real aT, const real aP )  { return 0.0; };
            inline real hdep0( const real aT ) { return 0.0; };
            inline real cpdep0( const real aT ) { return 0.0; };
            inline real sdep0( const real aT ) { return 0.0; };
            inline real dsdepdT0( const real aT ) { return 0.0; };

//------------------------------------------------------------------------------
// Forbidden funcitons
//------------------------------------------------------------------------------

            void
            eval_critical_point( real & aT, real & aP, real & aV );

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
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
        };

    }
}
#endif //BELFEM_CL_GM_EOS_TABLEGAS_HPP
