//
// Created by Christian Messe on 10.09.19.
//

#ifndef BELFEM_CL_GM_EOS_IDGAS_HPP
#define BELFEM_CL_GM_EOS_IDGAS_HPP

#include "typedefs.hpp"
#include "cl_GM_EoS.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * this is the ideal gas equaiton of state
         *
         * \f$ p \, v = R \, T \f$
         */
        class EoS_Idgas : public EoS
        {
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            EoS_Idgas( Gas & aParent );

//----------------------------------------------------------------------------

            ~EoS_Idgas() = default;

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

            real
            dpdT( const real aT, const real aV );

//----------------------------------------------------------------------------

            real
            d2pdT2( const real aT, const real aV );

//----------------------------------------------------------------------------

            real
            dpdv( const real aT, const real aV );

//----------------------------------------------------------------------------

            real
            d2pdv2( const real aT, const real aV );

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            real
            alpha( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            beta( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            kappa( const real aT, const real aP );

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

            real
            hdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            cpdep( const real aT, const real aP );

//------------------------------------------------------------------------------

            real
            sdep( const real aT, const real aP );

//------------------------------------------------------------------------------

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
//------------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */
#endif //BELFEM_CL_GM_EOS_IDGAS_HPP
