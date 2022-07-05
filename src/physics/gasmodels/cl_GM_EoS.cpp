//
// Created by Christian Messe on 09.09.19.
//

#include "cl_GM_EoS.hpp"
#include "cl_GM_Statevals.hpp"
#include "cl_Gas.hpp"

namespace belfem
{
    namespace gasmodels
    {
//------------------------------------------------------------------------------

        EoS::EoS( Gas & aParent ) :
                mParent( aParent ),
                mStatevals( aParent.statevals() ),
                mR( aParent.statevals().get( BELFEM_STATEVAL_R ) ),
                mM( aParent.statevals().get( BELFEM_STATEVAL_M ) )
        {

        }

//------------------------------------------------------------------------------

        void
        EoS::remix()
        {
            BELFEM_ERROR( false, "EoS::remix not implemented for this EoS" );
        }

//------------------------------------------------------------------------------

        real
        EoS::alpha( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::alpha not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::beta( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::beta not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::kappa( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::kappa not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::dvdT( const real & aT, const real & aP )
        {
            return - this->dpdT( aT, aP ) / this->dpdv( aT, aP );
        }

//------------------------------------------------------------------------------

        real
        EoS::h( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::h not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::hvap( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::hvap not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::cv( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::cv not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::cp( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::cp not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::w( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::w not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::s( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::s not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::dsdT( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::dsdT not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::dsdp( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::dsdp not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }
        
        
//------------------------------------------------------------------------------

        real
        EoS::pi( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::pi not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::dpidp( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::dpidp not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::d2pdT2( const real & aT, const real & aV )
        {
            BELFEM_ERROR( false, "EoS::d2pdT2 not implemented for this EoS");
            return BELFEM_QUIET_NAN ;
        }


//------------------------------------------------------------------------------

        real
        EoS::d2pdv2( const real & aT, const real & aV )
        {
            BELFEM_ERROR( false, "EoS::d2pdv2 not implemented for this EoS");
            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

        real
        EoS::hdep( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::hdep not implemented for this EoS");
            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        real
        EoS::cpdep( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::cpdep not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::sdep( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::sdep not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::dhdepdp( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::dhdepdp not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        // temperature derivative of entropy departure
        real
        EoS::dsdepdT( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::dsdepdT not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        // temperature derivative of entropy departure
        real
        EoS::dsdepdp( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "EoS::dsdepdp not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::hdep0( const real & aT )
        {
            BELFEM_ERROR( false, "EoS::hdep0 not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::cpdep0( const real & aT )
        {
            BELFEM_ERROR( false, "EoS::cpdep0 not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS::sdep0( const real & aT )
        {
            BELFEM_ERROR( false, "EoS::sdep0 not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        // temperature derivative of entropy departure
        real
        EoS::dsdepdT0( const real & aT )
        {
            BELFEM_ERROR( false, "EoS::dsdepdT0 not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

        real
        EoS::v( const uint & aIndex, const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "component-wise EoS::v not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

        real
        EoS::hdep( const uint & aIndex, const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "component-wise EoS::hdep not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

        real
        EoS::cpdep( const uint & aIndex, const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "component-wise EoS::cpdep not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }
        
//------------------------------------------------------------------------------

        real
        EoS::p_vap( const real & aT )
        {
            BELFEM_ERROR( false, "p_vap not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }


        real
        EoS::T_vap( const real & aP )
        {
            BELFEM_ERROR( false, "T_vap not implemented for this EoS" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */