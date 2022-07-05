//
// Created by Christian Messe on 2019-03-10.
//

#include "cl_Gas.hpp"
#include "cl_GM_EoS_Idgas.hpp"


namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        EoS_Idgas::EoS_Idgas( Gas & aParent ) :
                EoS( aParent )
        {

        }

//----------------------------------------------------------------------------

        void
        EoS_Idgas::remix()
        {
            /* do nothing */
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::p( const real & aT, const real & aV )
        {
            return mR * aT / aV;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::v( const real & aT, const real & aP )
        {
            return mR * aT / aP;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::T( const real & aP, const real & aV )
        {
            return aP * aV / mR;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::dpdT( const real & aT, const real & aV )
        {
            return mR / aV;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::d2pdT2( const real & aT, const real & aV )
        {
            return 0.0;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::dpdv( const real & aT, const real & aV )
        {
            return -(aT*mR)/std::pow( aV, 2 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::d2pdv2( const real & aT, const real & aV )
        {
            return 2.0*aT*mR/std::pow( aV, 3 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::alpha( const real & aT, const real & aP )
        {
            mStatevals.update_Tp( aT, aP );

            if( ! mStatevals.test( BELFEM_STATEVAL_KAPPA ) )
            {
                // 10.18419/opus-9381 ( 2.8 )
                mStatevals.set( BELFEM_STATEVAL_ALPHA,
                                1.0 / aT );
            }

            return mStatevals.get( BELFEM_STATEVAL_ALPHA );
        }
//------------------------------------------------------------------------------

        real
        EoS_Idgas::beta( const real & aT, const real & aP )
        {
            mStatevals.update_Tp( aT, aP );

            if( ! mStatevals.test( BELFEM_STATEVAL_BETA ) )
            {
                // 10.18419/opus-9381 ( 2.6 )
                mStatevals.set( BELFEM_STATEVAL_BETA,
                                1.0 / aT );
            }

            return mStatevals.get( BELFEM_STATEVAL_BETA );
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::kappa( const real & aT, const real & aP )
        {
            mStatevals.update_Tp( aT, aP );

            if( ! mStatevals.test( BELFEM_STATEVAL_KAPPA ) )
            {

                // 10.18419/opus-9381 ( 2.7 )
                mStatevals.set( BELFEM_STATEVAL_KAPPA,
                                1.0 / aP );
            }

            return mStatevals.get( BELFEM_STATEVAL_KAPPA );
        }

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

        real
        EoS_Idgas::hdep( const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::cpdep( const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::sdep( const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::dsdepdT( const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::hdep0( const real & aT )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::cpdep0( const real & aT )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::dhdepdp( const real & aT, const real & aP )
        {
            return 0.0;
        }


//------------------------------------------------------------------------------

        real
        EoS_Idgas::sdep0( const real & aT )
        {
            return 0.0;
        }


//------------------------------------------------------------------------------

        real
        EoS_Idgas::dsdepdT0( const real & aT )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::dsdepdp( const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        void
        EoS_Idgas::eval_critical_point( real & aT, real & aP, real & aV )
        {
            // an ideal gas has no critical point!
            aT = BELFEM_QUIET_NAN;
            aP = BELFEM_QUIET_NAN;
            aV = BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::v( const uint & aIndex, const real & aT, const real & aP )
        {
            return mParent.data( aIndex )->R() * aT / aP ;
        }

//------------------------------------------------------------------------------
        real
        EoS_Idgas::hdep( const uint & aIndex, const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::cpdep( const uint & aIndex, const real & aT, const real & aP )
        {
            return 0.0;
        }

//------------------------------------------------------------------------------
    }
}