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
        EoS_Idgas::p( const real T, const real v ) const
        {
            return mR * T / v;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::v( const real T, const real p ) const
        {
            return mR * T / p;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::T( const real p, const real v ) const
        {
            return p * v / mR;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::dpdT( const real T, const real v ) const
        {
            return mR / v;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::d2pdT2( const real T, const real v ) const
        {
            return 0.0;
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::dpdv( const real T, const real v ) const
        {
            return -(T*mR)/std::pow( v, 2 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::d2pdv2( const real T, const real v ) const
        {
            return 2.0*T*mR/std::pow( v, 3 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Idgas::alpha( const real T, const real p ) const
        {
            mStatevals.update_Tp( T, p );

            if( ! mStatevals.test( BELFEM_STATEVAL_ALPHA ) )
            {
                // 10.18419/opus-9381 ( 2.8 )
                mStatevals.set( BELFEM_STATEVAL_ALPHA,
                                1.0 / T );
            }

            return mStatevals.get( BELFEM_STATEVAL_ALPHA );
        }
//------------------------------------------------------------------------------

        real
        EoS_Idgas::beta( const real T, const real p ) const
        {
            mStatevals.update_Tp( T, p );

            if( ! mStatevals.test( BELFEM_STATEVAL_BETA ) )
            {
                // 10.18419/opus-9381 ( 2.6 )
                mStatevals.set( BELFEM_STATEVAL_BETA,
                                1.0 / T );
            }

            return mStatevals.get( BELFEM_STATEVAL_BETA );
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::kappa( const real T, const real p ) const
        {
            mStatevals.update_Tp( T, p );

            if( ! mStatevals.test( BELFEM_STATEVAL_KAPPA ) )
            {

                // 10.18419/opus-9381 ( 2.7 )
                mStatevals.set( BELFEM_STATEVAL_KAPPA,
                                1.0 / p );
            }

            return mStatevals.get( BELFEM_STATEVAL_KAPPA );
        }

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

        real
        EoS_Idgas::hdep( const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::cpdep( const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::sdep( const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::dsdepdT( const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::hdep0( const real T ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::cpdep0( const real T ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::dhdepdp( const real T, const real p ) const
        {
            return 0.0;
        }


//------------------------------------------------------------------------------

        real
        EoS_Idgas::sdep0( const real T ) const
        {
            return 0.0;
        }


//------------------------------------------------------------------------------

        real
        EoS_Idgas::dsdepdT0( const real T ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::dsdepdp( const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        void
        EoS_Idgas::eval_critical_point( real & T, real & p, real & v ) const
        {
            // an ideal gas has no critical point!
            T = BELFEM_QUIET_NAN;
            p = BELFEM_QUIET_NAN;
            v = BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::v( const uint aIndex, const real T, const real p ) const
        {
            return mParent.data( aIndex )->R() * T / p ;
        }

//------------------------------------------------------------------------------
        real
        EoS_Idgas::hdep( const uint aIndex, const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        EoS_Idgas::cpdep( const uint aIndex, const real T, const real p ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------
    }
}