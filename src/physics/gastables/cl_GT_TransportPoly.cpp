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

#include "cl_GT_TransportPoly.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        TransportPoly::TransportPoly(
            const enum TransportPolyType aType,
            const real aTmin,
            const real aTmax,
            const Vector <real> & aCoefficients,
            const TransportPolyKind aKind ):
                mType( aType ),
                mKind( aKind ),
                mScale( ( aType == TransportPolyType::VISCOSITY ) ? 1e-7 : 1e-4 ),
                mTmin( aTmin ),
                mTmax( aTmax ),
                mCoefficients( aCoefficients )
        {

        }

//------------------------------------------------------------------------------

        real
        TransportPoly::rawpoly( const belfem::real T ) const
        {
            return mCoefficients( 0 ) * std::log( T ) +
                   ( mCoefficients( 1 ) + mCoefficients( 2 )/T )/T
                   + mCoefficients( 3 );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::drawpoly( const belfem::real T ) const
        {
            return ( (  mCoefficients( 0 )*T - mCoefficients( 1 ))*T
                    - 2.0*mCoefficients( 2 ))* std::pow( T, -3 );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::ddrawpoly( const belfem::real T ) const
        {
            return  ( 6.0*mCoefficients( 2 ) + T*( 2.0*mCoefficients( 1 )
                    - T*mCoefficients( 0 ) ) )  * std::pow( T, -4 );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::eval( const real T ) const
        {
            return std::exp( this->rawpoly( T ) ) * mScale;
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::deval( const real T ) const
        {
            return this->eval( T ) * this->drawpoly( T );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::ddeval( const real T ) const
        {
            return this->eval( T ) * (
                    std::pow( this->drawpoly( T ) , 2 )
                    + this->ddrawpoly( T ) );
        }

//------------------------------------------------------------------------------

        void
        TransportPoly::set_T_min( const real aTmin )
        {
            mTmin = aTmin;
        }

//------------------------------------------------------------------------------

        void
        TransportPoly::set_T_max( const real aTmax )
        {
            mTmax = aTmax;
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */