//
// Created by Christian Messe on 25.08.19.
//

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
                const Vector <real> & aCoefficients ) :
                mType( aType ),
                mScale( ( aType == TransportPolyType::VISCOSITY ) ? 1e-7 : 1e-4 ),
                mTmin( aTmin ),
                mTmax( aTmax ),
                mCoefficients( aCoefficients )
        {

        }

//------------------------------------------------------------------------------

        real
        TransportPoly::rawpoly( const belfem::real & aT ) const
        {
            return mCoefficients( 0 ) * std::log( aT ) +
                   ( mCoefficients( 1 ) + mCoefficients( 2 )/aT )/aT
                   + mCoefficients( 3 );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::drawpoly( const belfem::real & aT ) const
        {
            return ( (  mCoefficients( 0 )*aT - mCoefficients( 1 ))*aT
                    - 2.0*mCoefficients( 2 ))* std::pow( aT, -3 );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::ddrawpoly( const belfem::real & aT ) const
        {
            return  ( 6.0*mCoefficients( 2 ) + aT*( 2.0*mCoefficients( 1 )
                    - aT*mCoefficients( 0 ) ) )  * std::pow( aT, -4 );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::eval( const real aT ) const
        {
            return std::exp( this->rawpoly( aT ) ) * mScale;
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::deval( const real aT ) const
        {
            return this->eval( aT ) * this->drawpoly( aT );
        }

//------------------------------------------------------------------------------

        real
        TransportPoly::ddeval( const real aT ) const
        {
            return this->eval( aT ) * (
                    std::pow( this->drawpoly( aT ) , 2 )
                    + this->ddrawpoly( aT ) );
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