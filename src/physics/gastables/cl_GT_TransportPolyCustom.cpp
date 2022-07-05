//
// Created by Christian Messe on 26.08.19.
//


#include "cl_GT_TransportPolyCustom.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        TransportPolyCustom::TransportPolyCustom(
                const enum TransportPolyType aType,
                const real aTmin,
                const real aTmax,
                const Vector <real>  & aCoefficients,
                const Vector< real > & aExponents ) :
                TransportPoly( aType, aTmin, aTmax, aCoefficients ),
                mExponents( aExponents ),
                mNumberOfCoeffs( aExponents.length() )
        {
            BELFEM_ASSERT( aCoefficients.length() == aExponents.length(),
                          "Length of polynomial coeffients and exponents does not match" );

        }

//------------------------------------------------------------------------------

        real
        TransportPolyCustom::rawpoly( const real & aT ) const
        {
            real aValue = 0.0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 )
                {
                    aValue += mCoefficients( k ) * std::pow( aT, mExponents( k ));
                }
                else
                {
                    aValue += mCoefficients( k );
                }
            }
            return aValue;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyCustom::drawpoly( const real & aT ) const
        {
            real aValue = 0.0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 )
                {
                    if ( mExponents( k ) == 1.0 )
                    {
                        aValue += mCoefficients( k );
                    }
                    else
                    {
                        aValue += mExponents( k )
                                  * mCoefficients( k ) * std::pow( aT, mExponents( k ) - 1.0 );
                    }
                }
            }

            return aValue;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyCustom::ddrawpoly( const real & aT ) const
        {
            real aValue = 0;
            for ( uint k = 0; k < mNumberOfCoeffs; ++k )
            {
                if ( mExponents( k ) != 0 && mExponents( k ) != 1 )
                {
                    if ( mExponents( k ) == 2.0 )
                    {
                        aValue += 2.0 * mCoefficients( k );
                    }
                    else
                    {
                        aValue += ( mExponents( k ) - 1.0 ) * mExponents( k )
                                  * mCoefficients( k )
                                  * std::pow( aT, mExponents( k ) - 2.0 );
                    }
                }
            }

            return aValue;
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */