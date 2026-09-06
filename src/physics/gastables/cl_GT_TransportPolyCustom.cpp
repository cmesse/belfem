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
                TransportPoly( aType, aTmin, aTmax, aCoefficients, TransportPolyKind::CUSTOM ),
                mExponents( aExponents ),
                mNumberOfCoeffs( aExponents.length() )
        {
            BELFEM_ASSERT( aCoefficients.length() == aExponents.length(),
                          "Length of polynomial coeffients and exponents does not match" );

        }

//------------------------------------------------------------------------------

        real
        TransportPolyCustom::rawpoly( const real T ) const
        {
            real aValue = 0.0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 )
                {
                    aValue += mCoefficients( k ) * std::pow( T, mExponents( k ));
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
        TransportPolyCustom::drawpoly( const real T ) const
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
                                  * mCoefficients( k ) * std::pow( T, mExponents( k ) - 1.0 );
                    }
                }
            }

            return aValue;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyCustom::ddrawpoly( const real T ) const
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
                                  * std::pow( T, mExponents( k ) - 2.0 );
                    }
                }
            }

            return aValue;
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */