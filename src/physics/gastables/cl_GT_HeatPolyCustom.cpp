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

#include "assert.hpp"
#include "constants.hpp"
#include "cl_GT_HeatPolyCustom.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        HeatPolyCustom::HeatPolyCustom(
                const real             aTmin,
                const real             aTmax,
                const real             aEnthalpyConstant,
                const real             aEntropyConstant,
                const Vector< real > & aCoefficients,
                const Vector< real > & aExponents ) :
                    HeatPoly(
                            aTmin,
                            aTmax,
                            aEnthalpyConstant,
                            aEntropyConstant,
                            aCoefficients ),
                    mExponents( aExponents ),
                    mNumberOfCoeffs( aExponents.length() )
                {
                    BELFEM_ASSERT( aCoefficients.length() == aExponents.length(),
                            "Length of polynomial coeffients and exponents does not match" );
                }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::Cp( const real T ) const
        {
            real cp = 0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0.0 )
                {
                    cp += mCoefficients( k ) * std::pow( T, mExponents( k ));
                }
                else
                {
                    cp += mCoefficients( k );
                }
            }
            return cp*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::H( const real T ) const
        {
            real h = mEnthalpyConstant;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) == -1.0 )
                {
                    h +=  mCoefficients( k ) * std::log( T );
                }
                else
                {
                    real tE =  mExponents( k ) + 1.0;

                    h += mCoefficients( k ) * std::pow( T, tE ) / tE;
                }
            }
            return h*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::S( const real T ) const
        {
            real s = mEntropyConstant;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) == 0.0 )
                {
                    s +=  mCoefficients( k ) * std::log( T );
                }
                else
                {
                    s += mCoefficients( k ) * std::pow( T, mExponents( k ) )
                            / mExponents( k );
                }
            }
            return s*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::dSdT( const real T ) const
        {
            real adSdT = 0.0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) == 1.0 )
                {
                    adSdT += mCoefficients( k );
                }
                else
                {
                    adSdT += mCoefficients( k ) * std::pow( T, mExponents( k ) - 1.0 );
                }
            }

            return adSdT*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::dCpdT( const real T ) const
        {
            real cp = 0.0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 )
                {
                    if( mExponents( k ) == 1.0 )
                    {
                        cp += mCoefficients( k );
                    }
                    else
                    {
                        cp += mExponents( k ) * mCoefficients( k )
                               * std::pow( T, mExponents( k ) - 1.0 );
                    }
                }
            }

            return cp*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::d2CpdT2( const real T ) const
        {
            real cp = 0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 && mExponents( k ) != 1 )
                {
                    if( mExponents( k ) == 2.0 )
                    {
                        cp += 2.0 * mCoefficients( k );
                    }
                    else
                    {
                        cp += ( mExponents( k ) - 1.0 ) *mExponents( k )
                                * mCoefficients( k )
                                * std::pow( T, mExponents( k ) - 2.0 );
                    }
                }
            }

            return cp*constant::Rm;
        }

//------------------------------------------------------------------------------
    }
}