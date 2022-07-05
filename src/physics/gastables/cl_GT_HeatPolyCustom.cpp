//
// Created by Christian Messe on 2019-08-19.
//

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
        HeatPolyCustom::Cp( const real & aT ) const
        {
            real aCp = 0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0.0 )
                {
                    aCp += mCoefficients( k ) * std::pow( aT, mExponents( k ));
                }
                else
                {
                    aCp += mCoefficients( k );
                }
            }
            return aCp*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::H( const real & aT ) const
        {
            real aH = mEnthalpyConstant;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) == -1.0 )
                {
                    aH +=  mCoefficients( k ) * std::log( aT );
                }
                else
                {
                    real tE =  mExponents( k ) + 1.0;

                    aH += mCoefficients( k ) * std::pow( aT, tE ) / tE;
                }
            }
            return aH*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::S( const real & aT ) const
        {
            real aS = mEntropyConstant;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) == 0.0 )
                {
                    aS +=  mCoefficients( k ) * std::log( aT );
                }
                else
                {
                    aS += mCoefficients( k ) * std::pow( aT, mExponents( k ) )
                            / mExponents( k );
                }
            }
            return aS*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::dSdT( const real & aT ) const
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
                    adSdT += mCoefficients( k ) * std::pow( aT, mExponents( k ) - 1.0 );
                }
            }

            return adSdT*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::dCpdT( const real & aT ) const
        {
            real aCp = 0.0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 )
                {
                    if( mExponents( k ) == 1.0 )
                    {
                        aCp += mCoefficients( k );
                    }
                    else
                    {
                        aCp += mExponents( k ) * mCoefficients( k )
                               * std::pow( aT, mExponents( k ) - 1.0 );
                    }
                }
            }

            return aCp*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyCustom::d2CpdT2( const real & aT ) const
        {
            real aCp = 0;
            for( uint k=0; k<mNumberOfCoeffs; ++k )
            {
                if( mExponents( k ) != 0 && mExponents( k ) != 1 )
                {
                    if( mExponents( k ) == 2.0 )
                    {
                        aCp += 2.0 * mCoefficients( k );
                    }
                    else
                    {
                        aCp += ( mExponents( k ) - 1.0 ) *mExponents( k )
                                * mCoefficients( k )
                                * std::pow( aT, mExponents( k ) - 2.0 );
                    }
                }
            }

            return aCp*constant::Rm;
        }

//------------------------------------------------------------------------------
    }
}