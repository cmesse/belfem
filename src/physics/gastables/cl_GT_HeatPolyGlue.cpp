//
// Created by Christian Messe on 2019-08-22.
//

#include "assert.hpp"
#include "constants.hpp"
#include "cl_GT_HeatPolyGlue.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        HeatPolyGlue::HeatPolyGlue(
                const real             aTmin,
                const real             aTmax,
                const real             aEnthalpyConstant,
                const real             aEntropyConstant,
                const Vector< real > & aCoefficients ) :
                    HeatPoly(
                            aTmin,
                            aTmax,
                            aEnthalpyConstant,
                            aEntropyConstant,
                            aCoefficients )
        {

        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::Cp( const real & aT ) const
        {
        	return  ( ( ( (
                             mCoefficients( 0 )   * aT
                          +  mCoefficients( 1 ) ) * aT
                          +  mCoefficients( 2 ) ) * aT
                          +  mCoefficients( 3 ) ) * aT
                          +  mCoefficients( 4 ) )
                          * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::H( const real & aT ) const
        {
			return  ( ( ( ( ( (
                             mCoefficients( 0 )*0.2  ) * aT
                          +  mCoefficients( 1 )*0.25 ) * aT
                          +  mCoefficients( 2 )/3.0  ) * aT
                          +  mCoefficients( 3 )*0.5  ) * aT
                          +  mCoefficients( 4 )      ) * aT + mEnthalpyConstant )
                          *   constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::S( const real & aT ) const
        {
        	return  ( ( ( ( ( (
                               mCoefficients( 0 )*0.25  ) * aT
                            +  mCoefficients( 1 )/3.0  ) * aT
                            +  mCoefficients( 2 )*0.5  ) * aT
                            +  mCoefficients( 3 )  ) * aT
                            +  mCoefficients( 4 ) * std::log( aT ) )
                            + mEntropyConstant )
                        *   constant::Rm ;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::dSdT( const real & aT ) const
        {
            return (((   mCoefficients( 0 )   * aT
                       + mCoefficients( 1 ) ) * aT
                       + mCoefficients( 2 ) ) * aT
                       + mCoefficients( 3 )
                       + mCoefficients( 4 ) / aT ) * constant::Rm;
        }
//------------------------------------------------------------------------------

        real
        HeatPolyGlue::dCpdT( const real & aT ) const
        {
			return ( ( ( (
                          mCoefficients( 0 )* 4.0 ) * aT
                       +  mCoefficients( 1 )* 3.0 ) * aT
                       +  mCoefficients( 2 )* 2.0 ) * aT
                       +  mCoefficients( 3 ) )
                       * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::d2CpdT2( const real & aT ) const
        {
			return ( ( (
                         mCoefficients( 0 )* 12.0 ) * aT
                      +  mCoefficients( 1 )*  6.0 ) * aT
                      +  mCoefficients( 2 )*  2.0 )
                       * constant::Rm;
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */