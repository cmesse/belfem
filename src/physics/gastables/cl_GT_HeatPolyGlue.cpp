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
        HeatPolyGlue::Cp( const real T ) const
        {
        	return  ( ( ( (
                             mCoefficients( 0 )   * T
                          +  mCoefficients( 1 ) ) * T
                          +  mCoefficients( 2 ) ) * T
                          +  mCoefficients( 3 ) ) * T
                          +  mCoefficients( 4 ) )
                          * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::H( const real T ) const
        {
			return  ( ( ( ( ( (
                             mCoefficients( 0 )*0.2  ) * T
                          +  mCoefficients( 1 )*0.25 ) * T
                          +  mCoefficients( 2 )/3.0  ) * T
                          +  mCoefficients( 3 )*0.5  ) * T
                          +  mCoefficients( 4 )      ) * T + mEnthalpyConstant )
                          *   constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::S( const real T ) const
        {
        	return  ( ( ( ( ( (
                               mCoefficients( 0 )*0.25  ) * T
                            +  mCoefficients( 1 )/3.0  ) * T
                            +  mCoefficients( 2 )*0.5  ) * T
                            +  mCoefficients( 3 )  ) * T
                            +  mCoefficients( 4 ) * std::log( T ) )
                            + mEntropyConstant )
                        *   constant::Rm ;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::dSdT( const real T ) const
        {
            return (((   mCoefficients( 0 )   * T
                       + mCoefficients( 1 ) ) * T
                       + mCoefficients( 2 ) ) * T
                       + mCoefficients( 3 )
                       + mCoefficients( 4 ) / T ) * constant::Rm;
        }
//------------------------------------------------------------------------------

        real
        HeatPolyGlue::dCpdT( const real T ) const
        {
			return ( ( ( (
                          mCoefficients( 0 )* 4.0 ) * T
                       +  mCoefficients( 1 )* 3.0 ) * T
                       +  mCoefficients( 2 )* 2.0 ) * T
                       +  mCoefficients( 3 ) )
                       * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPolyGlue::d2CpdT2( const real T ) const
        {
			return ( ( (
                         mCoefficients( 0 )* 12.0 ) * T
                      +  mCoefficients( 1 )*  6.0 ) * T
                      +  mCoefficients( 2 )*  2.0 )
                       * constant::Rm;
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */