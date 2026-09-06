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

#include "constants.hpp"
#include "cl_GT_HeatPoly.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        HeatPoly::HeatPoly(
                const real             aTmin,
                const real             aTmax,
                const real             aEnthalpyConstant,
                const real             aEntropyConstant,
                const Vector< real > & aCoefficients ) :
                    mTmin( aTmin ),
                    mTmax( aTmax ),
                    mEnthalpyConstant( aEnthalpyConstant ),
                    mEntropyConstant( aEntropyConstant ),
                    mCoefficients( aCoefficients )
        {

        }

//------------------------------------------------------------------------------

        real
        HeatPoly::Cp( const real T ) const
        {
            return    (  (     mCoefficients( 0 )/T
                      +        mCoefficients( 1 ) )/T
                      +        mCoefficients( 2 )
                      + T * ( mCoefficients( 3 )
                      + T * ( mCoefficients( 4 )
                      + T * ( mCoefficients( 5 )
                      + T *   mCoefficients( 6 ) ))))*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::H( const real T ) const
        {

            return   ( mEnthalpyConstant
                       -        mCoefficients( 0 )/T
                       +        mCoefficients( 1 ) * std::log( T )
                       + T * ( mCoefficients( 2 )
                       + T * ( mCoefficients( 3 )*0.5
                       + T * ( mCoefficients( 4 )/3.0
                       + T * ( mCoefficients( 5 )*0.25
                       + T *   mCoefficients( 6 )*0.2 ) ) ) ) )*constant::Rm;
        }

//------------------------------------------------------------------------------


        real
        HeatPoly::S( const real T ) const
        {
            return   (        mEntropyConstant
                     -      ( mCoefficients( 0 ) * 0.5 / T
                     +        mCoefficients( 1 ) ) / T
                     +        mCoefficients( 2 ) * std::log( T )
                     + T * ( mCoefficients( 3 )
                     + T * ( mCoefficients( 4 ) * 0.5
                     + T * ( mCoefficients( 5 ) / 3.0
                     + T *   mCoefficients( 6 ) * 0.25 ) ) ) ) * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::dSdT( const real T ) const
        {
            return (        ( mCoefficients( 0 )
                     + T * ( mCoefficients( 1 )
                     + T *   mCoefficients( 2 ) ) ) / ( T*T*T )
                     +        mCoefficients( 3 )
                     + T * ( mCoefficients( 4 )
                     + T * ( mCoefficients( 5 )
                     + T *   mCoefficients( 6 ) ) ) ) * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::dCpdT( const real T ) const
        {
            return  ((- 2.0 *  mCoefficients( 0 ) / T
                      -        mCoefficients( 1 ) ) / ( T*T )
                      +        mCoefficients( 3 )
                      + T * ( mCoefficients( 4 ) * 2.0
                      + T * ( mCoefficients( 5 ) * 3.0
                      + T * ( mCoefficients( 6 ) * 4.0 ) ) ) ) * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::d2CpdT2( const real T ) const
        {
            return 2.0 * ( (
                              3.0 *mCoefficients( 0 )
                                 + mCoefficients( 1 )*T )/( T*T*T*T )
                                 + mCoefficients( 4 )
                    + T * ( 3.0 * mCoefficients( 5 )
                    + T *   6.0 * mCoefficients( 6 ) ) ) * constant::Rm ;
        }

//------------------------------------------------------------------------------

        void
        HeatPoly::set_T_min( const real aTmin )
        {
            mTmin = aTmin;
        }

//------------------------------------------------------------------------------

        void
        HeatPoly::set_T_max( const real aTmax )
        {
            mTmax = aTmax;
        }

//------------------------------------------------------------------------------

        void
        HeatPoly::set_enthalpy_constant( const real & aEnthalpyConstant )
        {
            mEnthalpyConstant = aEnthalpyConstant / constant::Rm;
        }

//------------------------------------------------------------------------------

        void
        HeatPoly::set_entropy_constant( const real & aEntropyConstant )
        {
            mEntropyConstant = aEntropyConstant / constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::enthalpy_constant() const
        {
            return mEnthalpyConstant * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::entropy_constant() const
        {
            return mEntropyConstant * constant::Rm;
        }

//------------------------------------------------------------------------------
    }
}