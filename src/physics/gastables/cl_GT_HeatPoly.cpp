//
// Created by Christian Messe on 2019-08-19.
//

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
        HeatPoly::Cp( const real & aT ) const
        {
            return    (  (     mCoefficients( 0 )/aT
                      +        mCoefficients( 1 ) )/aT
                      +        mCoefficients( 2 )
                      + aT * ( mCoefficients( 3 )
                      + aT * ( mCoefficients( 4 )
                      + aT * ( mCoefficients( 5 )
                      + aT *   mCoefficients( 6 ) ))))*constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::H( const real & aT ) const
        {

            return   ( mEnthalpyConstant
                       -        mCoefficients( 0 )/aT
                       +        mCoefficients( 1 ) * std::log( aT )
                       + aT * ( mCoefficients( 2 )
                       + aT * ( mCoefficients( 3 )*0.5
                       + aT * ( mCoefficients( 4 )/3.0
                       + aT * ( mCoefficients( 5 )*0.25
                       + aT *   mCoefficients( 6 )*0.2 ) ) ) ) )*constant::Rm;
        }

//------------------------------------------------------------------------------


        real
        HeatPoly::S( const real & aT ) const
        {
            return   (        mEntropyConstant
                     -      ( mCoefficients( 0 ) * 0.5 / aT
                     +        mCoefficients( 1 ) ) / aT
                     +        mCoefficients( 2 ) * std::log( aT )
                     + aT * ( mCoefficients( 3 )
                     + aT * ( mCoefficients( 4 ) * 0.5
                     + aT * ( mCoefficients( 5 ) / 3.0
                     + aT *   mCoefficients( 6 ) * 0.25 ) ) ) ) * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::dSdT( const real & aT ) const
        {
            return          ( mCoefficients( 0 )
                     + aT * ( mCoefficients( 1 )
                     + aT *   mCoefficients( 2 ) ) ) / ( aT*aT*aT )
                     +        mCoefficients( 3 )
                     + aT * ( mCoefficients( 4 )
                     + aT * ( mCoefficients( 5 )
                     + aT *   mCoefficients( 6 ) ) ) * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::dCpdT( const real & aT ) const
        {
            return  ((- 2.0 *  mCoefficients( 0 ) / aT
                      -        mCoefficients( 1 ) ) / ( aT*aT )
                      +        mCoefficients( 3 )
                      + aT * ( mCoefficients( 4 ) * 2.0
                      + aT * ( mCoefficients( 5 ) * 3.0
                      + aT * ( mCoefficients( 6 ) * 4.0 ) ) ) ) * constant::Rm;
        }

//------------------------------------------------------------------------------

        real
        HeatPoly::d2CpdT2( const real & aT ) const
        {
            return 2.0 * ( (
                              3.0 *mCoefficients( 0 )
                                 + mCoefficients( 1 )*aT )*std::pow( aT, -4.0 )
                                 + mCoefficients( 4 )
                    + aT * ( 3.0 * mCoefficients( 5 )
                    + aT *   6.0 * mCoefficients( 6 ) ) ) * constant::Rm ;
        }

//------------------------------------------------------------------------------

        void
        HeatPoly::set_T_min( const real & aTmin )
        {
            mTmin = aTmin;
        }

//------------------------------------------------------------------------------

        void
        HeatPoly::set_T_max( const real & aTmax )
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