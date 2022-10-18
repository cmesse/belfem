//
// Created by christian on 10/17/22.
//
#include "constants.hpp"
#include "cl_BhExtrapolationOde.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "assert.hpp"

namespace belfem
{
//-----------------------------------------------------------------------------

    BhExtrapolationOde::BhExtrapolationOde( const real aBref, const real adHdBref,const real ad2HdB2ref, const real aBmax ) :
        ODE( 1 )
    {
        // check sanity of data
        BELFEM_ERROR( adHdBref < constant::nu0, "Error at last point of bh-curve must be < 1/mu0!" );

        // assume first value for switch value
        real tBold = 0 ;

        mBswitch = aBmax ;

        uint tCount = 0 ;
        real tLogNu = std::log( constant::nu0 );

        while( std::abs( tBold - mBswitch ) > BELFEM_EPSILON )
        {
            // compute the coefficients
            create_beam_poly( aBref,
                              std::log( adHdBref ),
                              ad2HdB2ref / adHdBref,
                              mBswitch,
                              tLogNu,
                              0, mCoeffs );

            // grab the coefficients for the first derivative
            real tA = 3 * mCoeffs( 0 );
            real tB = 2 * mCoeffs( 1 );
            real tC = mCoeffs( 2 );
            real tD = std::sqrt( tB * tB - 4.0 * tA * tC );

            // solve the quadratic equation
            real tX0 = ( -tB - tD ) / ( 2.0 * tA );
            real tX1 = ( -tB + tD ) / ( 2.0 * tA );

            tBold = mBswitch;
            mBswitch = tX0 < tX1 ? tX0 : tX1;

            BELFEM_ERROR( tCount++ < 100, "too many iterations" );
        }

        BELFEM_ERROR( mBswitch > aBref , "error while extrapolating b-h function" );
    }

//-----------------------------------------------------------------------------

    void
    BhExtrapolationOde::compute(
            const real           & aT,
            const Vector< real > & aY,
            Vector< real > & adYdT )
    {
        adYdT( 0 ) =
                aT < mBswitch ?
                std::exp( polyval( mCoeffs, aT ) ) : constant::nu0 ;
    }

//-----------------------------------------------------------------------------

}