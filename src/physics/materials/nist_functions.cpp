//
// Created by christian on 4/21/22.
//

#include "nist_functions.hpp"

#include "assert.h"
#include "fn_create_beam_poly.hpp"

namespace belfem
{
    namespace nist
    {
//----------------------------------------------------------------------------

        real
        cp_janaf( const Vector< real > & aC, const real aT )
        {
            return ( aC( 0 ) + aC( 1 ) * aT ) / ( aT * aT )
                         + aC( 2 )
                         + aT * ( aC( 3 )
                         + aT * ( aC( 4 )
                         + aT * ( aC( 5 )
                         + aT *   aC( 6 ) ) ) ) ;

        }


//----------------------------------------------------------------------------

        real
        dcpdT_janaf( const Vector< real > & aC, const real aT )
        {
            return aC( 3 )
                   - ( 2.0 * aC( 0 ) + aC( 1 ) * aT ) / ( aT * aT * aT )
                   + aT * ( 2.0 * aC( 4 )
                   + aT * ( 3.0 * aC( 5 )
                   + aT * ( 4.0 * aC( 6 ) ) ) ) ;
        }


//----------------------------------------------------------------------------

        real
        proppoly( const Vector< real > & aC, const real aT  )
        {
            int tNumCoeffs = aC.length() ;

            BELFEM_ASSERT( aT != 0.0, "do not pass aT=0.0 into nist_poly!" );

            BELFEM_ASSERT( tNumCoeffs > 2,
                           "need at least three coefficients for NIST polynomial" );

            real tX = std::log10( aT );
            real tY = aC( tNumCoeffs-1 );


            for( int k=tNumCoeffs-2; k>=0; k-- )
            {
                tY *= tX ;
                tY += aC( k );
            }

            return std::pow( 10., tY );
        }

//----------------------------------------------------------------------------

        real
        dproppoly( const Vector< real > & aC, const real aT )
        {
            int tNumCoeffs = aC.length() ;

            BELFEM_ASSERT( aT != 0.0, "do not pass aT=0.0 into nist_poly!" );

            BELFEM_ASSERT( tNumCoeffs > 2,
                           "need at leat three coefficients for NIST polynomial" );

            real tX = std::log10( aT );
            real tdY = ( tNumCoeffs-1.) * aC( tNumCoeffs-1 );
            real tY = aC( tNumCoeffs-1 );

            for( int k=tNumCoeffs-2; k>=1; k-- )
            {
                tdY *= tX ;
                tdY += k * aC( k );
                tY  *= tX ;
                tY  += aC( k ) ;
            }
            tY *= tX ;
            tY += aC( 0 );

            return std::pow( 10., tY ) * tdY / aT ;
        }

//----------------------------------------------------------------------------

        real
        rho0( const Vector< real > & aP, const real aRRR, const real aT )
        {
            real tRho0 = aP( 0 ) / aRRR ;
            real tRhoi = aP( 1 ) * std::pow( aT, aP( 2 ) ) /
                         ( 1. + aP( 1 )*aP( 3 ) * std::pow( aT, aP( 2 ) - aP( 4 ) ) *
                                std::exp( -std::pow( aP( 5 ) / aT, aP( 6 ) ) ) ) ;

            return tRho0 + tRhoi + aP( 7 ) * tRhoi * tRho0 / ( tRhoi + tRho0 );
        }

//----------------------------------------------------------------------------

        real
        lambda0( const Vector< real > & aP, const real aBeta, const real aT )
        {
            // help constants
            real a = aP( 1 )*std::pow(aT, aP( 2 ) ) ;
            real c = aP( 1 )*aP( 3 )*std::pow(aT,aP( 2 )+aP( 4 )) ;
            real d = std::exp( -std::pow(aP( 5 )/aT,aP( 6 ) ) ) ;
            real b = 1. + c*d ;

            real w0  = aBeta / aT ;
            real wi  = a/b;
            real wi0 = aP( 7 ) * ( wi*w0/(wi+w0) );

            return 1./( w0 + wi + wi0 );
        }

//----------------------------------------------------------------------------

        real
        dinvlambda0dT( const Vector< real > & aP, const real aBeta, const real aT )
        {
            // help constants
            real a = aP( 1 )*std::pow(aT, aP( 2 ) ) ;
            real c = aP( 1 )*aP( 3 )*std::pow(aT,aP( 2 )+aP( 4 )) ;
            real d = std::exp( -std::pow(aP( 5 )/aT,aP( 6 ) ) ) ;
            real b = 1. + c*d ;

            // derivatives of help constants
            real da = a * aP( 2 )/aT ;
            real dc = c * (aP( 2 )+aP( 4 ))/aT ;
            real dd = d * aP( 6 ) * std::pow(aP( 5 )/aT,aP( 6 ) )/ aT ;
            real db = c*dd + dc*d ;

            real w0 = aBeta / aT ;
            real dw0 = -w0/ aT ;
            real wi = a/b;
            real dwi = (b*da - a*db)/(b*b);
            //real wi0 = aP7 * ( wi*w0/(wi+w0) );
            real dwi0 = aP( 7 )/((wi+w0)*(wi+w0))*(wi*wi*dw0 + w0*w0*dwi);

            return dw0 + dwi + dwi0 ;
        }

//----------------------------------------------------------------------------

        real
        extend_kohler( const Vector< real >  & aA,
                    const real                 aKmin,
                    const real                 aXmin0,
                    Vector< real >           & aB )
        {
            real tK = BELFEM_REAL_MAX ;
            real tdKdX ;
            uint c =0 ;

            real aXmin = aXmin0 ;

            while( std::abs( tK ) > BELFEM_EPS )
            {
                tK =  std::pow( 10, polyval( aA, std::log10( aXmin ) ) ) ;
                tdKdX = tK * dpolyval( aA, std::log10( aXmin ) ) / aXmin ;
                tK -= aKmin ;
                aXmin -= tK / tdKdX ;
                BELFEM_ERROR( c++ < 100, "too many iterations" );
            }

            tK =  std::pow( 10, polyval( aA, std::log10( aXmin ) ) ) ;
            tdKdX = tK * dpolyval( aA, std::log10( aXmin ) ) / aXmin ;

            aB.set_size( 4 );
            create_beam_poly( 0, 0.0, 0.0, aXmin, tK, tdKdX, aB );

            return aXmin ;
        }

//-----------------------------------------------------------------------

        void
        extend_resistivity( const Vector< real > & aPoly1,
                            const   real     aTswitch,
                            Vector< real > & aPoly0 )
        {
            real tF  = polyval( aPoly1, aTswitch );
            real tdF = dpolyval( aPoly1, aTswitch );

            aPoly0.set_size( 3, 0 );
            aPoly0( 0 ) = tF - tdF * aTswitch * 0.5 ;
            aPoly0( 2 ) = tdF /( 2.0 * aTswitch );
        }

//-----------------------------------------------------------------------
    }
}