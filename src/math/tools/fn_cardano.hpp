//
// Created by Christian Messe on 13.09.19.
//

#ifndef BELFEM_FN_CARDANO_HPP
#define BELFEM_FN_CARDANO_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "fn_sort.hpp"
#include "fn_sign.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    /**
     * solve a cubic equation a*x^3 + b*x^2 + c*d + d = 0 to x
     */
    inline void
    cardano( const  Vector< real > & aA, Vector< real > & aX )
    {
        if ( std::abs( aA( 0 ) ) < BELFEM_EPSILON )
        {
            if( std::abs( aA( 1 ) ) < BELFEM_EPSILON )
            {
                if( std::abs( aA( 2 ) ) > BELFEM_EPSILON )
                {
                    aX.set_size( 1 );

                    aX = -aA( 3 ) / aA( 2 );
                }
                else
                {
                    aX.set_size( 0 );
                }
            }
            else
            {
                real tD = std::pow( aA( 2 ), 2 ) - 4.0 * aA( 1 ) * aA( 3 );

                if( tD >= 0.0 )
                {
                    tD = std::sqrt( tD );
                    aX.set_size( 2 );
                    aX( 0 ) = -aA( 2 ) - tD;
                    aX( 1 ) = -aA( 2 ) + tD;
                    aX /= 2.0 * aA( 1 );
                }
                else
                {
                    aX.set_size( 0 );
                }
            }
        }
        else
        {

            real tP = ( 9.0*aA(0)*aA(2)-3.0*std::pow( aA(1),2) )/(9.0*std::pow( aA(0), 2 ));
            real tQ = ( 2.0* std::pow( aA(1), 3) - 9.0*aA(0)*aA(1)*aA(2) + 27.0*std::pow(aA(0),2)*aA(3))/(27.0*std::pow( aA(0),3));
            real tR  = aA(1)/(3*aA(0));
            real tD = 0.25 * std::pow( tQ, 2 ) + std::pow( tP, 3)/27.0;

            if ( tD > 0.0 )
            {
                aX.set_size( 1 );
                real tU = -0.5*tQ + std::sqrt( tD );
                real tV = -0.5*tQ - std::sqrt( tD );
                tU = sign( tU )*std::pow( std::abs( tU ), 1.0/3.0 );
                tV = sign( tV )*std::pow( std::abs( tV ), 1.0/3.0 );

                cplx tX = tU + tV;



                if( std::abs( std::imag( tX ) ) < BELFEM_EPSILON )
                {
                    aX( 0 ) = std::real( tX ) - tR;
                }
                else
                {
                    const cplx tF1( -0.5, 0.5 * std::sqrt( 3.0 ));
                    const cplx tF2( -0.5, -0.5 * std::sqrt( 3.0 ));

                    tX = tF1 * tU + tF2 * tV;
                    if( std::abs( std::imag( tX ) ) < BELFEM_EPSILON )
                    {
                        aX( 0 ) = std::real( tX ) - tR;
                    }
                    else
                    {
                        tX = tF2*tU + tF1 * tV;
                        if( std::abs( std::imag( tX ) ) < BELFEM_EPSILON )
                        {
                            aX( 0 ) = std::real( tX ) - tR;
                        }
                        else
                        {
                            BELFEM_ERROR( false, "Something went wrong while trying to solve cubic equaition" );
                        }
                    }

                }

            }
            else if( tD < 0 )
            {
                cplx tA = std::sqrt( -4.0/3.0*tP );
                cplx tB = std::acos( -0.5*tQ*std::sqrt( -27.0/std::pow( tP, 3 ) ) )/3.0;
                const cplx tC = 2.0*std::acos( 0.0 )/3.0;

                aX.set_size( 3 );
                aX( 0 ) = std::real(  tA * std::cos( tB ) ) -tR;
                aX( 1 ) = std::real( -tA * std::cos( tB + tC ) ) -tR;
                aX( 2 ) = std::real( -tA * std::cos( tB - tC ) ) -tR;
                sort( aX );
            }
            else
            {
                aX.set_size( 1 );
                aX( 0 ) = -tR;
            }
        }

    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_CARDANO_HPP
