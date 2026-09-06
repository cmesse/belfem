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

#include <gastables/GT_globals.hpp>
#include "cl_GM_EoS_AlphaFunction.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        AlphaFunction::AlphaFunction( ) :
                mTcrit( 0.0 ),
                mInvTcrit( 0.0 ),
                mC1( 0.0  ),
                mC2( 0.0  ),
                mC3( 0.0  )
        {

        }
//----------------------------------------------------------------------------

        AlphaFunction::AlphaFunction( const real aTcrit,
                                      const Vector< real > & aCoeffs ) :
                mTcrit( aTcrit ),
                mInvTcrit( 1.0/aTcrit ),
                mC1( aCoeffs(0 ) ),
                mC2( aCoeffs( 1) ),
                mC3( aCoeffs( 2) )
        {

        }


//----------------------------------------------------------------------------

        AlphaFunction::AlphaFunction(
                const real aTcrit,
                const real & c1,
                const real & c2,
                const real & c3 ) :
                    mTcrit( aTcrit ),
                    mInvTcrit( 1.0/aTcrit ),
                    mC1( c1 ),
                    mC2( c2 ),
                    mC3( c3 )
        {

        }


//----------------------------------------------------------------------------

        real
        AlphaFunction::alpha( const real T ) const
        {
            BELFEM_ASSERT( T > 0, "Temperature must be positive for alpha function" );
            this->eval( T, 0 );
            return mAlpha[ 0 ];
        }

//----------------------------------------------------------------------------

        real
        AlphaFunction::dalphadT( const real T ) const
        {
            BELFEM_ASSERT( T > 0, "Temperature must be positive for alpha function" );
            this->eval( T, 1 );
            return mAlpha[ 1 ];
        }

//----------------------------------------------------------------------------

        real
        AlphaFunction::d2alphadT2( const real T ) const
        {
            BELFEM_ASSERT( T > 0, "Temperature must be positive for alpha function" );
            this->eval( T, 2 );
            return mAlpha[ 2 ];
        }

//----------------------------------------------------------------------------

        void
        AlphaFunction::eval( const real T, const int aDeriv ) const
        {
           BELFEM_ERROR( false, "this function must not be called");
        }

//----------------------------------------------------------------------------

        void
        AlphaFunction::find_minimum()
        {
            mTmin = gastables::gTmax + 1.0;

            real tT0 = 10.0;
            real tT1 = 6000;

            uint tCount = 0;

            real tF0 = this->dalphadT( tT0 );
            real tF = this->dalphadT( tT1 );

            real tT = tT1 ;

            // first run with bisection
                while( std::abs( tF ) > 0.0001 )
                {
                    tT = 0.5 * ( tT0 + tT1 );
                    tF = this->dalphadT( tT );

                    if ( tF0 * tF  > 0 )
                    {
                        tT0 = tT;
                        tF0 = tF;
                    }
                    else
                    {
                        tT1 = tT;
                    }

                    ++tCount;

                BELFEM_ERROR( tCount < 1000, "Created infinite loop while trying to find mimumum of alpha function" );

            }

            tCount = 0;

            // final run with newton
            tT0 = -1.0;

            while( std::abs( tT - tT0 ) > 1e-12 )
            {
                tT0 = tT;
                tT -= this->dalphadT( tT ) / this->d2alphadT2( tT );

                ++tCount;

                BELFEM_ERROR( tCount < 1000, "Created infinite loop while trying to find mimumum of alpha function" );
            }

            mTmin = tT;
            mAlphaMin = this->alpha( tT );
        }

//----------------------------------------------------------------------------

        AlphaFunction_Empty::AlphaFunction_Empty() :
                AlphaFunction()
        {
        }


//----------------------------------------------------------------------------
        real
        AlphaFunction_Empty::alpha( const real T ) const
        {
            return 0.0;
        }

//----------------------------------------------------------------------------

        real
        AlphaFunction_Empty::dalphadT( const real T ) const
        {
            return 0.0;
        }

//----------------------------------------------------------------------------

        real
        AlphaFunction_Empty::d2alphadT2( const real T ) const
        {
            return 0.0;
        }

//----------------------------------------------------------------------------

        AlphaFunction_Classic::AlphaFunction_Classic(
                const real aTcrit,
                const real & c1,
                const real & c2 ) :
                AlphaFunction( aTcrit,
                               c1,
                               c2,
                               c1*c2)
        {
        }

//----------------------------------------------------------------------------

        void
        AlphaFunction_Classic::eval( const real T, const int aDeriv ) const
        {
            // calculate function itself
            if ( mT[ 0 ] != T )
            {
                mWork[ 0 ] = std::sqrt( T );

                // remember this temperature
                mT[ 0 ] = T;

                mAlpha[ 0 ] = std::pow( 1.0 + mC1 * ( 1.0 - mC2 * mWork[ 0 ] ), 2 );

            }

            // calculate first derivative
            if ( aDeriv >= 1 )
            {
                // first derivative
                if ( mT[ 1 ] != T )
                {
                    // remember this temperature
                    mT[ 1 ] = T;

                    mAlpha[ 1 ] = ( mC3 * ( mC1 * ( mWork[ 0 ] * mC2 - 1.0 ) - 1.0 )) / mWork[ 0 ];
                }

                if ( aDeriv >= 2 )
                {
                    // second derivative
                    if ( mT[ 2 ] != T )
                    {
                        // remember this temperature
                        mT[ 2 ] = T;

                        mAlpha[ 2 ] = 0.5 * ( mC3 * mC3 - mAlpha[ 1 ] ) / T;
                    }
                }
            }
        }

//----------------------------------------------------------------------------

        AlphaFunction_MC::AlphaFunction_MC(
                const real aTcrit,
                const real & c1,
                const real & c2,
                const real & c3):
                AlphaFunction( aTcrit,
                               c1,
                               c2,
                               c3)
        {
            this->find_minimum();
        }

//----------------------------------------------------------------------------

        void
        AlphaFunction_MC::eval( const real T, const int aDeriv ) const
        {

            if( mT[ 0 ] != T )
            {
                // remember this temperature
                mT[ 0 ] = T;
                mX[ 0 ] = std::sqrt( T*mInvTcrit );
                mX[ 1 ] = 1.0 - mX[ 0 ];

                if( T < mTmin )
                {
                    if ( T < mTcrit )
                    {
                        mF[ 0 ] = 1.0 + (( mC3 * mX[ 1 ] + mC2 ) * mX[ 1 ] + mC1 ) * mX[ 1 ];
                    }
                    else
                    {
                        mF[ 0 ] = 1.0 + mC1 * mX[ 1 ];
                    }

                    mAlpha[ 0 ] = mF[ 0 ] * mF[ 0 ];
                }
                else
                {
                    mAlpha[ 0 ] = mAlphaMin;
                }
            }

            if( aDeriv >= 1 )
            {
                if( mT[ 1 ] != T )
                {
                    // remember this temperature
                    mT[ 1 ] = T;
                    mX[ 2 ] = - 1.0 / ( 2.0 * mTcrit * mX[ 0 ] );
                    if( T < mTmin )
                    {
                        if ( T < mTcrit )
                        {
                            mF[ 1 ] = ( mC1 + mX[ 1 ] * ( 2.0 * mC2 + 3.0 * mC3 * mX[ 1 ] )) * mX[ 2 ];
                        }
                        else
                        {
                            mF[ 1 ] = mC1 * mX[ 2 ];
                        }

                        mAlpha[ 1 ] = 2.0 * mF[ 0 ] * mF[ 1 ];
                    }
                    else
                    {
                        mAlpha[ 1 ] = 0.0;
                    }
                }

                if( aDeriv >= 2 )
                {
                    // second derivative
                    if ( mT[ 2 ] != T )
                    {
                        // remember this temperature
                        mT[ 2 ] = T;
                        mX[ 3 ] = 0.25 / ( T * mTcrit * mX[ 0 ] );
                        if( T < mTmin )
                        {
                            if ( T < mTcrit )
                            {
                                mF[ 2 ] = 2.0 * ( mC2 + 3.0 * mC3 * mX[ 1 ] ) * mX[ 2 ] * mX[ 2 ]
                                          + ( mC1 + mX[ 1 ] * ( 2.0 * mC2 + 3.0 * mC3 * mX[ 1 ] )) * mX[ 3 ];
                            }
                            else
                            {
                                mF[ 2 ] = mC1 * mX[ 3 ];
                            }

                            mAlpha[ 2 ] = 2.0 * ( mF[ 1 ] * mF[ 1 ] + mF[ 0 ] * mF[ 2 ] );
                        }
                        else
                        {
                            mAlpha[ 2 ] = 0.0;
                        }
                    }
                }
            }
        }

//----------------------------------------------------------------------------


        AlphaFunction_CCR::AlphaFunction_CCR(
                const real aTcrit,
                const real & c1,
                const real & c2,
                const real & c3):
            AlphaFunction( aTcrit,
                           c1,
                           c2,
                           c3)
        {
        }

//----------------------------------------------------------------------------

        void
        AlphaFunction_CCR::eval( const real T, const int aDeriv ) const
        {

            if( mT[ 0 ] != T )
            {
                // remember this temperature
                mT[ 0 ] = T;

                mF[ 0 ] = mC1 * ( 1.0 - T*mInvTcrit );

                if( T < mTcrit )
                {
                    mWork[ 0 ] = std::sqrt( T*mInvTcrit );
                    mWork[ 1 ] = 1.0 - mWork[ 0 ];

                    mG[ 0 ] = 1.0 + ( mC2 + mC3 * mWork[ 1 ] ) \
                        * mWork[ 1 ] * mWork[ 1 ];

                    mAlpha[ 0 ] = std::exp( mF[ 0 ] ) * mG [ 0 ] * mG[ 0 ];

                }
                else
                {
                    mAlpha[ 0 ] = std::exp( mF[ 0 ] );
                }
            }

            if( aDeriv >= 1 )
            {
                if( mT[ 1 ] != T )
                {
                    // remember this temperature
                    mT[ 1 ] = T;

                    mF[ 1 ] = -mC1 * mInvTcrit;

                    if ( T < mTcrit )
                    {
                        // d( 1 - sqrt( T/T_crit ) )/dT
                        mWork[ 2 ] = -0.5 * mInvTcrit / mWork[ 0 ];

                        mG[ 1 ] = mWork[ 1 ]
                                  * ( 2.0 * mC2 + 3.0 * mC3 * mWork[ 1 ] )
                                  * mWork[ 2 ];

                        mAlpha[ 1 ] = mG[ 0 ] * std::exp( mF[ 0 ] ) *
                                      ( mG[ 0 ] * mF[ 1 ] + 2.0 * mG[ 1 ] );
                    }
                    else
                    {
                        mAlpha[ 1 ] = mF[ 1 ] * mAlpha[ 0 ];
                    }
                }

                if( aDeriv >= 2 )
                {
                    // second derivative
                    if ( mT[ 2 ] != T )
                    {
                        // remember this temperature
                        mT[ 2 ] = T;

                        if ( T < mTcrit )
                        {
                            // d2( 1 - sqrt( T/T_crit ) )/dT2
                            mWork[ 3 ] = 0.25 / ( T * mTcrit * mWork[ 0 ] ); // ok

                            mG[ 2 ] = 2.0 * ( mC2 + 3.0 * mC3 * mWork[ 1 ] )
                                      * mWork[ 2 ] * mWork[ 2 ]
                                      + mWork[ 1 ] * ( 2.0 * mC2 + 3.0 * mC3 *mWork[ 1 ] )
                                      * mWork[ 3 ]; // ok


                            mAlpha[ 2 ] = std::exp( mF[ 0 ] ) *
                                    ( mG[0]*mG[0]*mF[1]*mF[1]
                                     + 4.0 * mG[0] * mF[1]*mG[1]
                                     + 2.0 * ( mG[1]*mG[1] + mG[0]*mG[2])); // ok

                        }
                        else
                        {
                            mAlpha[ 2 ] = mF[ 1 ] * mAlpha[ 1 ];
                        }
                    }
                }
            }
        }

//----------------------------------------------------------------------------

        AlphaFunction_PM::AlphaFunction_PM(
                    const real             aTcrit,
                    const Vector< real > & aCoeffs ):
                AlphaFunction( aTcrit, aCoeffs)
        {
        }
//----------------------------------------------------------------------------

        void
        AlphaFunction_PM::eval( const real T, const int aDeriv ) const
        {
            // calculate function itself
            if( mT[ 0 ] != T )
            {
                // remember this temperature
                mT[ 0 ] = T;

                mX      = std::sqrt( T * mInvTcrit );

                mY[ 0 ] = 1.0 - mX;

                mF[ 0 ] = mY[ 0 ] * ( mC1 + mY[ 0 ] * ( mC2 + mY[ 0 ]* mC3 ) );

                mAlpha[ 0 ] = std::exp( mF[ 0 ] );
            }

            // calculate first derivative
            if( aDeriv >= 1 )
            {
                if( mT[ 1 ] != T )
                {
                    // remember this temperature
                    mT[ 1 ] = T;

                    // derivative of help functions
                    mY[ 1 ] = -0.5 * mInvTcrit / mX;

                    mF[ 1 ] = ( mC1 + mY[ 0 ] * ( 2.0 * mC2 + 3.0*mC3 * mY[ 0 ] ) ) * mY[ 1 ];

                    mAlpha[ 1 ] = mAlpha[ 0 ] * mF[ 1 ];
                }

                if( aDeriv == 2 )
                {
                    if( mT[ 2 ] != T )
                    {
                        // remember this temperature
                        mT[ 2 ] = T;

                        // derivative of help functions
                        mY[ 2 ] = - 0.5 * mY[ 1 ] / T ;

                        mF[ 2 ] =  ( 2.0 * mC2 + 6.0 * mC3 * mY[ 0 ] ) * std::pow( mY[ 1 ], 2 )
                                   + ( mC1 + mY[ 0 ] * ( 2.0 * mC2 + 3.0 * mC3 * mY[ 0 ] ) ) * mY[ 2 ];

                        mAlpha[ 2 ] = mAlpha[ 1 ] * mF[ 1 ]
                                      +mAlpha[ 0 ] *  mF[ 2 ];
                    }
                }
            }
        }

//----------------------------------------------------------------------------
    }
}