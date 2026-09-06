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

#include <cmath>

#include "assert.hpp"
#include "constants.hpp"
#include "cl_Gas.hpp"
#include "cl_GM_HelmholtzTransport_LemmonJacobsen.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        HelmholtzTransport_LemmonJacobsen::HelmholtzTransport_LemmonJacobsen(
                Gas & aParent,
                const HelmholtzModel aModel ) :
            HelmholtzTransport( aParent ),
            mModel( aModel )
        {
            this->init_tables() ;

            /* Eq. ( 2 ) prefactor. M is in g/mol and sigma in nm, the units the
             * paper writes the equation in, and the result is micro Pa s */
            mConstEta0 = 0.0266958 * std::sqrt( mM ) / ( mSigma * mSigma ) ;

            // twice the critical temperature, Table I footnote
            mTref = 2.0 * mTcrit ;
        }

//----------------------------------------------------------------------------

        void
        HelmholtzTransport_LemmonJacobsen::init_tables()
        {
            // Table II, common to every fluid of the correlation
            mB = { 0.431, -0.4623, 0.08406, 0.005341, -0.00331 } ;

            switch( mModel )
            {
                case( HelmholtzModel::Nitrogen ) :
                {
                    // Table I
                    mTcrit     = 126.192 ;
                    mPcrit     = 3.3958e6 ;
                    mM         = 28.01348 ;
                    mRhocrit   = 11.1839e3 * mM * 1.0e-3 ;   // mol/dm^3 -> kg/m^3
                    mEpsilonKb = 98.94 ;
                    mSigma     = 0.3656 ;

                    // Table III
                    mEtaN = {  10.72,  0.03989, 0.001208, -7.402, 4.620 } ;
                    mEtaT = {   0.1,   0.25,    3.2,        0.9,  0.3   } ;
                    mEtaD = {   2.0,  10.0,    12.0,        2.0,  1.0   } ;
                    mEtaL = {   0.0,   1.0,     1.0,        2.0,  3.0   } ;

                    // Table IV, rows 1 to 3
                    mLambdaN1 =  1.511 ;
                    mLambdaN2 =  2.117 ;  mLambdaT2 = -1.0 ;
                    mLambdaN3 = -3.332 ;  mLambdaT3 = -0.7 ;

                    // Table IV, from row 4
                    mLambdaN = { 8.862, 31.11, -73.13, 20.03, -0.7096, 0.2672 } ;
                    mLambdaT = { 0.0,    0.03,   0.2,   0.8,   0.6,    1.9    } ;
                    mLambdaD = { 1.0,    2.0,    3.0,   4.0,   8.0,   10.0    } ;
                    mLambdaL = { 0.0,    0.0,    1.0,   2.0,   2.0,    2.0    } ;

                    // Table I, critical enhancement
                    mXi0    = 0.17e-9 ;
                    mGammaC = 0.055 ;
                    mQd     = 0.40e-9 ;

                    break ;
                }
                case( HelmholtzModel::Oxygen ) :
                {
                    // Table I
                    mTcrit     = 154.581 ;
                    mPcrit     = 5.043e6 ;
                    mM         = 31.9988 ;
                    mRhocrit   = 13.63e3 * mM * 1.0e-3 ;     // mol/dm^3 -> kg/m^3
                    mEpsilonKb = 118.5 ;
                    mSigma     = 0.3428 ;

                    // Table III
                    mEtaN = { 17.67, 0.4042, 0.0001077, 0.3510, -13.67 } ;
                    mEtaT = {  0.05, 0.0,    2.10,      0.0,      0.5  } ;
                    mEtaD = {  1.0,  5.0,   12.0,       8.0,      1.0  } ;
                    mEtaL = {  0.0,  0.0,    0.0,       1.0,      2.0  } ;

                    // Table IV, rows 1 to 3
                    mLambdaN1 =  1.036 ;
                    mLambdaN2 =  6.283 ;  mLambdaT2 = -0.9 ;
                    mLambdaN3 = -4.262 ;  mLambdaT3 = -0.6 ;

                    // Table IV, from row 4
                    mLambdaN = { 15.31, 8.898, -0.7336, 6.728, -4.374, -0.4747 } ;
                    mLambdaT = {  0.0,  0.0,    0.3,    4.3,    0.5,    1.8    } ;
                    mLambdaD = {  1.0,  3.0,    4.0,    5.0,    7.0,   10.0    } ;
                    mLambdaL = {  0.0,  0.0,    0.0,    2.0,    2.0,    2.0    } ;

                    // Table I, critical enhancement
                    mXi0    = 0.24e-9 ;
                    mGammaC = 0.055 ;
                    mQd     = 0.51e-9 ;

                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false,
                        "Lemmon-Jacobsen transport is implemented for nitrogen "
                        "and oxygen only" );
                    break ;
                }
            }
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::omega( const real T ) const
        {
            // collision integral, Table II. T* = T / ( epsilon / k )
            const real tLnTstar = std::log( T / mEpsilonKb ) ;

            real tSum = mB( 4 ) ;

            for( int k = 3; k >= 0; --k )
            {
                tSum = tSum * tLnTstar + mB( k ) ;
            }

            return std::exp( tSum ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::eta_0( const real T ) const
        {
            // Eq. ( 2 ), micro Pa s
            return mConstEta0 * std::sqrt( T ) / this->omega( T ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::eta_r() const
        {
            // Eq. ( 3 ), micro Pa s. c_i is zero for l_i = 0 and one otherwise
            const real tau   = this->get( BELFEM_LJTRANS_TAU );
            const real delta = this->get( BELFEM_LJTRANS_DELTA );

            real aValue = 0.0 ;

            for( uint k = 0; k < mEtaN.length(); ++k )
            {
                real tTerm =   mEtaN( k )
                             * std::pow( tau,   mEtaT( k ) )
                             * std::pow( delta, mEtaD( k ) ) ;

                if( mEtaL( k ) > 0.0 )
                {
                    tTerm *= std::exp( -std::pow( delta, mEtaL( k ) ) ) ;
                }

                aValue += tTerm ;
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::lambda_0( const real T ) const
        {
            /* Eq. ( 5 ), mW / ( m K ). The first term is the dilute gas
             * viscosity in micro Pa s divided by one micro Pa s, so the
             * unconverted value of eta_0 goes in here */
            const real tau = this->get( BELFEM_LJTRANS_TAU );

            return   mLambdaN1 * this->eta_0( T )
                   + mLambdaN2 * std::pow( tau, mLambdaT2 )
                   + mLambdaN3 * std::pow( tau, mLambdaT3 ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::lambda_r() const
        {
            // Eq. ( 6 ), mW / ( m K )
            const real tau   = this->get( BELFEM_LJTRANS_TAU );
            const real delta = this->get( BELFEM_LJTRANS_DELTA );

            real aValue = 0.0 ;

            for( uint k = 0; k < mLambdaN.length(); ++k )
            {
                real tTerm =   mLambdaN( k )
                             * std::pow( tau,   mLambdaT( k ) )
                             * std::pow( delta, mLambdaD( k ) ) ;

                if( mLambdaL( k ) > 0.0 )
                {
                    tTerm *= std::exp( -std::pow( delta, mLambdaL( k ) ) ) ;
                }

                aValue += tTerm ;
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::chi( const real T, const real v ) const
        {
            /* Eq. ( 11 ), symmetrised compressibility, dimensionless.
             *
             *   chi = pc * rho / rho_c^2 * ( d rho / d p )_T
             *
             * with rho = 1 / v the derivative is
             *
             *   ( d rho / d p )_T = -1 / ( v^2 * ( d p / d v )_T ) */
            const real tRho = 1.0 / v ;

            const real tDrhoDp = -1.0 / ( v * v * mEoS.dpdv( T, v ) ) ;

            return mPcrit * tRho / ( mRhocrit * mRhocrit ) * tDrhoDp ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::lambda_c(
                const real T, const real p ) const
        {
            // Olchowy and Sengers, Eqs. ( 7 ) to ( 11 ). Result in W / ( m K )
            constexpr real tR0    = 1.01 ;
            constexpr real tNu    = 0.63 ;
            constexpr real tGamma = 1.2415 ;

            const real v   = this->get( BELFEM_LJTRANS_V );
            const real rho = 1.0 / v ;

            /* Eq. ( 10 ). The reference state is taken at the same density but
             * at Tref, and the bracket is set to zero where it turns negative,
             * which happens at high temperature */
            const real tDeltaChi =
                    this->chi( T, v ) - this->chi( mTref, v ) * mTref / T ;

            if( tDeltaChi <= 0.0 )
            {
                return 0.0 ;
            }

            const real tXi = mXi0 * std::pow( tDeltaChi / mGammaC, tNu / tGamma ) ;

            const real tCp = mEoS.cp( T, p );
            const real tCv = mEoS.cv( T, p );

            const real tXiQd = tXi / mQd ;

            // Eq. ( 8 )
            const real tOmega = 2.0 / constant::pi
                    * (   ( tCp - tCv ) / tCp * std::atan( tXiQd )
                        + tCv / tCp * tXiQd ) ;

            // Eq. ( 9 )
            const real tOmega0 = 2.0 / constant::pi
                    * ( 1.0 - std::exp( -1.0 /
                        (   1.0 / tXiQd
                          + tXiQd * tXiQd / 3.0
                            * ( mRhocrit * mRhocrit ) / ( rho * rho ) ) ) ) ;

            /* Eq. ( 7 ). Everything here is SI -- rho in kg/m^3, cp in
             * J/(kg K), the viscosity in Pa s and xi in m -- so the result is
             * W / ( m K ) directly and must not be scaled again */
            const real tEta = this->mu( T, p ) ;

            /* constant::kB is the exact SI-2019 value; the paper quotes
             * 1.380658e-23. The 6.5e-6 relative difference is far below the
             * digits of the paper's own verification table */
            return rho * tCp * constant::kB * tR0 * T
                   / ( 6.0 * constant::pi * tEta * tXi )
                   * ( tOmega - tOmega0 ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::mu( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            if( ! this->test( BELFEM_LJTRANS_ETA ) )
            {
                // Eq. ( 1 ), micro Pa s converted to Pa s
                this->set( BELFEM_LJTRANS_ETA,
                        ( this->eta_0( T ) + this->eta_r() ) * 1.0e-6 );
            }

            return this->get( BELFEM_LJTRANS_ETA );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_LemmonJacobsen::lambda( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            if( ! this->test( BELFEM_LJTRANS_LAMBDA ) )
            {
                /* Eq. ( 4 ). The dilute gas and residual parts are in
                 * mW / ( m K ) and are scaled here; the critical enhancement
                 * already comes back in W / ( m K ) */
                this->set( BELFEM_LJTRANS_LAMBDA,
                        ( this->lambda_0( T ) + this->lambda_r() ) * 1.0e-3
                        + this->lambda_c( T, p ) );
            }

            return this->get( BELFEM_LJTRANS_LAMBDA );
        }

//----------------------------------------------------------------------------
    }
}
