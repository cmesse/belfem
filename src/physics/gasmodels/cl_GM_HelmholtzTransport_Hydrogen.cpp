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
#include "cl_GM_HelmholtzTransport_Hydrogen.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        HelmholtzTransport_Hydrogen::HelmholtzTransport_Hydrogen(
                Gas & aParent,
                const HelmholtzModel aModel ) :
            HelmholtzTransport( aParent ),
            mModel( aModel )
        {
            // Table 2, coefficients of Eq. ( 4 )
            mA = {  2.09630e-1, -4.55274e-1, 1.43602e-1,
                   -3.35325e-2,  2.76981e-3 } ;

            // Table 3, coefficients of Eq. ( 7 )
            mB = { -0.1870, 2.4871, 3.71513, -11.0972,
                    9.09655, -3.8292, 0.5166 } ;

            // Table 4, coefficients of Eq. ( 9 )
            mC = { 6.43449673, 4.56334068e-02, 2.32797868e-01,
                   9.58326120e-01, 1.27941189e-01, 3.63576595e-01 } ;

            // Eq. ( 3 ). M in g/mol and sigma in nm, giving micro Pa s
            mConstEta0 = 0.021357 * std::sqrt( mM ) / ( mSigma * mSigma ) ;

            /* Eq. ( 6 ) with the erratum's Avogadro factor. sigma is converted
             * from nm to m and the molar mass from g/mol to kg/mol, so the
             * second viscosity virial comes out in m^3/kg and multiplies a
             * density in kg/m^3 directly.
             *
             * The erratum quotes N_A = 6.022137e23 /mol as the value current
             * when the correlation was fitted; constant::NA is the exact
             * SI-2019 value and differs by 1.2e-7 relative, far below the
             * digits of the erratum's own test values. */
            mConstBeta = std::pow( mSigma * 1.0e-9, 3 ) * constant::NA
                         / ( mM * 1.0e-3 ) ;

            /* Thermal conductivity, Assael et al. The reduced variables of
             * Eqs. ( 2 ) and ( 3 ) are formed with the critical point of the
             * equation of state, which is Leachman et al., the same equation
             * of state the correlation was fitted against */
            mTcritEoS   = aParent.component( 0 )->data()->T_crit() ;
            mPcritEoS   = aParent.component( 0 )->data()->p_crit() ;
            mRhocritEoS = aParent.component( 0 )->data()->rho_crit() ;

            mTrefC = 1.5 * mTcritEoS ;

            if( mModel == HelmholtzModel::ParaHydrogen )
            {
                // Assael Table 5
                mLambdaA1 = { -1.24500e0,  3.10212e2, -3.31004e2,  2.46016e2,
                              -6.57810e1,  1.08260e1, -5.19659e-1, 1.43979e-2 } ;
                mLambdaA2 = {  1.42304e4, -1.93922e4,  1.58379e4, -4.81812e3,
                               7.28639e2, -3.57365e1,  1.00000e0 } ;
                mLambdaB1 = {  2.65975e-2, -1.33826e-3, 1.30219e-2,
                              -5.67678e-3, -9.23380e-5 } ;
                mLambdaB2 = { -1.21727e-3,  3.66663e-3, 3.88715e-3,
                              -9.21055e-3,  4.00723e-3 } ;
            }
            else
            {
                /* Assael Table 2. Orthohydrogen has no correlation of its own
                 * and is served by the normal hydrogen tables, as REFPROP does */
                mLambdaA1 = { -3.40976e-1, 4.58820e0, -1.45080e0, 3.26394e-1,
                               3.16939e-3, 1.90592e-4, -1.13900e-6 } ;
                mLambdaA2 = {  1.38497e2, -2.21878e1, 4.57151e0, 1.00000e0 } ;
                mLambdaB1 = {  3.63081e-2, -2.07629e-2, 3.14810e-2,
                              -1.43097e-2,  1.74980e-3 } ;
                mLambdaB2 = {  1.83370e-3, -8.86716e-3, 1.58260e-2,
                              -1.06283e-2,  2.80673e-3 } ;
            }
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::lambda_0( const real T ) const
        {
            // Assael Eq. ( 2 ), a rational function of T / T_crit, W/(m K)
            const real x = T / mTcritEoS ;

            real tNum = 0.0 ;

            for( int k = mLambdaA1.length() - 1; k >= 0; --k )
            {
                tNum = tNum * x + mLambdaA1( k ) ;
            }

            real tDen = 0.0 ;

            for( int k = mLambdaA2.length() - 1; k >= 0; --k )
            {
                tDen = tDen * x + mLambdaA2( k ) ;
            }

            return tNum / tDen ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::lambda_e(
                const real T, const real rho ) const
        {
            // Assael Eq. ( 3 ), W/(m K)
            const real x = T / mTcritEoS ;
            const real d = rho / mRhocritEoS ;

            real aValue = 0.0 ;
            real tPow   = 1.0 ;

            for( uint k = 0; k < 5; ++k )
            {
                tPow   *= d ;
                aValue += ( mLambdaB1( k ) + mLambdaB2( k ) * x ) * tPow ;
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::chi( const real T, const real v ) const
        {
            // the bracket of Assael Eq. ( 7 ), dimensionless
            const real tRho    = 1.0 / v ;
            const real tDrhoDp = -1.0 / ( v * v * mEoS.dpdv( T, v ) ) ;

            return mPcritEoS * tRho / ( mRhocritEoS * mRhocritEoS ) * tDrhoDp ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::lambda_c(
                const real T, const real p, const real rho ) const
        {
            /* Simplified Olchowy-Sengers crossover, Assael Eqs. ( 4 ) to ( 7 ).
             * Same shape as the Lemmon-Jacobsen enhancement, but this paper
             * writes the cut-off as a wavenumber, so the products are qD * xi
             * where that one has xi / qD */
            constexpr real tRD    = 1.01 ;
            constexpr real tNu    = 0.63 ;
            constexpr real tGamma = 1.2415 ;

            const real v = 1.0 / rho ;

            const real tDeltaChi =
                    this->chi( T, v ) - this->chi( mTrefC, v ) * mTrefC / T ;

            if( tDeltaChi <= 0.0 )
            {
                return 0.0 ;
            }

            const real tXi = mXi0 * std::pow( tDeltaChi / mGammaC, tNu / tGamma ) ;

            const real tCp = mEoS.cp( T, p );
            const real tCv = mEoS.cv( T, p );

            const real tQdXi = mQd * tXi ;

            // Eq. ( 5 )
            const real tOmega = 2.0 / constant::pi
                    * (   ( tCp - tCv ) / tCp * std::atan( tQdXi )
                        + tCv / tCp * tQdXi ) ;

            // Eq. ( 6 )
            const real tOmega0 = 2.0 / constant::pi
                    * ( 1.0 - std::exp( -1.0 /
                        (   1.0 / tQdXi
                          + tQdXi * tQdXi / 3.0
                            * ( mRhocritEoS * mRhocritEoS ) / ( rho * rho ) ) ) ) ;

            // Eq. ( 4 ), SI throughout, so the result is W/(m K)
            return rho * tCp * tRD * constant::kB * T
                   / ( 6.0 * constant::pi * this->mu( T, p ) * tXi )
                   * ( tOmega - tOmega0 ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::lambda( const real T, const real p ) const
        {
            const real tRho = 1.0 / mEoS.v( T, p ) ;

            // Assael Eq. ( 1 ), every term already in W/(m K)
            return   this->lambda_0( T )
                   + this->lambda_e( T, tRho )
                   + this->lambda_c( T, p, tRho ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::eta_0( const real T ) const
        {
            // Eq. ( 4 ), reduced effective cross section
            const real tLnTstar = std::log( T / mEpsilonKb ) ;

            real tSum = mA( 4 ) ;

            for( int k = 3; k >= 0; --k )
            {
                tSum = tSum * tLnTstar + mA( k ) ;
            }

            // Eq. ( 3 ), micro Pa s
            return mConstEta0 * std::sqrt( T ) / std::exp( tSum ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::eta_1( const real T ) const
        {
            const real tTstar = T / mEpsilonKb ;

            /* Eq. ( 7 ). The exponent on T* is NEGATIVE; the paper as printed
             * has it positive and the erratum corrects it */
            real tBstar = mB( 0 ) ;
            real tPow   = 1.0 ;

            for( uint k = 1; k < 7; ++k )
            {
                tPow   /= tTstar ;
                tBstar += mB( k ) * tPow ;
            }

            // Eqs. ( 5 ) and ( 6 )
            return tBstar * mConstBeta * this->eta_0( T ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Hydrogen::mu( const real T, const real p ) const
        {
            // the correlation is written on temperature and density
            const real tRho = 1.0 / mEoS.v( T, p ) ;

            const real tTr   = T / mTcrit ;
            const real tRhoR = tRho / mRhoScale ;

            const real tRhoR2 = tRhoR * tRhoR ;

            // Eq. ( 9 ), micro Pa s
            const real aValue =
                      this->eta_0( T )
                    + this->eta_1( T ) * tRho
                    + mC( 0 ) * tRhoR2 * std::exp(
                            mC( 1 ) * tTr
                          + mC( 2 ) / tTr
                          + mC( 3 ) * tRhoR2 / ( mC( 4 ) + tTr )
                          + mC( 5 ) * tRhoR2 * tRhoR2 * tRhoR2 ) ;

            // micro Pa s -> Pa s
            return aValue * 1.0e-6 ;
        }

//----------------------------------------------------------------------------
    }
}
