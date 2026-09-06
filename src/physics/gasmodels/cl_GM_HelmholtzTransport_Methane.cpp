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
#include "cl_GM_HelmholtzTransport_Methane.hpp"


namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        HelmholtzTransport_Methane::HelmholtzTransport_Methane( Gas & aParent ) :
            HelmholtzTransport( aParent ),
            mTcrit( aParent.component( 0 )->data()->T_crit() ),
            mPcrit( aParent.component( 0 )->data()->p_crit() ),
            mM( aParent.component( 0 )->data()->M() ),
            mVcrit( 1.0 / aParent.component( 0 )->data()->rho_crit() ),
            mR( constant::Rm / aParent.component( 0 )->data()->M() ),
            mZcrit( aParent.component( 0 )->data()->Z_crit() )
        {
            this->init_tables();
            this->init_constants() ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::mu( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            if ( ! this->test( BELFEM_METHANE_MU ) )
            {
                this->set( BELFEM_METHANE_MU, this->eta_0() + this->eta_ex() );
            }
            return this->get( BELFEM_METHANE_MU );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            if ( ! this->test( BELFEM_METHANE_LAMBDA ) )
            {
                this->set( BELFEM_METHANE_LAMBDA,
                           this->lambda_0() + this->lambda_cr( T, p ) + this->lambda_ex() );
            }
            return this->get( BELFEM_METHANE_LAMBDA );;
        }

//----------------------------------------------------------------------------

        void
        HelmholtzTransport_Methane::init_tables()
        {
            mOmega.set_size( 9 );

            mOmega( 0 ) = -3.0328138281 ;
            mOmega( 1 ) =  16.918880086 ;
            mOmega( 2 ) = -37.189364917 ;
            mOmega( 3 ) =  41.288861858 ;
            mOmega( 4 ) = -24.615921140 ;
            mOmega( 5 ) =   8.9488430959 ;
            mOmega( 6 ) = -1.8739245042 ;
            mOmega( 7 ) =  0.2096610139 ;
            mOmega( 8 ) = -0.0096570437074 ;

            mF.set_size( 2 );
            mF( 0 ) = 1.45885 ;
            mF( 1 ) = -0.4377162 ;

            // Table 9, values for eta_Ex
            mA.set_size( 11 );
            mA(  0 ) = 1. ;
            mA(  1 ) = 1. ;
            mA(  2 ) = 2. ;
            mA(  3 ) = 2. ;
            mA(  4 ) = 2. ;
            mA(  5 ) = 3. ;
            mA(  6 ) = 3. ;
            mA(  7 ) = 4. ;
            mA(  8 ) = 4. ;
            mA(  9 ) = 1. ;
            mA( 10 ) = 1. ;

            mB.set_size( 11 );
            mB(  0 ) = 0.0 ;
            mB(  1 ) = 1.0 ;
            mB(  2 ) = 0.0 ;
            mB(  3 ) = 1.0 ;
            mB(  4 ) = 1.5 ;
            mB(  5 ) = 0.0 ;
            mB(  6 ) = 2.0 ;
            mB(  7 ) = 0.0 ;
            mB(  8 ) = 1.0 ;
            mB(  9 ) = 0.0 ;
            mB( 10 ) = 1.0 ;

            mG.set_size( 11 );
            mG(	 0	) =	 0.41250137	;
            mG(	 1	) =	-0.14390912	;
            mG(	 2	) =	 0.10366993	;
            mG(	 3	) =	 0.40287464	;
            mG(	 4	) =	-0.24903524	;
            mG(	 5	) =	-0.12953131	;
            mG(	 6	) =	 0.06575776	;
            mG(	 7	) =	 0.02566628	;
            mG(	 8	) =	-0.03716526	;
            mG(	 9	) =	-0.38798341	;
            mG(	10	) =	 0.03533815	;

            // Table 9, Values for lambda_ex
            mC.set_size( 7 );
            mC( 0 ) =  1. ;
            mC( 1 ) =  3. ;
            mC( 2 ) =  4. ;
            mC( 3 ) =  4. ;
            mC( 4 ) =  5. ;
            mC( 5 ) =  5. ;
            mC( 6 ) =  2. ;

            mD.set_size( 7 );
            mD( 0 ) =  0. ;
            mD( 1 ) =  0. ;
            mD( 2 ) =  0. ;
            mD( 3 ) =  1. ;
            mD( 4 ) =  0. ;
            mD( 5 ) =  1. ;
            mD( 6 ) =  0. ;

            mJ.set_size( 7 );
            mJ( 0 ) =   2.4149207 ;
            mJ( 1 ) =   0.55166331 ;
            mJ( 2 ) =  -0.52837734 ;
            mJ( 3 ) =   0.073809553 ;
            mJ( 4 ) =   0.24465507 ;
            mJ( 5 ) =  -0.047613626 ;
            mJ( 6 ) =   1.5554612 ;
        }

//----------------------------------------------------------------------------

        void
        HelmholtzTransport_Methane::init_constants()
        {
            // for Eq. ( 10 )
            mConstEta0 = std::sqrt( mM * constant::kB /
                                    ( constant::NA * constant::pi  ) )
                    * 5.0 / ( 16.0 * mSigma * mSigma ) ;


            // for Eq. ( 15 )
            mConstEtaEx =   std::pow( mPcrit, 2.0/ 3.0 )
                          * std::sqrt( mM / constant::NA ) /
                            std::pow( mTcrit * constant::kB, 1.0 / 6.0 );


            /* for Eq. ( 18 ). Friend fits Lambda* = 2.235e9 1/m against his own
             * critical point ( Tc = 190.551 K, pc = 4.5992 MPa, rhoc = 162.660 kg/m3 ).
             * We run on the Setzmann & Wagner critical point instead; Tc cancels in
             * lambda_cr ( Tc^2 / tau^2 = T^2 ), and pc and rhoc are identical, so the
             * fitted value applies unchanged. */
            mConstLambdaCr = 2.235e9 * constant::kB * mR * mR * mTcrit * mTcrit /
                    ( 6.0 * constant::pi * mPcrit * mVcrit * mVcrit );

            mConstExpLambdaCr = ( 1.190 - 0.633 ) / 1.190 ;

            // for Eq. ( 17 )
            mConstLambdaEx =   std::pow( mPcrit, 2.0 / 3.0 )
                             * std::pow( constant::kB, 5.0 / 6.0 ) /
                    ( std::pow( mTcrit, 1.0 / 6.0 ) * std::sqrt( mM / constant::NA ) );

        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::eta_0() const
        {
            if ( ! this->test( BELFEM_METHANE_MU0 ) )
            {
                const real T = this->get( BELFEM_METHANE_T );

                this->set( BELFEM_METHANE_MU0, mConstEta0 * std::sqrt( T ) / this->omega() );
            }

            // Eq. ( 10 )
            return this->get( BELFEM_METHANE_MU0 );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::eta_ex() const
        {
            const real & tTau = this->get( BELFEM_METHANE_TAU );

            const real & tDelta = this->get( BELFEM_METHANE_DELTA );

            real tA = 0.0 ;
            for( uint k=0; k<9; ++k )
            {
                tA += mG( k ) * std::pow( tDelta, mA( k ) ) * std::pow( tTau, mB( k ) );
            }

            real tB = 1.0 ;
            for( uint k=9; k<11; ++k )
            {
                tB += mG( k ) * std::pow( tDelta, mA( k ) ) * std::pow( tTau, mB( k ) );
            }

            return mConstEtaEx * tA / tB ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::omega() const
        {
            // Eq. ( 12 )
            real aInvOmega = 0.0 ;

            real tTred = this->get( BELFEM_METHANE_T ) / mEpsilonKb ;

            for ( uint k=0; k<9; ++k )
            {
                aInvOmega += mOmega( k ) * std::pow( tTred, k / 3.0 - 1.0 ) ;
            }

            return 1.0 / aInvOmega ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda_0() const
        {

            const real & tTau = this->get( BELFEM_METHANE_TAU );
            real tFint = mF( 0 ) + mF( 1 ) * mEpsilonKb / this->get( BELFEM_METHANE_T ) ;

            return mR * this->eta_0() * ( 3.75 - tFint * ( tTau * tTau * mEoS.phi0_tt() + 1.5 ) );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda_cr( const real T, const real p  ) const
        {
            const real & tTau   = this->get( BELFEM_METHANE_TAU );
            const real & tDelta = this->get( BELFEM_METHANE_DELTA );

            /* Eq. ( 18 ): the susceptibility factor chi^((gamma-nu)/gamma) and the
             * crossover function F confine the enhancement to the critical region;
             * without them the term is a spurious background at all conditions.
             * chi ~ (dp/drho)^-1 can turn negative in metastable states, where
             * pow() would return NaN; clamp to zero, which kills the term. */
            const real tChi = this->chi() ;

            if( tChi <= 0.0 )
            {
                return 0.0 ;
            }

            return mConstLambdaCr * std::pow(
                    1.0 + tDelta * ( mEoS.phir_d() - tTau * mEoS.phir_dt() ), 2 ) /
                    ( tTau * tTau * this->mu( T, p ) )
                    * std::pow( tChi, mConstExpLambdaCr ) * this->f() ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda_ex() const
        {

            const real & tT = this->get( BELFEM_METHANE_T );
            const real & tV = this->get( BELFEM_METHANE_V );
            const real & tTau = this->get( BELFEM_METHANE_TAU );
            const real & tDelta = this->get( BELFEM_METHANE_DELTA );


            real aValue;
            if ( ( tT < mTcrit ) && ( tV > mVcrit ) )
            {
                // Eq. ( 16 )
                real tDeltaSigma = mVcrit / mEoS.v( tT, mEoS.p_vap( tT ) );

                aValue = mJ( 6 ) * tDelta * tDelta / tDeltaSigma ;
            }
            else
            {
                aValue = mJ( 6 ) * tDelta * tDelta ;
            }

            for ( int k=0; k<6; ++k )
            {
                aValue += mJ( k ) * std::pow( tDelta, mC( k ) ) * std::pow( tTau, mD( k ) );
            }

            return aValue * mConstLambdaEx ;
         }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::chi() const
        {
            // Eq. ( 19 )
            const real & tTau   = this->get( BELFEM_METHANE_TAU);
            const real & tDelta = this->get( BELFEM_METHANE_DELTA );

            return mZcrit * tDelta * tTau /
                ( 1 + tDelta * ( 2 * mEoS.phir_d() + tDelta * mEoS.phir_dd() ) );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::f() const
        {
            real tTstar   = std::abs( 1.0 - 1.0 / this->get( BELFEM_METHANE_TAU ) );
            real tRhostar = 1 - this->get( BELFEM_METHANE_DELTA );

            // Eq. ( 20 ); T* and rho* are Eqs. ( 21 ) and ( 22 )
            return std::exp( - ( 2.646 * std::sqrt( tTstar )
             + tRhostar * ( tRhostar * 2.678 - 0.637 ) ) );
        }

//----------------------------------------------------------------------------

    }
}