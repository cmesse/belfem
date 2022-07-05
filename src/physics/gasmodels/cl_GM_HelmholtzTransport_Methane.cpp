//
// Created by Christian Messe on 23.11.20.
//
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
        HelmholtzTransport_Methane::mu( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            if ( ! this->test( BELFEM_METHANE_MU ) )
            {
                this->set( BELFEM_METHANE_MU, this->eta_0() + this->eta_ex() );
            }
            return this->get( BELFEM_METHANE_MU );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            if ( ! this->test( BELFEM_METHANE_LAMBDA ) )
            {
                this->set( BELFEM_METHANE_LAMBDA,
                           this->lambda_0() + this->lambda_cr( aT, aP ) + this->lambda_ex() );
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


            // for Eq. ( 18 )
            // the original paper fits the data to Lambda* = 2.235e9
            // however, we have a slightly different critical point.
            // Therefore, the value is corrected to  2.231293e9
            mConstLambdaCr = 2.231293e9 * constant::kB * mR * mR * mTcrit * mTcrit /
                    ( 6.0 * constant::pi * mPcrit * mVcrit * mVcrit );

            mConstExpLambdaCr = ( 1.190 - 0.633 ) / 1.190 ;

            // for Eq. ( 17 )
            mConstLambdaEx =   std::pow( mPcrit, 2.0 / 3.0 )
                             * std::pow( constant::kB, 5.0 / 6.0 ) /
                    ( std::pow( mTcrit, 1.0 / 6.0 ) * std::sqrt( mM / constant::NA ) );

        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::eta_0()
        {
            if ( ! this->test( BELFEM_METHANE_MU0 ) )
            {
                const real & aT = this->get( BELFEM_METHANE_T );

                this->set( BELFEM_METHANE_MU0, mConstEta0 * std::sqrt( aT ) / this->omega() );
            }

            // Eq. ( 10 )
            return this->get( BELFEM_METHANE_MU0 );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::eta_ex()
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
        HelmholtzTransport_Methane::omega()
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
        HelmholtzTransport_Methane::lambda_0()
        {

            const real & tTau = this->get( BELFEM_METHANE_TAU );
            real tFint = mF( 0 ) + mF( 1 ) * mEpsilonKb / this->get( BELFEM_METHANE_T ) ;

            return mR * this->eta_0() * ( 3.75 - tFint * ( tTau * tTau * mEoS.phi0_tt() + 1.5 ) );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda_cr( const real & aT, const real & aP  )
        {
            const real & tTau   = this->get( BELFEM_METHANE_TAU );
            const real & tDelta = this->get( BELFEM_METHANE_DELTA );
            return mConstLambdaCr * std::pow(
                    1.0 + tDelta * ( mEoS.phir_d() - tTau * mEoS.phir_dt() ), 2 ) /
                    ( tTau * tTau * this->mu( aT, aP ) ) ;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::lambda_ex()
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
        HelmholtzTransport_Methane::chi()
        {
            // Eq. ( 19 )
            const real & tTau   = this->get( BELFEM_METHANE_TAU);
            const real & tDelta = this->get( BELFEM_METHANE_DELTA );

            return mZcrit * tDelta * tTau /
                ( 1 + tDelta * ( 2 * mEoS.phir_d() + tDelta * mEoS.phir_dd() ) );
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport_Methane::f()
        {
            real tTstar   = std::abs( 1.0 - 1.0 / this->get( BELFEM_METHANE_TAU ) );
            real tRhostar = 1 - this->get( BELFEM_METHANE_DELTA );

            // Eq. ( 22 )
            return std::exp( - ( 2.646 * std::sqrt( tTstar )
             + tRhostar * ( tRhostar * 2.678 - 0.637 ) ) );
        }

//----------------------------------------------------------------------------

    }
}