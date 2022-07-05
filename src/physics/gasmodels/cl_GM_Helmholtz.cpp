//
// Created by Christian Messe on 17.08.20.
//

#include "cl_GM_Helmholtz.hpp"
#include "assert.hpp"
#include "cl_Gas.hpp"
#include "fn_linspace.hpp"
#include "fn_polyval.hpp"
#include "fn_polyfit.hpp"
#include "cl_GT_RefGas.hpp"
namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        Helmholtz::Helmholtz( Gas & aParent, const string & aLabel ) :
                EoS( aParent ),
                mLabel( aLabel )
        {
            mCubicEoS = new EoS_Cubic( aParent, GasModel::SRK );
            mCubicEoS->remix() ;
        }

//----------------------------------------------------------------------------

        void
        Helmholtz::delete_cubic_eos()
        {
            delete mCubicEoS ;
        }

//----------------------------------------------------------------------------

        void
        Helmholtz::set_reference_point()
        {
            // synchronize offset
            this->set_reference_point( gastables::gTref, gastables::gPref );

            // grab reference gas
            gastables::RefGas * tRefGas = mParent.component( 0 );
            BELFEM_ERROR( tRefGas->label() == mLabel,
                         "Parent has wrong label ( is %s, expect %s )",
                         tRefGas->label().c_str(), mLabel.c_str() );

            // grab data object
            const gastables::GasData * tData = tRefGas->data() ;

            // offset
            mH0 += tData->Href() / tData->M() ;
            mU0 += tData->Href() / tData->M() ;
            mS0 += tData->Sref() / tData->M() ;
        }

//----------------------------------------------------------------------------

        void
        Helmholtz::set_critical_point_in_data_object()
        {
            // get data object
            gastables::GasData * tData = mParent.component( 0 )->data() ;

            tData->set_t_crit( mTcrit );
            tData->set_p_crit( mPcrit );
            tData->set_rho_crit( mRhocrit );
            tData->set_molar_mass( mM );
        }

//----------------------------------------------------------------------------

        // initialize the offsets for enthalpy and entropy
        void
        Helmholtz::set_reference_point( const real aTref, const real aPref )
        {

            mS0 = 0.0 ;
            mH0 = 0.0 ;
            mU0 = 0.0 ;

            real tH0 = this->h( aTref, aPref ) ;
            real tS0 = this->s( aTref, aPref );

            mH0 = -tH0 ;
            mU0 = mH0 - aPref * this->v( aTref, aPref );
            mS0 = -tS0 ;
        }

//----------------------------------------------------------------------------

        void
        Helmholtz::remix()
        {
            mCubicEoS->remix();
        }

//----------------------------------------------------------------------------

        void
        Helmholtz::eval_critical_point( real & aT, real & aP, real & aV )
        {
            aT = mTcrit ;
            aP = mPcrit ;
            aV = mVcrit ;
        }

//----------------------------------------------------------------------------
        void
        Helmholtz::init_Tvap_poly()
        {
            uint tN = 100 ;
            Vector< real > tT( tN );
            Vector< real > tX( tN );

            linspace( mTtriple, mTcrit, tN, tT );

            for( uint k=0; k<tN; ++k )
            {
                tX( k ) = std::log( this->p_vap( tT( k ) ) / mPcrit );
            }

            polyfit( tX, tT, 3, mTvap );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::phi( const real & aT, const real & aV )
        {
            this->update_Tv( aT, aV );

            return this->phi0() + this->phir();
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::p( const real & aT, const real & aV )
        {
            this->update_Tv( aT, aV );

            return mR * aT * ( 1.0 + mDelta * this->phir_d() ) / aV ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dpdv( const real & aT, const real & aV )
        {
            this->update_Tv( aT, aV );

            if( ! this->test( BELFEM_HELMHOLTZ_DPDV ) )
            {
                this->set( BELFEM_HELMHOLTZ_DPDV ,
                     - ( mR * aT / mVcrit ) *
                     ( mDelta * ( mDelta * this->phir_dd() + 2.0 * this->phir_d() ) + 1.0 )
                     * mDelta / aV  );
            }

            return this->get( BELFEM_HELMHOLTZ_DPDV );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dpdT( const real & aT, const real & aV )
        {
            this->update_Tv( aT, aV );

            if( ! this->test( BELFEM_HELMHOLTZ_DPDT ) )
            {
                this->set( BELFEM_HELMHOLTZ_DPDT,
                ( mR / aV ) * ( ( 1.0 + mDelta * this->phir_d() )
                    - mTau * mDelta * this->phir_dt() ) );

            }
            return  this->get( BELFEM_HELMHOLTZ_DPDT );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dvdT( const real & aT, const real & aV )
        {
            this->update_Tv( aT, aV );

            if( ! this->test( BELFEM_HELMHOLTZ_DVDT ) )
            {
                this->set( BELFEM_HELMHOLTZ_DVDT,
                           - this->dpdT( aT, aV ) / this->dpdv( aT, aV ) );
            }
            return  this->get( BELFEM_HELMHOLTZ_DVDT );
        }
//----------------------------------------------------------------------------

        real
        Helmholtz::alpha( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );
            return aP * this->beta( aT, aP ) * this->kappa( aT, aP );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::beta( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );
            const real & tV = mHelmholtzVals[ BELFEM_HELMHOLTZ_V ];

            return this->dpdT( aT, tV ) / aP ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::kappa( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );
            const real & tV = mHelmholtzVals[ BELFEM_HELMHOLTZ_V ];
            return -1.0 / ( tV * this->dpdv( aT, tV ) ) ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::v( const real & aT, const real & aP )
        {
            real aV ;

            // relaxation factor
            real tOmega ;

            // set the liquid flag of the parent
            mParent.set_liquid_flag( this->is_liquid( aT, aP ) );

            if( this->is_liquid( aT, aP ) )
            {
                aV = polyval( mVliq, aT ) ;
                tOmega = 0.5 ;
            }
            else
            {
                aV = mCubicEoS->v( aT, aP );

                // empirical factor
                tOmega = std::max( 1.0 - 2.0 * std::exp( - aT / mTcrit ), 0.1 );
            }

            real tF = BELFEM_REAL_MAX ;
            real tDeltaV ;

            uint tCount = 0;

            while ( std::abs( tF / aP ) > 1e-8 )
            {
                tF = this->p( aT, aV ) - aP ;

                tDeltaV = tF / this->dpdv( aT, aV );

                // limit step
                if ( tDeltaV > aV )
                {
                    tDeltaV = aV ;
                }

                aV -= tOmega * tDeltaV ;

                BELFEM_ERROR( tCount++ < 200,
                             "Too many iterations for T=%8.3f K, p=%8.3f bar, rho=%8.3f kg/m^3, relax=%8.3f",
                             ( double ) aT, ( double ) aP * 1e-5, ( double ) 1.0 / aV, ( double ) tOmega
                             );
            }

            return aV ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::T( const real & aP, const real & aV )
        {
            real aT ;

            // check if we might know the solution already
            real tA = std::abs( aP - this->get( BELFEM_HELMHOLTZ_P ) ) / aP ;
            real tB = std::abs( aV - this->get( BELFEM_HELMHOLTZ_V ) ) / aV ;

            // use data from memory as initial guess
            if( tA < 0.25 && tB < 0.25 )
            {
                aT = this->get( BELFEM_HELMHOLTZ_T );
            }
            else
            {
                aT = std::max( aP * aV / mR, mTtriple );
            }

            real tF = BELFEM_REAL_MAX ;
            real tdF ;
            real tdT ;
            uint tCount = 0 ;
            real tOmega = 0.99 ;

            while( std::abs( tF / aP ) > 1e-8 )
            {

                tF  = this->p( aT, aV ) - aP ;
                tdF = this->dpdT( aT, aV ) ;

                tdT = tF / tdF ;
                if( tdT > ( aT - mTtriple ) )
                {
                    tdT = aT - mTtriple ;
                    tOmega = 0.3 ;
                }

                aT -= tOmega * tdT ;

                BELFEM_ERROR( tCount++ < 100,
                             "Iteration failure in T(p,v) for p=%8.3f Pa v=%10.4f m^3/kg",
                             aP, aV );
            }

            return aT ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::u( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            return mTau * ( this->phi0_t() + this->phir_t() )  * mR * aT + mU0 ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::h( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            return ( 1.0 + mTau * ( this->phi0_t() + this->phir_t() )
                     + mDelta * this->phir_d() ) * mR * aT  + mH0 ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::s( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            return ( mTau * ( this->phi0_t() + this->phir_t() )
                  - this->phi0()  - this->phir() ) * mR + mS0 ;

        }

//----------------------------------------------------------------------------

        real
        Helmholtz::cv( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );
            return - mTau * mTau * ( this->phi0_tt() + this->phir_tt() ) * mR ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::cp( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            return ( std::pow( ( 1.0 + mDelta * ( this->phir_d() - mTau * this->phir_dt() ) ), 2 ) /
                 ( 1.0 + mDelta * ( 2.0 * this->phir_d() + mDelta * this->phir_dd()))
                   - mTau * mTau * ( this->phi0_tt() + this->phir_tt()  ) ) * mR ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::w( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            real tA = 1.0 + mDelta * ( 2.0 * this->phir_d() + mDelta * this->phir_dd());
            real tB = 1.0 + mDelta * ( this->phir_d() - mTau * this->phir_dt());
            real tC = mTau * mTau * ( this->phi0_tt() + this->phir_tt());

            return std::sqrt(( tA - tB * tB / tC ) * ( mR * aT ));
        }
//----------------------------------------------------------------------------

        real
        Helmholtz::dsdT( const real & aT, const real & aP )
        {
            return this->cp( aT, aP ) / aT ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dsdp( const real & aT, const real & aP )
        {
            this->update_Tp( aT, aP );

            return - this->dvdT( aT, mHelmholtzVals[ BELFEM_HELMHOLTZ_V ] );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phi0()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHI0 ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHI0,
                           this->compute_phi0() ) ;
            }

            return this->get( BELFEM_HELMHOLTZ_PHI0 );
        }


//----------------------------------------------------------------------------

        const real &
        Helmholtz::phi0_t()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHI0_T ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHI0_T,
                           this->compute_phi0_t() ) ;
            }

            return this->get( BELFEM_HELMHOLTZ_PHI0_T );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phi0_tt()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHI0_TT ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHI0_TT,
                           this->compute_phi0_tt() ) ;
            }

            return this->get( BELFEM_HELMHOLTZ_PHI0_TT );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR, this->compute_phir() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_d()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_D ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_D, this->compute_phir_d() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_D );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_dd()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_DD ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_DD, this->compute_phir_dd() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_DD );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_t()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_T ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_T, this->compute_phir_t() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_T );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_tt()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_TT ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_TT, this->compute_phir_tt() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_TT );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_dt()
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_DT ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_DT, this->compute_phir_dt() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_DT );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phi0()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phi0() not implemented for this Helmholtz model" );
            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir() not implemented for this Helmholtz model" );
            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phi0_t()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phi0_t() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phi0_tt()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phi0_tt() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_d()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_d() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_dd()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_dd() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_t()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_t() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_tt()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_tt() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_dt()
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_dt() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::p_vap( const real & aT )
        {
            return std::exp( this->pi_vap( aT ) ) * mPcrit ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::T_vap( const real & aP )
        {
            // compute initial guess
            real aT = polyval( mTvap, std::log( aP / mPcrit ) );

            real tDeltaT = BELFEM_REAL_MAX ;
            real tPvap ;
            real tPivap ;
            while( std::abs( tDeltaT ) > BELFEM_EPSILON_T )
            {
                tPivap = this->pi_vap( aT );
                tPvap = std::exp( tPivap ) * mPcrit ;
                tDeltaT = ( tPvap - aP ) / this->dpvap_dT( aT, tPvap, tPivap ) ;

                aT -= 0.99 * tDeltaT ;
            }

            return aT ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::hvap( const real & aT, const real & aP )
        {
            if( aT < mTcrit )
            {
                real tTvap = this->T_vap( aP );
                real aHliq = this->h( tTvap - BELFEM_EPSILON_T, aP );
                real aHgas = this->h( tTvap + BELFEM_EPSILON_T, aP );
                return aHgas - aHliq ;
            }
            else
            {
                return 0.0 ;
            }
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dpvap_dT( const real & aT, const real aPvap, const real aPiVap )
        {
            // real tG = aPvap / mPcrit ;
            real tdG = - ( aPiVap + this->psi_vap( aT ) )/ aT;

            return ( aT > mTcrit ) ? -aPvap * tdG : aPvap * tdG ;
        }

//----------------------------------------------------------------------------

        bool
        Helmholtz::is_liquid( const real & aT, const real & aP )
        {
            if( aT < mTtriple )
            {
                return true ;
            }
            else if ( aT > mTcrit )
            {
                return false ;
            }
            else
            {
                return aP > this->p_vap( aT );
            }
        }

//----------------------------------------------------------------------------
    }
}