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

            // h_raw - u_raw = p*v identically ( both derive from the same phi ),
            // so u = h - p*v requires the same offset for both
            mH0 = -tH0 ;
            mU0 = mH0 ;
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
        Helmholtz::eval_critical_point( real & T, real & p, real & v ) const
        {
            T = mTcrit ;
            p = mPcrit ;
            v = mVcrit ;
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
        Helmholtz::phi( const real T, const real v ) const
        {
            this->update_Tv( T, v );

            return this->phi0() + this->phir();
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::p( const real T, const real v ) const
        {
            this->update_Tv( T, v );

            return mR * T * ( 1.0 + mDelta * this->phir_d() ) / v ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dpdv( const real T, const real v ) const
        {
            this->update_Tv( T, v );

            if( ! this->test( BELFEM_HELMHOLTZ_DPDV ) )
            {
                this->set( BELFEM_HELMHOLTZ_DPDV ,
                     - ( mR * T / mVcrit ) *
                     ( mDelta * ( mDelta * this->phir_dd() + 2.0 * this->phir_d() ) + 1.0 )
                     * mDelta / v  );
            }

            return this->get( BELFEM_HELMHOLTZ_DPDV );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dpdT( const real T, const real v ) const
        {
            this->update_Tv( T, v );

            if( ! this->test( BELFEM_HELMHOLTZ_DPDT ) )
            {
                this->set( BELFEM_HELMHOLTZ_DPDT,
                ( mR / v ) * ( ( 1.0 + mDelta * this->phir_d() )
                    - mTau * mDelta * this->phir_dt() ) );

            }
            return  this->get( BELFEM_HELMHOLTZ_DPDT );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dvdT( const real T, const real v ) const
        {
            this->update_Tv( T, v );

            if( ! this->test( BELFEM_HELMHOLTZ_DVDT ) )
            {
                this->set( BELFEM_HELMHOLTZ_DVDT,
                           - this->dpdT( T, v ) / this->dpdv( T, v ) );
            }
            return  this->get( BELFEM_HELMHOLTZ_DVDT );
        }
//----------------------------------------------------------------------------

        real
        Helmholtz::alpha( const real T, const real p ) const
        {
            this->update_Tp( T, p );
            return p * this->beta( T, p ) * this->kappa( T, p );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::beta( const real T, const real p ) const
        {
            this->update_Tp( T, p );
            const real & tV = mHelmholtzVals[ BELFEM_HELMHOLTZ_V ];

            return this->dpdT( T, tV ) / p ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::kappa( const real T, const real p ) const
        {
            this->update_Tp( T, p );
            const real & tV = mHelmholtzVals[ BELFEM_HELMHOLTZ_V ];
            return -1.0 / ( tV * this->dpdv( T, tV ) ) ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::v( const real T, const real p ) const
        {
            /* always-active validity guard: every (T,p) evaluation funnels
             * through here — Helmholtz::update_Tp as well as the transport
             * classes' own update_Tp. Outside the published pressure range
             * the equation of state returns plausible wrong numbers, and if
             * a gas fails, the model fails too. */
            BELFEM_ERROR( p <= mPmax,
                    "Pressure %f MPa exceeds the validity limit %f MPa of the equation of state",
                    ( float ) ( p * 1e-6 ), ( float ) ( mPmax * 1e-6 ) );

            real v ;

            // relaxation factor
            real tOmega ;

            // set the liquid flag of the parent
            mParent.set_liquid_flag( this->is_liquid( T, p ) );

            if( this->is_liquid( T, p ) )
            {
                v = polyval( mVliq, T ) ;
                tOmega = 0.5 ;
            }
            else
            {
                v = mCubicEoS->v( T, p );

                // empirical factor
                tOmega = std::max( 1.0 - 2.0 * std::exp( - T / mTcrit ), 0.1 );
            }

            real tF = BELFEM_REAL_MAX ;
            real tDeltaV ;

            uint tCount = 0;

            while ( std::abs( tF / p ) > 1e-8 )
            {
                tF = this->p( T, v ) - p ;

                tDeltaV = tF / this->dpdv( T, v );

                // limit step
                if ( tDeltaV > v )
                {
                    tDeltaV = v ;
                }

                v -= tOmega * tDeltaV ;

                BELFEM_ERROR( tCount++ < 200,
                             "Too many iterations for T=%8.3f K, p=%8.3f bar, rho=%8.3f kg/m^3, relax=%8.3f",
                             ( double ) T, ( double ) p * 1e-5, ( double ) 1.0 / v, ( double ) tOmega
                             );
            }

            return v ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::T( const real p, const real v ) const
        {
            real T ;

            // check if we might know the solution already
            real tA = std::abs( p - this->get( BELFEM_HELMHOLTZ_P ) ) / p ;
            real tB = std::abs( v - this->get( BELFEM_HELMHOLTZ_V ) ) / v ;

            // use data from memory as initial guess
            if( tA < 0.25 && tB < 0.25 )
            {
                T = this->get( BELFEM_HELMHOLTZ_T );
            }
            else
            {
                T = std::max( p * v / mR, mTtriple );
            }

            real tF = BELFEM_REAL_MAX ;
            real tdF ;
            real tdT ;
            uint tCount = 0 ;
            real tOmega = 0.99 ;

            while( std::abs( tF / p ) > 1e-8 )
            {

                tF  = this->p( T, v ) - p ;
                tdF = this->dpdT( T, v ) ;

                tdT = tF / tdF ;
                if( tdT > ( T - mTtriple ) )
                {
                    tdT = T - mTtriple ;
                    tOmega = 0.3 ;
                }

                T -= tOmega * tdT ;

                BELFEM_ERROR( tCount++ < 100,
                             "Iteration failure in T(p,v) for p=%8.3f Pa v=%10.4f m^3/kg",
                             p, v );
            }

            return T ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::u( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            return mTau * ( this->phi0_t() + this->phir_t() )  * mR * T + mU0 ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::h( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            return ( 1.0 + mTau * ( this->phi0_t() + this->phir_t() )
                     + mDelta * this->phir_d() ) * mR * T  + mH0 ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::s( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            return ( mTau * ( this->phi0_t() + this->phir_t() )
                  - this->phi0()  - this->phir() ) * mR + mS0 ;

        }

//----------------------------------------------------------------------------

        real
        Helmholtz::cv( const real T, const real p ) const
        {
            this->update_Tp( T, p );
            return - mTau * mTau * ( this->phi0_tt() + this->phir_tt() ) * mR ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::cp( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            return ( std::pow( ( 1.0 + mDelta * ( this->phir_d() - mTau * this->phir_dt() ) ), 2 ) /
                 ( 1.0 + mDelta * ( 2.0 * this->phir_d() + mDelta * this->phir_dd()))
                   - mTau * mTau * ( this->phi0_tt() + this->phir_tt()  ) ) * mR ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::w( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            real tA = 1.0 + mDelta * ( 2.0 * this->phir_d() + mDelta * this->phir_dd());
            real tB = 1.0 + mDelta * ( this->phir_d() - mTau * this->phir_dt());
            real tC = mTau * mTau * ( this->phi0_tt() + this->phir_tt());

            return std::sqrt(( tA - tB * tB / tC ) * ( mR * T ));
        }
//----------------------------------------------------------------------------

        real
        Helmholtz::dsdT( const real T, const real p ) const
        {
            return this->cp( T, p ) / T ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dsdp( const real T, const real p ) const
        {
            this->update_Tp( T, p );

            return - this->dvdT( T, mHelmholtzVals[ BELFEM_HELMHOLTZ_V ] );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phi0() const
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
        Helmholtz::phi0_t() const
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
        Helmholtz::phi0_tt() const
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
        Helmholtz::phir() const
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR, this->compute_phir() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_d() const
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_D ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_D, this->compute_phir_d() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_D );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_dd() const
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_DD ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_DD, this->compute_phir_dd() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_DD );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_t() const
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_T ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_T, this->compute_phir_t() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_T );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_tt() const
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_TT ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_TT, this->compute_phir_tt() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_TT );
        }

//----------------------------------------------------------------------------

        const real &
        Helmholtz::phir_dt() const
        {
            if( ! this->test( BELFEM_HELMHOLTZ_PHIR_DT ) )
            {
                this->set( BELFEM_HELMHOLTZ_PHIR_DT, this->compute_phir_dt() ) ;
            }
            return this->get( BELFEM_HELMHOLTZ_PHIR_DT );
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phi0() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phi0() not implemented for this Helmholtz model" );
            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir() not implemented for this Helmholtz model" );
            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phi0_t() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phi0_t() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phi0_tt() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phi0_tt() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_d() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_d() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_dd() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_dd() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_t() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_t() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_tt() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_tt() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::compute_phir_dt() const
        {
            BELFEM_ERROR(
                    false,
                    "compute_phir_dt() not implemented for this Helmholtz model" );

            return BELFEM_QUIET_NAN ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::p_vap( const real T ) const
        {
            return std::exp( this->pi_vap( T ) ) * mPcrit ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::T_vap( const real p ) const
        {
            // compute initial guess
            real T = polyval( mTvap, std::log( p / mPcrit ) );

            real tDeltaT = BELFEM_REAL_MAX ;
            real tPvap ;
            real tPivap ;
            uint tCount = 0 ;
            while( std::abs( tDeltaT ) > BELFEM_EPSILON_T )
            {
                tPivap = this->pi_vap( T );
                tPvap = std::exp( tPivap ) * mPcrit ;
                tDeltaT = ( tPvap - p ) / this->dpvap_dT( T, tPvap, tPivap ) ;

                T -= 0.99 * tDeltaT ;

                BELFEM_ERROR( tCount++ < 100,
                        "T_vap did not converge for %s at p = %e Pa",
                        mLabel.c_str(), ( double ) p );
            }

            return T ;
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::hvap( const real T, const real p ) const
        {
            if( T < mTcrit )
            {
                real tTvap = this->T_vap( p );
                real aHliq = this->h( tTvap - BELFEM_EPSILON_T, p );
                real aHgas = this->h( tTvap + BELFEM_EPSILON_T, p );
                return aHgas - aHliq ;
            }
            else
            {
                return 0.0 ;
            }
        }

//----------------------------------------------------------------------------

        real
        Helmholtz::dpvap_dT( const real T, const real aPvap, const real aPiVap ) const
        {
            real tdG = - ( aPiVap + this->psi_vap( T ) )/ T;

            return ( T > mTcrit ) ? -aPvap * tdG : aPvap * tdG ;
        }

//----------------------------------------------------------------------------

        bool
        Helmholtz::is_liquid( const real T, const real p ) const
        {
            if( T < mTtriple )
            {
                return true ;
            }
            else if ( T > mTcrit )
            {
                return false ;
            }
            else
            {
                return p > this->p_vap( T );
            }
        }

//----------------------------------------------------------------------------
    }
}