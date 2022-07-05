//
// Created by Christian Messe on 18.08.20.
//

#ifndef BELFEM_CL_GM_EOS_OXYGEN_VAPOR_HPP

#include "en_Helmholtz.hpp"
#include "cl_GM_Helmholtz.hpp"
#include "cl_Vector.hpp"
#include "constants.hpp"

namespace belfem
{
    namespace gasmodels
    {
        /**
         * see https://doi.org/10.1016/0378-3812(85)87016-3
         */
        class EoS_Oxygen : public Helmholtz
        {

            // reference conditions for enthalpy and entropy
            const real mT0     = BELFEM_TREF ;
            const real mP0     = BELFEM_PREF ;
                  real mTau0   = 0.0 ;
                  real mDelta0 = 0.0 ;

            // coefficient for ideal gas contribution
            Vector< real > mK ;

            // coefficients for real gas contribution
            Vector< real > mN ;
            Vector< real > mD ; // actually r in paper
            Vector< real > mT ; // actually s in paper

            // container with evaluated values for real gas contribution
            Vector< real > mDeltaPowD ;
            Vector< real > mTauPowT ;

            // help function for idgas contribution
            Vector< real > mE ;

            // help function delta^A, entry 24 contains delta value
            Vector< real > mF ;

            // help function tau^B, entry 2 contains tau value
            Vector< real > mG ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

             EoS_Oxygen( Gas & aParent ) ;

            ~EoS_Oxygen();

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            real
            pi_vap( const real & aT );

//----------------------------------------------------------------------------

            real
            psi_vap( const real & aT );

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            init_tables() ;

            void
            update_delta_pow_a() ;

            void
            update_tau_pow_b() ;

            void
            update_e() ;

            void
            update_f() ;

            void
            update_g() ;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            real
            compute_phi0() ;

//----------------------------------------------------------------------------

            real
            compute_phi0_t() ;

//----------------------------------------------------------------------------

            real
            compute_phi0_tt() ;

//----------------------------------------------------------------------------

            real
            compute_phir() ;

//----------------------------------------------------------------------------

            real
            compute_phir_d() ;

//----------------------------------------------------------------------------

            real
            compute_phir_dd() ;

//----------------------------------------------------------------------------

            real
            compute_phir_t() ;

//----------------------------------------------------------------------------

            real
            compute_phir_tt() ;

//----------------------------------------------------------------------------

            real
            compute_phir_dt() ;

//----------------------------------------------------------------------------

        };

//----------------------------------------------------------------------------

        inline real
        EoS_Oxygen::pi_vap( const real & aT )
        {
            real tTheta =  1.0 - aT / mTcrit ;

            return ( (   mNvap( 0 )   * tTheta
                       + mNvap( 1 ) ) * tTheta
                       + mNvap( 2 ) ) * tTheta * mTcrit / aT ;
        }

//----------------------------------------------------------------------------

        inline real
        EoS_Oxygen::psi_vap( const real & aT )
        {
            real tTheta =  1.0 - aT / mTcrit ;

            return      (   3.0 * mNvap( 0 )   * tTheta
                          + 2.0 * mNvap( 1 ) ) * tTheta
                          +       mNvap( 2 ) ;
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_delta_pow_a()
        {
            if( mDelta != mDeltaPowD( 32 ) )
            {
                for( uint k=0; k<32; ++k )
                {
                    mDeltaPowD( k ) = std::pow( mDelta, mD( k ) );
                }
                mDeltaPowD( 32 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_tau_pow_b()
        {
            if( mTau != mTauPowT( 32 ) )
            {
                for( uint k=0; k<32; ++k )
                {
                    mTauPowT( k ) = std::pow( mTau, mT( k ) );
                }
                mTauPowT( 32 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_e()
        {
            if( mTau != mE( 2 ) )
            {

                mE( 0 ) = std::sqrt( mTau );
                mE( 1 ) = std::exp(  mK( 6 ) * mTau );
                mE( 2 ) = 2./3. * std::exp( -mK( 7 ) * mTau );
                mE( 3 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_f()
        {
            this->update_delta_pow_a() ;
            this->update_tau_pow_b() ;

            if( mDelta != mF( 32 ) || mTau != mF( 33 ) )
            {
                for( uint k=0; k<32; ++k )
                {
                    mF( k ) = mN( k ) * mDeltaPowD( k ) * mTauPowT( k );
                }
                mF( 32 ) = mDelta ;
                mF( 33 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_g()
        {
            if( mDelta != mG( 2 ) )
            {
                mG( 0 ) = std::exp( -mDelta * mDelta );
                mG( 1 ) = std::exp( -mDelta * mDelta * mDelta * mDelta );
                mG( 2 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_EOS_OXYGEN_VAPOR_HPP
