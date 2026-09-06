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

#ifndef BELFEM_CL_GM_EOS_OXYGEN_VAPOR_HPP
#define BELFEM_CL_GM_EOS_OXYGEN_VAPOR_HPP

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
                  real mDelta0 = 0.0 ;

            // coefficient for ideal gas contribution
            Vector< real > mK ;

            // coefficients for real gas contribution
            Vector< real > mN ;
            //! actually r in paper; integer, so the powers are built by repeated
            //! multiplication rather than std::pow ( see update_delta_pow_a )
            Vector< uint > mD ;
            Vector< real > mT ; // actually s in paper

            // container with evaluated values for real gas contribution
            mutable Vector< real > mDeltaPowD ;
            mutable Vector< real > mTauPowT ;

            // help function for idgas contribution
            mutable Vector< real > mE ;

            // residual terms N_k delta^d_k tau^t_k, entries 32 and 33 are the delta and tau stamps
            mutable Vector< real > mF ;

            // exp( -delta^2 ) and exp( -delta^4 ), entry 2 is the delta stamp
            mutable Vector< real > mG ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

             EoS_Oxygen( Gas & aParent ) ;

            ~EoS_Oxygen();

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            real
            pi_vap( const real T ) const;

//----------------------------------------------------------------------------

            real
            psi_vap( const real T ) const;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            init_tables() ;

            void
            update_delta_pow_a() const ;

            void
            update_tau_pow_b() const ;

            void
            update_e() const ;

            void
            update_f() const ;

            void
            update_g() const ;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            real
            compute_phi0() const ;

//----------------------------------------------------------------------------

            real
            compute_phi0_t() const ;

//----------------------------------------------------------------------------

            real
            compute_phi0_tt() const ;

//----------------------------------------------------------------------------

            real
            compute_phir() const ;

//----------------------------------------------------------------------------

            real
            compute_phir_d() const ;

//----------------------------------------------------------------------------

            real
            compute_phir_dd() const ;

//----------------------------------------------------------------------------

            real
            compute_phir_t() const ;

//----------------------------------------------------------------------------

            real
            compute_phir_tt() const ;

//----------------------------------------------------------------------------

            real
            compute_phir_dt() const ;

//----------------------------------------------------------------------------

        };

//----------------------------------------------------------------------------

        inline real
        EoS_Oxygen::pi_vap( const real T ) const
        {
            real tTheta =  1.0 - T / mTcrit ;

            return ( (   mNvap( 0 )   * tTheta
                       + mNvap( 1 ) ) * tTheta
                       + mNvap( 2 ) ) * tTheta * mTcrit / T ;
        }

//----------------------------------------------------------------------------

        inline real
        EoS_Oxygen::psi_vap( const real T ) const
        {
            real tTheta =  1.0 - T / mTcrit ;

            return      (   3.0 * mNvap( 0 )   * tTheta
                          + 2.0 * mNvap( 1 ) ) * tTheta
                          +       mNvap( 2 ) ;
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_delta_pow_a() const
        {
            if( mDelta != mDeltaPowD( 32 ) )
            {
                // integer exponents: repeated multiplication instead of std::pow
                for( uint k=0; k<32; ++k )
                {
                    real tPow = 1.0 ;

                    for( uint i=0; i<mD( k ); ++i )
                    {
                        tPow *= mDelta ;
                    }

                    mDeltaPowD( k ) = tPow ;
                }
                mDeltaPowD( 32 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_tau_pow_b() const
        {
            if( mTau != mTauPowT( 32 ) )
            {
                /* every tabulated tau exponent is a multiple of one half, so
                 * 2*mT is an integer and hpow() replaces std::pow. The doubled
                 * exponent is derived from mT rather than tabulated a second
                 * time, so the two can never disagree. */
                const real tSqrtTau = std::sqrt( mTau ) ;

                for( uint k=0; k<32; ++k )
                {
                    const int tN = static_cast< int >( std::lround( 2.0 * mT( k ) ) );

                    BELFEM_ASSERT( std::abs( 2.0 * mT( k ) - tN ) < 1e-9,
                            "tau exponent %f of term %u is not a multiple of one half",
                            ( float ) mT( k ), ( unsigned int ) k );

                    mTauPowT( k ) = hpow( mTau, tSqrtTau, tN );
                }
                mTauPowT( 32 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_e() const
        {
            // mE( 3 ) is the stamp; mE( 0-2 ) hold values
            if( mTau != mE( 3 ) )
            {

                mE( 0 ) = std::sqrt( mTau );
                mE( 1 ) = std::exp(  mK( 6 ) * mTau );
                mE( 2 ) = 2./3. * std::exp( -mK( 7 ) * mTau );
                mE( 3 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Oxygen::update_f() const
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
        EoS_Oxygen::update_g() const
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
