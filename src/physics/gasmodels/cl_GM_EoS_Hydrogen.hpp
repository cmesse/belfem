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

#ifndef BELFEM_CL_GM_HELMHOLTZ_HYDROGEN_HPP
#define BELFEM_CL_GM_HELMHOLTZ_HYDROGEN_HPP
#include "en_Helmholtz.hpp"
#include "cl_GM_Helmholtz.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gasmodels
    {
        /**
         * Fundamental Helmholtz equation of state for hydrogen, after
         * Leachman, Jacobsen, Penoncello and Lemmon 2009,
         * doi 10.1063/1.3160306.
         *
         * Covers all three spin variants - parahydrogen, normal hydrogen and
         * orthohydrogen - selected by the HelmholtzModel passed to the
         * constructor, which swaps the whole coefficient set including the
         * critical point and the vapor pressure ancillary. Normal hydrogen is
         * the one that carries CAS 1333-74-0 and is what the gas tables resolve
         * by default.
         *
         * Valid from the triple point to 1000 K and to 2000 MPa.
         *
         * Unlike its methane and oxygen siblings the tau exponents of this
         * formulation do not sit on a half integer grid, so the powers are
         * evaluated with std::pow rather than the hpow helper of the base.
         *
         * @ingroup grp_physics_gasmodels
         * @see @ref physics_gasmodels_gasmodels_usage_guide
         */
        class EoS_Hydrogen : public Helmholtz
        {
            //! parameter a from Table 4, needed for alpha0, eq. (31)
            Vector< real > mA ;

            //! parameter b from Table 4, needed for alpha0, eq. (31)
            Vector< real > mB ;

            //! number of a/b term pairs from Table 4, needed for alpha0, eq. (31)
            real mNab ;

            //! parameter N from Table 5, needed for alpha_r, eq. (32)
            Vector< real > mN ;

            //! parameter t from Table 5, needed for alpha_r, eq. (32)
            Vector< real > mT ;

            //! parameter d from Table 5, needed for alpha_r, eq. (32)
            //! delta exponents; integer, so the powers are built by repeated
            //! multiplication rather than std::pow ( see update_delta_pow_d )
            Vector< uint > mD ;

            //! parameter phi from Table 6, needed for alpha_r, eq. (32)
            Vector< real > mAlpha ;

            //! parameter beta from Table 6, needed for alpha_r, eq. (32)
            Vector< real > mBeta ;

            //! parameter gamma from Table 6, needed for alpha_r, eq. (32)
            Vector< real > mGamma ;

            //! parameter D from Table 6, needed for alpha_r, eq. (32)
            Vector< real > mPsi ;

            //! parameter N from Table 8, needed for vapor pressure
            //Vector< real > mVapN ;

            //! parameter k from Table 8, needed for vapor pressure
            //Vector< real > mVapK ;


            // container for expression delta^d
            // entry 14 is the stamp and contains the value of delta
            mutable Vector< real > mDeltaPowD ;

            // container for expression tau^theta
            // entry 14 is the stamp and contains the value of tau
            mutable Vector< real > mTauPowT ;

            // help function for phi0
            mutable Vector< real > mE ;

            // help function
            // entry 14 contains value of tau
            // entry 15 contains value of delta
            mutable Vector< real > mF ;

            // help function
            // entry 14 contains value of tau
            // entry 15 contains value of delta
            mutable Vector< real > mG ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            EoS_Hydrogen( Gas & aParent, const HelmholtzModel aModel );

            ~EoS_Hydrogen();

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            /**
             * sets the lookup table, called by constructor
             * @param aModel
             */
            void
            select_table( const HelmholtzModel aModel  );

//----------------------------------------------------------------------------
// help functions
//----------------------------------------------------------------------------

            void
            update_e() const ;

            void
            update_f() const ;

//----------------------------------------------------------------------------

            void
            update_g() const;

//----------------------------------------------------------------------------

            void
            update_delta_pow_d() const;

//----------------------------------------------------------------------------

            void
            update_tau_pow_theta() const;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            real
            compute_phi0() const ;

 //----------------------------------------------------------------------------

            real
            compute_phir() const ;

//----------------------------------------------------------------------------

            real
            compute_phi0_t() const ;

//----------------------------------------------------------------------------

            real
            compute_phi0_tt() const ;

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
// help functions for real gas contribution to helmholtz function
//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_e() const
        {
            // the e-values depend on tau only; mE( 0 ) is the stamp
            if ( mTau != mE( 0 ) )
            {
                mE( 0 ) = mTau ;
                for( uint k=2; k<mNab; ++k )
                {
                    mE( k ) = std::exp( mB( k ) * mTau );
                }
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_f() const
        {
            this->update_tau_pow_theta() ;
            this->update_delta_pow_d() ;

            if( mTau != mF( 14 ) || mDelta != mF( 15 ) )
            {
                mF(  0 ) = mN(  0 ) * mDeltaPowD(  0 ) * mTauPowT(  0 ) ;
                mF(  1 ) = mN(  1 ) * mDeltaPowD(  1 ) * mTauPowT(  1 ) ;
                mF(  2 ) = mN(  2 ) * mDeltaPowD(  2 ) * mTauPowT(  2 ) ;
                mF(  3 ) = mN(  3 ) * mDeltaPowD(  3 ) * mTauPowT(  3 ) ;
                mF(  4 ) = mN(  4 ) * mDeltaPowD(  4 ) * mTauPowT(  4 ) ;
                mF(  5 ) = mN(  5 ) * mDeltaPowD(  5 ) * mTauPowT(  5 ) ;
                mF(  6 ) = mN(  6 ) * mDeltaPowD(  6 ) * mTauPowT(  6 ) ;
                mF(  7 ) = mN(  7 ) * mDeltaPowD(  7 ) * mTauPowT(  7 ) ;
                mF(  8 ) = mN(  8 ) * mDeltaPowD(  8 ) * mTauPowT(  8 ) ;
                mF(  9 ) = mN(  9 ) * mDeltaPowD(  9 ) * mTauPowT(  9 ) ;
                mF( 10 ) = mN( 10 ) * mDeltaPowD( 10 ) * mTauPowT( 10 ) ;
                mF( 11 ) = mN( 11 ) * mDeltaPowD( 11 ) * mTauPowT( 11 ) ;
                mF( 12 ) = mN( 12 ) * mDeltaPowD( 12 ) * mTauPowT( 12 ) ;
                mF( 13 ) = mN( 13 ) * mDeltaPowD( 13 ) * mTauPowT( 13 ) ;

                mF( 14 ) = mTau ;
                mF( 15 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_g() const
        {
            if( mTau != mG( 14 ) || mDelta != mG( 15 ) )
            {
                mG( 7 ) = std::exp( -mDelta );
                mG( 8 ) = mG( 7 );

                mG(  9 ) = std::exp(
                        mAlpha(  9 ) * std::pow( mDelta - mPsi(  9 ), 2 )
                        +   mBeta(  9 ) * std::pow( mTau - mGamma(  9 ), 2 ) );

                mG( 10 ) = std::exp(
                        mAlpha( 10 ) * std::pow( mDelta - mPsi( 10 ), 2 )
                        +   mBeta( 10 ) * std::pow( mTau - mGamma( 10 ), 2 ) );

                mG( 11 ) = std::exp(
                        mAlpha( 11 ) * std::pow( mDelta - mPsi( 11 ), 2 )
                        +   mBeta( 11 ) * std::pow( mTau - mGamma( 11 ), 2 ) );

                mG( 12 ) = std::exp(
                        mAlpha( 12 ) * std::pow( mDelta - mPsi( 12 ), 2 )
                        +   mBeta( 12 ) * std::pow( mTau - mGamma( 12 ), 2 ) );


                mG( 13 ) = std::exp(
                        mAlpha( 13 ) * std::pow( mDelta - mPsi( 13 ), 2 )
                        +   mBeta( 13 ) * std::pow( mTau - mGamma( 13 ), 2 ) );


                mG( 14 ) = mTau ;
                mG( 15 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_delta_pow_d() const
        {
            if( mDelta != mDeltaPowD( 14 ) )
            {
                // integer exponents: repeated multiplication instead of std::pow
                for( uint k=0; k<14; ++k )
                {
                    real tPow = 1.0 ;

                    for( uint i=0; i<mD( k ); ++i )
                    {
                        tPow *= mDelta ;
                    }

                    mDeltaPowD( k ) = tPow ;
                }

                mDeltaPowD( 14 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_tau_pow_theta() const
        {
            if( mTau != mTauPowT( 14 ) )
            {
                mTauPowT(  0 ) = std::pow( mTau, mT(  0 ) );
                mTauPowT(  1 ) = std::pow( mTau, mT(  1 ) );
                mTauPowT(  2 ) = std::pow( mTau, mT(  2 ) );
                mTauPowT(  3 ) = std::pow( mTau, mT(  3 ) );
                mTauPowT(  4 ) = std::pow( mTau, mT(  4 ) );
                mTauPowT(  5 ) = std::pow( mTau, mT(  5 ) );
                mTauPowT(  6 ) = std::pow( mTau, mT(  6 ) );
                mTauPowT(  7 ) = std::pow( mTau, mT(  7 ) );
                mTauPowT(  8 ) = std::pow( mTau, mT(  8 ) );
                mTauPowT(  9 ) = std::pow( mTau, mT(  9 ) );
                mTauPowT( 10 ) = std::pow( mTau, mT( 10 ) );
                mTauPowT( 11 ) = std::pow( mTau, mT( 11 ) );
                mTauPowT( 12 ) = std::pow( mTau, mT( 12 ) );
                mTauPowT( 13 ) = std::pow( mTau, mT( 13 ) );

                mTauPowT( 14 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HELMHOLTZ_HYDROGEN_HPP
