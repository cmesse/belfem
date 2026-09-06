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

#ifndef BELFEM_CL_GM_HELMHOLTZ_METHANE_HPP
#define BELFEM_CL_GM_HELMHOLTZ_METHANE_HPP

#include "en_Helmholtz.hpp"
#include "cl_GM_Helmholtz.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gasmodels
    {
        /**
         * see https://doi.org/10.1063/1.555898
         *
         * @ingroup grp_physics_gasmodels
         * @see @ref physics_gasmodels_gasmodels_usage_guide
         */
        class EoS_Methane : public Helmholtz
        {
            // ideal gas table
            Vector< real > mA ;
            Vector< real > mB ;  // theta in Paper

            // real gas table
            Vector< real > mN ;
            Vector< real > mT ;
            //! delta exponents; integer, so the powers are built by repeated
            //! multiplication rather than std::pow ( see update_delta_pow_d )
            Vector< uint > mD ;
            Vector< real > mAlpha ;
            Vector< real > mBeta ;
            Vector< real > mPsi ;
            Vector< real > mGamma ;

            //! exponent of delta inside the exponential of the residual terms;
            //! integer, see update_delta_pow_d
            Vector< uint > mC ;

            // containers
            mutable Vector< real > mDeltaPowD ;

            //! delta^mC, cached alongside mDeltaPowD; the derivative loops used
            //! to call std::pow for this on every evaluation
            mutable Vector< real > mDeltaPowC ;

            mutable Vector< real > mTauPowT ;
            mutable Vector< real > mE ;
            mutable Vector< real > mF ;
            mutable Vector< real > mG ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

             EoS_Methane( Gas & aParent ) ;

            ~EoS_Methane();

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
        private:
//----------------------------------------------------------------------------

            void
            init_tables() ;

//----------------------------------------------------------------------------
// help functions
//----------------------------------------------------------------------------

            void
            update_e() const ;

//----------------------------------------------------------------------------
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
        };

//----------------------------------------------------------------------------
// help functions for real gas contribution to helmholtz function
//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_e() const
        {
            if( mTau != mE( 0 ) )
            {
                mE( 0 ) = mTau ;
                mE( 3 ) = std::exp( mB( 3 ) * mTau ) ;
                mE( 4 ) = std::exp( mB( 4 ) * mTau ) ;
                mE( 5 ) = std::exp( mB( 5 ) * mTau ) ;
                mE( 6 ) = std::exp( mB( 6 ) * mTau ) ;
                mE( 7 ) = std::exp( mB( 7 ) * mTau ) ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_f() const
        {
            this->update_tau_pow_theta() ;
            this->update_delta_pow_d() ;

            if( mTau != mF( 40 ) || mDelta != mF( 41 ) )
            {
                for( uint k=0; k<40; ++k )
                {
                    mF( k ) = mN( k ) * mDeltaPowD( k ) * mTauPowT( k ) ;
                }

                mF( 40 ) = mTau ;
                mF( 41 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_g() const
        {
            if( mTau != mG( 40 ) || mDelta != mG( 41 ) )
            {
                if ( mDelta != mG( 41 ) )
                {
                    mG( 13 ) = std::exp( -mDelta );
                    mG( 14 ) = mG( 13 );
                    mG( 15 ) = mG( 13 );
                    mG( 16 ) = mG( 13 );
                    mG( 17 ) = mG( 13 );
                    mG( 18 ) = mG( 13 );
                    mG( 19 ) = mG( 13 );

                    mG( 20 ) = std::exp( -mDelta * mDelta );
                    mG( 21 ) = mG( 20 );
                    mG( 22 ) = mG( 20 );
                    mG( 23 ) = mG( 20 );
                    mG( 24 ) = mG( 20 );

                    mG( 25 ) = std::exp( -mDelta * mDelta * mDelta );
                    mG( 26 ) = mG( 25 );
                    mG( 27 ) = mG( 25 );
                    mG( 28 ) = mG( 25 );

                    mG( 29 ) = std::exp( -mDelta * mDelta * mDelta * mDelta );
                    mG( 30 ) = mG( 29 );
                    mG( 31 ) = mG( 29 );
                    mG( 32 ) = mG( 29 );
                    mG( 33 ) = mG( 29 );
                    mG( 34 ) = mG( 29 );
                    mG( 35 ) = mG( 29 );
                }

                for( uint k=36; k<40; ++k )
                {
                    mG( k ) = std::exp(
                              mAlpha( k )  * std::pow( mDelta - mPsi( k ) , 2 )
                            + mBeta( k ) * std::pow( mTau - mGamma( k ),  2 ) );
                }

                mG( 40 ) = mTau ;
                mG( 41 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_delta_pow_d() const
        {
            if( mDelta != mDeltaPowD( 40 ) )
            {
                /* the exponents are integers, so repeated multiplication
                 * replaces std::pow( real, real ). Note that mD is zero for the
                 * last three terms, hence the accumulator starts at one. */
                for( uint k=0; k<40; ++k )
                {
                    real tPow = 1.0 ;

                    for( uint i=0; i<mD( k ); ++i )
                    {
                        tPow *= mDelta ;
                    }

                    mDeltaPowD( k ) = tPow ;

                    // delta^c for the exponential of the residual terms
                    mDeltaPowC( k ) = ipow( mDelta, mC( k ) ) ;
                }

                mDeltaPowD( 40 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_tau_pow_theta() const
        {
            if( mTau != mTauPowT( 40 ) )
            {
                /* every tabulated tau exponent is a multiple of one half, so
                 * 2*mT is an integer and hpow() replaces std::pow. The doubled
                 * exponent is derived from mT rather than tabulated a second
                 * time, so the two can never disagree. */
                const real tSqrtTau = std::sqrt( mTau ) ;

                for( uint k=0; k<40; ++k )
                {
                    const int tN = static_cast< int >( std::lround( 2.0 * mT( k ) ) );

                    BELFEM_ASSERT( std::abs( 2.0 * mT( k ) - tN ) < 1e-9,
                            "tau exponent %f of term %u is not a multiple of one half",
                            ( float ) mT( k ), ( unsigned int ) k );

                    mTauPowT( k ) = hpow( mTau, tSqrtTau, tN );
                }

                mTauPowT( 40 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

    }
}

#endif //BELFEM_CL_GM_HELMHOLTZ_METHANE_HPP
