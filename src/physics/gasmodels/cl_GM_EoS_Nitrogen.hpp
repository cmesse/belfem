/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_GM_EOS_NITROGEN_HPP
#define BELFEM_CL_GM_EOS_NITROGEN_HPP

#include "en_Helmholtz.hpp"
#include "cl_GM_Helmholtz.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gasmodels
    {

        /**
         * 10.1063/1.1349047
         */
        class EoS_Nitrogen : public Helmholtz
        {
            const Vector< real > mA ;
            const Vector< real > mN ;
            const Vector< uint > mI ;
            const Vector< real > mJ ;
            const Vector< uint > mL ;
            const Vector< real > mBeta ;
            const Vector< real > mGamma ;
            const Vector< real > mPhi ;

            mutable Vector< real > mE ;
            mutable Vector< real > mF ;
            mutable Vector< real > mDeltaPowI ;
            mutable Vector< real > mDeltaPowL ;
            mutable Vector< real > mTauPowJ ;
            mutable Vector< real > mSwap ;
            mutable real mChi ; // <-- help parameter
        public:

            EoS_Nitrogen( Gas & aParent ) ;

            ~EoS_Nitrogen() ;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            real
            compute_phi0() const override ;

            real
            compute_phi0_t() const override;


            real
            compute_phi0_tt() const override;


            real
            compute_phir() const override ;


            real
            compute_phir_d() const override ;

            real
            compute_phir_dd() const override ;

            real
            compute_phir_t() const override ;

            real
            compute_phir_tt() const override ;

            real
            compute_phir_dt() const override ;

        private:

            /**
             * critical point, validity range, vapor pressure ancillary and
             * the initial guess polynomials of the base class
             */
            void
            init_tables() ;

            void
            update_e() const;

            void
            update_f() const;

            void
            update_delta_pow_i() const;

            void
            update_tau_pow_j() const;

        };

        inline void
        EoS_Nitrogen::update_e() const
        {
            if ( mE( 2 ) == mTau ) return ;

            mChi = std::exp( -mA( 7 ) * mTau ) ;

            mE( 0 ) = std::log( mTau );
            mE( 1 ) = 1. ;
            mE( 2 ) = mTau ;
            mE( 3 ) = 1./mTau ;
            mE( 4 ) = mE( 3 )*mE( 3 );
            mE( 5 ) = mE( 3 )*mE( 4 );
            mE( 6 ) = std::log( 1.- mChi );
            mE( 7 ) = 0. ;
        }

        inline void
        EoS_Nitrogen::update_f() const
        {
            this->update_delta_pow_i();
            this->update_tau_pow_j();

            bool tDelta = mF( 36 ) == mDelta ;
            bool tTau = mF( 37 ) == mTau ;

            if ( tDelta && tTau ) return ;

            if ( ! tDelta )
            {
                for ( uint k=6; k<32; ++k )
                {
                    uint n = mL( k );
                    mDeltaPowL( k ) = mDelta;
                    for ( uint l=1; l<n; ++l )
                    {
                        mDeltaPowL( k ) *= mDelta;
                    }
                    mF( k ) = std::exp( -mDeltaPowL( k ) );
                }
            }

            uint k = 32 ;
            real xi = mDelta - 1. ;
            xi *= xi ;

            for ( uint m=0; m<4; ++m, ++k )
            {
                real eta = mTau - mGamma( m );
                mF( k ) = std::exp( -mPhi( m ) * xi - mBeta( m ) * eta * eta );
            }

            mF( 36 ) = mDelta;
            mF( 37 ) = mTau;
        }
        inline void
        EoS_Nitrogen::update_delta_pow_i() const
        {
            if ( mDelta != mDeltaPowI( 36 ) )
            {

                mDeltaPowI.fill( 1. );

                for ( uint k=0; k<36; ++k )
                {
                    uint n = mI( k );
                    BELFEM_ASSERT( n > 0, "Invalid exponent" );
                    mDeltaPowI( k ) = mDelta;
                    for ( uint i=1; i<n; ++i )
                    {
                        mDeltaPowI( k ) *= mDelta;
                    }
                }
                mDeltaPowI( 36 ) = mDelta;
            }
        }

        inline void
        EoS_Nitrogen::update_tau_pow_j() const
        {
            if ( mTau != mTauPowJ( 36 ) )
            {
                for ( uint k=0; k<36; ++k )
                {
                    mTauPowJ( k ) = std::pow( mTau, mJ( k ) );
                }
                mTauPowJ( 36 ) = mTau;
            }
        }


    }
}

#endif //BELFEM_CL_GM_EOS_NITROGEN_HPP
