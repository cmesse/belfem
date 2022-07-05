//
// Created by Christian Messe on 18.08.20.
//

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
         */
        class EoS_Methane : public Helmholtz
        {
            // ideal gas table
            Vector< real > mA ;
            Vector< real > mB ;  // theta in Paper

            // real gas table
            Vector< real > mN ;
            Vector< real > mT ;
            Vector< real > mD ;
            Vector< real > mAlpha ;
            Vector< real > mBeta ;
            Vector< real > mPsi ;
            Vector< real > mGamma ;
            Vector< real > mC ;

            // containers
            Vector< real > mDeltaPowD ;
            Vector< real > mTauPowT ;
            Vector< real > mE ;
            Vector< real > mF ;
            Vector< real > mG ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

             EoS_Methane( Gas & aParent ) ;

            ~EoS_Methane();

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
        private:
//----------------------------------------------------------------------------

            void
            init_tables() ;

//----------------------------------------------------------------------------
// help functions
//----------------------------------------------------------------------------

            void
            update_e() ;

//----------------------------------------------------------------------------
            void
            update_f() ;

//----------------------------------------------------------------------------

            void
            update_g();

//----------------------------------------------------------------------------

            void
            update_delta_pow_d();

//----------------------------------------------------------------------------

            void
            update_tau_pow_theta();

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
// help functions for real gas contribution to helmholtz function
//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_e()
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
        EoS_Methane::update_f()
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
        EoS_Methane::update_g()
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
                mG( 41 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_delta_pow_d()
        {
            if( mDelta != mDeltaPowD( 40 ) )
            {
                for( uint k=0; k<40; ++k )
                {
                    mDeltaPowD( k ) = std::pow( mDelta, mD( k ) );
                }

                mDeltaPowD( 40 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Methane::update_tau_pow_theta()
        {
            if( mTau != mTauPowT( 40 ) )
            {
                for( uint k=0; k<40; ++k )
                {
                    mTauPowT( k ) = std::pow( mTau, mT( k ) );
                }

                mTauPowT( 40 ) = mTau ;
            }
        }

//----------------------------------------------------------------------------

    }
}

#endif //BELFEM_CL_GM_HELMHOLTZ_METHANE_HPP
