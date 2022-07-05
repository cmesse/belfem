//
// Created by Christian Messe on 17.08.20.
//

#ifndef BELFEM_CL_GM_HELMHOLTZ_HYDROGEN_HPP
#define BELFEM_CL_GM_HELMHOLTZ_HYDROGEN_HPP
#include "en_Helmholtz.hpp"
#include "cl_GM_Helmholtz.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gasmodels
    {
        class EoS_Hydrogen : public Helmholtz
        {
            //! parameter a from Table 4, needed for alpha0, eq. (31)
            Vector< real > mA ;

            //! parameter b from Table 4, needed for alpha0, eq. (31)
            Vector< real > mB ;

            //! parameter n from from Table 4, needed for alpha0, eq. (31)
            real mNab ;

            //! parameter N from Table 5, needed for alpha_r, eq. (32)
            Vector< real > mN ;

            //! parameter t from Table 5, needed for alpha_r, eq. (32)
            Vector< real > mT ;

            //! parameter d from Table 5, needed for alpha_r, eq. (32)
            Vector< real > mD ;

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
            // entry 14 contains value of delta
            Vector< real > mDeltaPowD ;

            // container for expression tau^theta
            // entry 14 contains value of delta
            Vector< real > mTauPowT ;

            // help funciton for phi0
            Vector< real > mE ;

            // help function
            // entry 14 contains value of tau
            // entry 15 contains value of delta
            Vector< real > mF ;

            // help function
            // entry 14 contains value of tau
            // entry 15 contains value of delta
            Vector< real > mG ;

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
            update_e() ;

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
        protected:
//----------------------------------------------------------------------------

            real
            compute_phi0() ;

 //----------------------------------------------------------------------------

            real
            compute_phir() ;

//----------------------------------------------------------------------------

            real
            compute_phi0_t() ;

//----------------------------------------------------------------------------

            real
            compute_phi0_tt() ;

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
// help functions for real gas contribution to helmholtz function
//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_e()
        {
            if ( mTau != mF( 14 ) || mDelta != mE( 0 ) )
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
        EoS_Hydrogen::update_f()
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
        EoS_Hydrogen::update_g()
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
        EoS_Hydrogen::update_delta_pow_d()
        {
            if( mDelta != mDeltaPowD( 14 ) )
            {
                mDeltaPowD(  0 ) = std::pow( mDelta, mD(  0 ) );
                mDeltaPowD(  1 ) = std::pow( mDelta, mD(  1 ) );
                mDeltaPowD(  2 ) = std::pow( mDelta, mD(  2 ) );
                mDeltaPowD(  3 ) = std::pow( mDelta, mD(  3 ) );
                mDeltaPowD(  4 ) = std::pow( mDelta, mD(  4 ) );
                mDeltaPowD(  5 ) = std::pow( mDelta, mD(  5 ) );
                mDeltaPowD(  6 ) = std::pow( mDelta, mD(  6 ) );
                mDeltaPowD(  7 ) = std::pow( mDelta, mD(  7 ) );
                mDeltaPowD(  8 ) = std::pow( mDelta, mD(  8 ) );
                mDeltaPowD(  9 ) = std::pow( mDelta, mD(  9 ) );
                mDeltaPowD( 10 ) = std::pow( mDelta, mD( 10 ) );
                mDeltaPowD( 11 ) = std::pow( mDelta, mD( 11 ) );
                mDeltaPowD( 12 ) = std::pow( mDelta, mD( 12 ) );
                mDeltaPowD( 13 ) = std::pow( mDelta, mD( 13 ) );

                mDeltaPowD( 14 ) = mDelta ;
            }
        }

//----------------------------------------------------------------------------

        inline void
        EoS_Hydrogen::update_tau_pow_theta()
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
