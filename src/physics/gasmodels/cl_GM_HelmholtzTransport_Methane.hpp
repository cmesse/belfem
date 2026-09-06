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

#ifndef BELFEM_CL_GM_HelmholtzTransportMETHANE_HPP
#define BELFEM_CL_GM_HelmholtzTransportMETHANE_HPP

#include "cl_GM_HelmholtzTransport.hpp"

#include "typedefs.hpp"
#include "cl_Gas.hpp"
#include "cl_Vector.hpp"
#include "cl_Bitset.hpp"

#define BELFEM_METHANE_T                 0
#define BELFEM_METHANE_P                 1
#define BELFEM_METHANE_V                 2
#define BELFEM_METHANE_TAU               3
#define BELFEM_METHANE_DELTA             4
#define BELFEM_METHANE_MU0               5
#define BELFEM_METHANE_MU                6
#define BELFEM_METHANE_LAMBDA            7
#define BELFEM_METHANE_NUMVALS           8
namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * a special class for trans properties of methane
         *
         * 10.1063/1.555828
         */
        class HelmholtzTransport_Methane : public HelmholtzTransport
        {
            // critical data, bound to the EoS data object at runtime.
            // Friend's own Table 1 values are Tc = 190.551 K, pc = 4.5992e6 Pa;
            // running on the EoS critical point instead is deliberate, see
            // init_constants()
            const real mTcrit   ;
            const real mPcrit   ;
            const real mM       ;  // = 16.043e-3 ;
            const real mVcrit   ;  // = 1.0 / ( mM * 10.139e3 ) ;

            const real mR       ; // = constant::Rm / mM ;
            const real mZcrit   ; // = mPcrit * mVcrit / ( mR * mTcrit );

            // Table 8
            Vector< real > mOmega ;
            Vector< real > mF ;

            // Table 9
            Vector< real > mA ;
            Vector< real > mB ;
            Vector< real > mG ;

            Vector< real > mC ;
            Vector< real > mD ;
            Vector< real > mJ ;

            // intermolecular parameters
            const real mEpsilonKb = 174.0;
            const real mSigma = 0.36652e-9 ;

            // constant for eta_0
            real mConstEta0 ;

            // constant for eta_ex
            real mConstEtaEx ;


            real mConstLambdaCr ;
            real mConstExpLambdaCr ;

            real mConstLambdaEx ;

            //! values with stored data
            //! logically const cache, written by the const evaluators
            mutable real mVals[ BELFEM_METHANE_NUMVALS ] = { 0.0 };
            mutable Bitset< BELFEM_METHANE_NUMVALS > mBits;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            HelmholtzTransport_Methane( Gas & aParent );

            ~HelmholtzTransport_Methane() = default ;

//----------------------------------------------------------------------------

            real
            mu( const real T, const real p ) const;

//----------------------------------------------------------------------------

            real
            lambda( const real T, const real p ) const;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            init_tables();

//----------------------------------------------------------------------------

            void
            init_constants();

//----------------------------------------------------------------------------

            real
            eta_0() const;

//----------------------------------------------------------------------------

            real
            eta_ex() const;

//----------------------------------------------------------------------------

            /**
             * Reduced collision integral of Friend Eq. ( 12 ), evaluated at
             * the reduced temperature T / ( epsilon / k ) of the state the
             * caller last set. Takes no arguments; the temperature comes from
             * the internal state cache.
             */
            real
            omega() const;

//----------------------------------------------------------------------------

            real
            lambda_0() const;

//----------------------------------------------------------------------------

            /**
             * lambda near critical point
             * @return
             */
            real
            lambda_cr( const real T, const real p ) const;

//----------------------------------------------------------------------------

            /**
             * lambda excess
             * @return
             */
            real
            lambda_ex() const;

//----------------------------------------------------------------------------

            /**
             * symmetrized compressibility chi_T*, Eq. ( 19 ), for lambda_cr
             * @return
             */
            real
            chi() const;

//----------------------------------------------------------------------------

            /**
             * crossover ( damping ) function F( T*, rho* ), Eq. ( 20 ),
             * for lambda_cr
             * @return
             */
            real
            f() const;


//----------------------------------------------------------------------------

            void
            update_Tp( const real T, const real p ) const;

//----------------------------------------------------------------------------

            bool
            test( const index_t aIndex ) const ;

//----------------------------------------------------------------------------

            void
            set( const index_t aIndex, const real aValue ) const ;

//----------------------------------------------------------------------------

            const real &
            get( const index_t aIndex ) const ;

//----------------------------------------------------------------------------
        } ;

//----------------------------------------------------------------------------

        inline void
        HelmholtzTransport_Methane::update_Tp( const real T, const real p ) const
        {
            BELFEM_ASSERT( T   > 0, "Invalid temprerature" );
            BELFEM_ASSERT( p > 0, "Invalid pressure" );
            if(               T != mVals[ BELFEM_METHANE_T ]
                              || p != mVals[ BELFEM_METHANE_P ] )
            {
                mBits.reset();

                mVals[ BELFEM_METHANE_T ]     = T ;
                mVals[ BELFEM_METHANE_P ]     = p ;
                mVals[ BELFEM_METHANE_V ]     = mEoS.v( T, p );

                mVals[ BELFEM_METHANE_TAU ]   = mTcrit / T ;
                mVals[ BELFEM_METHANE_DELTA ] = mVcrit / mVals[ BELFEM_METHANE_V ] ;


                mBits.set( BELFEM_METHANE_T );
                mBits.set( BELFEM_METHANE_P );
                mBits.set( BELFEM_METHANE_V );
                mBits.set( BELFEM_METHANE_TAU );
                mBits.set( BELFEM_METHANE_DELTA );
            }
        }

//----------------------------------------------------------------------------

        inline bool
        HelmholtzTransport_Methane::test( const index_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_METHANE_NUMVALS ,
                          "Invalid Methane transport state index: %u", ( unsigned int ) aIndex );

            return mBits.test( aIndex );
        }

//----------------------------------------------------------------------------

        inline void
        HelmholtzTransport_Methane::set( const index_t aIndex, const real aValue ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_METHANE_NUMVALS ,
                          "Invalid Methane transport state index: %u", ( unsigned int ) aIndex );

            // set value
            mVals[ aIndex ] = aValue;

            // update flag
            mBits.set( aIndex );
        }

//----------------------------------------------------------------------------

        inline const real &
        HelmholtzTransport_Methane::get( const index_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_METHANE_NUMVALS ,
                          "Invalid Helmholtz state index: %u", ( unsigned int ) aIndex );

            return mVals[ aIndex ];
        }

//----------------------------------------------------------------------------

    }
}

#endif //BELFEM_CL_GM_HelmholtzTransportMETHANE_HPP
