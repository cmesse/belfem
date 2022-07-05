//
// Created by Christian Messe on 23.11.20.
//

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
            // critical data, may be different to the ones in Helmholtz
            const real mTcrit   ;  // = 190.551 ;
            const real mPcrit   ;  // = 4.5992e6 ;
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

            // values with stored data
            real mVals[ BELFEM_METHANE_NUMVALS ] = { 0.0 };
            Bitset< BELFEM_METHANE_NUMVALS > mBits;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            HelmholtzTransport_Methane( Gas & aParent );

            ~HelmholtzTransport_Methane() = default ;

//----------------------------------------------------------------------------

            real
            mu( const real & aT, const real & aP );

//----------------------------------------------------------------------------

            real
            lambda( const real & aT, const real & aP );

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
            eta_0();

//----------------------------------------------------------------------------

            real
            eta_ex();

//----------------------------------------------------------------------------

            /**
             * reduced collision integral
             * @param aT : Temperature in K
             */
            real
            omega();

//----------------------------------------------------------------------------

            real
            lambda_0();

//----------------------------------------------------------------------------

            /**
             * lambda near critical point
             * @return
             */
            real
            lambda_cr( const real & aT, const real & aP );

//----------------------------------------------------------------------------

            /**
             * lambda excess
             * @return
             */
            real
            lambda_ex();

//----------------------------------------------------------------------------

            /**
             * damping function for lambda_cr
             * @return
             */
            real
            chi();

//----------------------------------------------------------------------------

            /**
             * crossover function
             * @return
             */
            real
            f();


//----------------------------------------------------------------------------

            void
            update_Tp( const real & aT, const real & aP );

//----------------------------------------------------------------------------

            bool
            test( const index_t & aIndex ) const ;

//----------------------------------------------------------------------------

            void
            set( const index_t & aIndex, const real aValue ) ;

//----------------------------------------------------------------------------

            const real &
            get( const index_t & aIndex ) const ;

//----------------------------------------------------------------------------
        } ;

//----------------------------------------------------------------------------

        inline void
        HelmholtzTransport_Methane::update_Tp( const real & aT, const real & aP )
        {
            BELFEM_ASSERT( aT   > 0, "Invalid temprerature" );
            BELFEM_ASSERT( aP > 0, "Invalid pressure" );
            if(               aT != mVals[ BELFEM_METHANE_T ]
                              || aP != mVals[ BELFEM_METHANE_P ] )
            {
                mBits.reset();

                mVals[ BELFEM_METHANE_T ]     = aT ;
                mVals[ BELFEM_METHANE_P ]     = aP ;
                mVals[ BELFEM_METHANE_V ]     = mEoS.v( aT, aP );

                mVals[ BELFEM_METHANE_TAU ]   = mTcrit / aT ;
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
        HelmholtzTransport_Methane::test( const index_t & aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_HELMHOLTZ_NUMVALS ,
                          "Invalid Helmholtz state index: %u", ( unsigned int ) aIndex );

            return mBits.test( aIndex );
        }

//----------------------------------------------------------------------------

        inline void
        HelmholtzTransport_Methane::set( const index_t & aIndex, const real aValue )
        {
            BELFEM_ASSERT( aIndex < BELFEM_HELMHOLTZ_NUMVALS ,
                          "Invalid Helmholtz state index", ( unsigned int ) aIndex );

            // set value
            mVals[ aIndex ] = aValue;

            // update flag
            mBits.set( aIndex );
        }

//----------------------------------------------------------------------------

        inline const real &
        HelmholtzTransport_Methane::get( const index_t & aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_METHANE_NUMVALS ,
                          "Invalid Helmholtz state index: %u", ( unsigned int ) aIndex );

            return mVals[ aIndex ];
        }

//----------------------------------------------------------------------------

    }
}

#endif //BELFEM_CL_GM_HelmholtzTransportMETHANE_HPP
