//
// Created by Christian Messe on 17.08.20.
//

#ifndef BELFEM_CL_GM_HELMHOLTZ_HPP
#define BELFEM_CL_GM_HELMHOLTZ_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Bitset.hpp"
#include "cl_Vector.hpp"
#include "cl_GM_EoS.hpp"
#include "cl_GM_EoS_Cubic.hpp"

#define BELFEM_HELMHOLTZ_T                 0
#define BELFEM_HELMHOLTZ_P                 1
#define BELFEM_HELMHOLTZ_V                 2
#define BELFEM_HELMHOLTZ_TAU               3
#define BELFEM_HELMHOLTZ_DELTA             4
#define BELFEM_HELMHOLTZ_PHI0              5
#define BELFEM_HELMHOLTZ_PHI0_T            6
#define BELFEM_HELMHOLTZ_PHI0_TT           7
#define BELFEM_HELMHOLTZ_PHIR              8
#define BELFEM_HELMHOLTZ_PHIR_T            9
#define BELFEM_HELMHOLTZ_PHIR_D           10
#define BELFEM_HELMHOLTZ_PHIR_TT          11
#define BELFEM_HELMHOLTZ_PHIR_DD          12
#define BELFEM_HELMHOLTZ_PHIR_DT          13
#define BELFEM_HELMHOLTZ_DPDV             14
#define BELFEM_HELMHOLTZ_DPDT             15
#define BELFEM_HELMHOLTZ_DVDT             16
#define BELFEM_HELMHOLTZ_NUMVALS          17

namespace belfem
{
    // forward declaration for parent
    class Gas;

    namespace gasmodels
    {
        // forward declaration for statevals
        class Statevals;

//----------------------------------------------------------------------------
        /**
         * a model for the helmholtz energy, specifically user for
         * cryogenic fluids
         */
        class Helmholtz : public EoS
        {
//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            //Gas        & mParent ;
            //Statevals  & mStatevals ;

            //const real & mR ;
            //const real & mM ;

            const      string mLabel ;

            // help model to find initial solution
            EoS_Cubic * mCubicEoS = nullptr ;

            Vector< real > mNvap ;
            Vector< real > mKvap ;
            Vector< real > mTvap ;  //! initial solution for inversion of pvap

            // initial guess for specific volume in liquid state
            Vector< real > mVliq ;

            // enthalpy offset
            real mH0 = 0.0 ;

            // entropy offset
            real mS0 = 0.0 ;

            // inner energy offset
            real mU0 = 0.0 ;

            friend void deriv_test( Helmholtz &, Vector< real > & aR2 );

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            /**
             * internal state variables
             */
            real mHelmholtzVals[ BELFEM_HELMHOLTZ_NUMVALS ] = { 0.0 };
            Bitset< BELFEM_HELMHOLTZ_NUMVALS > mHelmholtzBits;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            // shortcuts
            const real & mTau   = mHelmholtzVals[ BELFEM_HELMHOLTZ_TAU ];
            const real & mDelta = mHelmholtzVals[ BELFEM_HELMHOLTZ_DELTA ];

            //! critical temperature in K
            real mTcrit   = BELFEM_QUIET_NAN ;

            //! critical pressure in Pa
            real mPcrit   = BELFEM_QUIET_NAN ;

            //! critical density in kg/m^3
            real mRhocrit = BELFEM_QUIET_NAN ;

            //! critical specific volume in m^3/kg
            real mVcrit   = BELFEM_QUIET_NAN ;

            //! triple point in K
            real mTtriple = BELFEM_QUIET_NAN ;

            //! maximum temperature of gas model in K
            real mTmax = BELFEM_QUIET_NAN ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Helmholtz( Gas & aParent, const string & aLabel );

            virtual ~Helmholtz() = default ;

//----------------------------------------------------------------------------

            /**
             * minumum temperature for gas model ( triple point )
             */
             const real &
             T_min() const ;

//----------------------------------------------------------------------------

            /**
             * maximum temperature for gas model
             */
            const real &
            T_max() const ;

//----------------------------------------------------------------------------

            /**
             * critical temperature
             */
            const real &
            T_crit() const ;

//----------------------------------------------------------------------------

            /**
             * critical volume
             */
            const real &
            v_crit() const ;

//----------------------------------------------------------------------------

            /**
             * @param aT  vapor temperature in K
             * @return    vapor pressure in Pa
             */
            virtual real
            p_vap( const real aT );

//----------------------------------------------------------------------------

            /**
             * @param aP  vapor pressure in Pa
             * @return    vapor temperature in K
             */
            virtual real
            T_vap( const real aP );

//----------------------------------------------------------------------------

            /**
             * vaporization enthalpy in J/kg
             * @param aT
             * @param aP
             * @return
             */
            real
            hvap( const real aT, const real aP ) ;

//----------------------------------------------------------------------------

           /**
           * the main energy function,
           * where phi = f / ( R * T )
           * @param aT    : T
           * @param aV    : v
           * @return
           */
           real
           phi( const real aT, const real aV );

//----------------------------------------------------------------------------

            /**
             * pressure in Pa
             * @param aT
             * @param aV
             * @return
             */
           real
           p( const real aT, const real aV );

//----------------------------------------------------------------------------

           real
           v( const real aT, const real aP );

//----------------------------------------------------------------------------

           real
           T( const real aP, const real aV );

//----------------------------------------------------------------------------

           real
           dpdv( const real aT, const real aV );

//----------------------------------------------------------------------------

           real
           dpdT( const real aT, const real aV );

//----------------------------------------------------------------------------

           real
           dvdT( const real aT, const real aV );

//----------------------------------------------------------------------------
// Caloric Functions
//----------------------------------------------------------------------------

           real
           u( const real aT, const real aP );

           real
           h( const real aT, const real aP );

           real
           s( const real aT, const real aP );

           real
           cv( const real aT, const real aP );

           real
           cp( const real aT, const real aP );

           real
           w( const real aT, const real aP );

           real
           dsdT( const real aT, const real aP );

           real
           dsdp( const real aT, const real aP );

           real
           alpha( const real aT, const real aP );

           real
           beta( const real aT, const real aP );

           real
           kappa( const real aT, const real aP );

//----------------------------------------------------------------------------

            /**
             * this function remixes the cubic help gas
             */
            void
            remix();

//----------------------------------------------------------------------------

            /**
             * return the critical point data
             * @param aT
             * @param aP
             * @param aV
             */
            void
            eval_critical_point( real & aT, real & aP, real & aV );


//----------------------------------------------------------------------------

            /**
             * reset all state variales
             */
            void
            update_Tv( const real aT, const real aV ) ;

//----------------------------------------------------------------------------

            void
            update_Tp( const real aT, const real aP ) ;

//----------------------------------------------------------------------------

            /**
             * the main energy function, ideal gas contribution
             * @param aTau     : T / T_crit
             * @param aDelta   : v_crit / v
             * @return
             */
            const real &
            phi0();

            const real &
            phi0_d();

            const real &
            phi0_dd();

            const real &
            phi0_t();

            const real &
            phi0_tt();

//----------------------------------------------------------------------------

            /**
             * the main energy function, real gas contribution
             */
            const real &
            phir();

            const real &
            phir_d();

            const real &
            phir_dd();

            const real &
            phir_t();

            const real &
            phir_tt();

            const real &
            phir_dt();

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            /**
             * test if a value is up to date
             */
           bool
           test( const index_t aIndex ) const ;

//----------------------------------------------------------------------------

            /**
             * get a value from the memory
             */
            const real &
            get( const index_t aIndex ) const ;

//----------------------------------------------------------------------------

            /**
             * write a value into the memory
             */
            void
            set( const index_t aIndex, const real aValue ) ;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            // tidy up parent class
            void
            delete_cubic_eos();

//----------------------------------------------------------------------------

            // make data object consisient
            void
            set_critical_point_in_data_object();

//----------------------------------------------------------------------------

            // make sure that parent has the right type
            void
            check_parent() ;

//----------------------------------------------------------------------------

            // initialize the offsets for enthalpy and entropy
            // with respect to CEA table
            void
            set_reference_point();

//----------------------------------------------------------------------------
            // initialize the offsets for enthalpy and entropy
            void
            set_reference_point( const real aTref, const real aPref );

//----------------------------------------------------------------------------
            // initialize lookup table for initial guess of Tvap
            void
            init_Tvap_poly();

//----------------------------------------------------------------------------

            // help function for vapor pressure
            virtual real
            pi_vap( const real aT ) ;

//----------------------------------------------------------------------------

            // help function for vapor pressure
            virtual real
            psi_vap( const real aT ) ;

//----------------------------------------------------------------------------

            // derivative for vapor pressure
            virtual real
            dpvap_dT( const real aT, const real aPvap, const real aPiVap );

//----------------------------------------------------------------------------

            // check of this state is liquid or gaseous
            bool
            is_liquid( const real aT, const real aP ) ;

//----------------------------------------------------------------------------

            virtual real
            compute_phi0() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phi0_t() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phi0_tt() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_d() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_dd() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_t() ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_tt() ;


//----------------------------------------------------------------------------

            virtual real
            compute_phir_dt() ;

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

        inline const real &
        Helmholtz::T_min() const
        {
            return mTtriple ;
        }

//----------------------------------------------------------------------------

        inline const real &
        Helmholtz::T_max() const
        {
            return mTmax ;
        }

//----------------------------------------------------------------------------

        inline const real &
        Helmholtz::T_crit() const
        {
            return mTcrit ;
        }

//----------------------------------------------------------------------------

        inline const real &
        Helmholtz::v_crit() const
        {
            return mVcrit ;
        }

//----------------------------------------------------------------------------

        inline bool
        Helmholtz::test( const index_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_HELMHOLTZ_NUMVALS ,
                          "Invalid Helmholtz state index: %u", ( unsigned int ) aIndex );

            return mHelmholtzBits.test( aIndex );
        }

//----------------------------------------------------------------------------

        inline const real &
        Helmholtz::get( const index_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_HELMHOLTZ_NUMVALS ,
                          "Invalid Helmholtz state index: %u", ( unsigned int ) aIndex );

            return mHelmholtzVals[ aIndex ];
        }

//----------------------------------------------------------------------------

        inline void
        Helmholtz::set( const index_t aIndex, const real aValue )
        {
            BELFEM_ASSERT( aIndex < BELFEM_HELMHOLTZ_NUMVALS ,
                          "Invalid Helmholtz state index", ( unsigned int ) aIndex );

            // set value
            mHelmholtzVals[ aIndex ] = aValue;

            // update flag
            mHelmholtzBits.set( aIndex );
        }

//----------------------------------------------------------------------------

        inline void
        Helmholtz::update_Tv( const real aT, const real aV )
        {
            BELFEM_ASSERT( aT   > 0, "Invalid Temprerature" );
            BELFEM_ASSERT( aV > 0, "Invalid specific volume" );
            if(               aT != mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]
                           || aV != mHelmholtzVals[ BELFEM_HELMHOLTZ_V ] )
            {
                mHelmholtzBits.reset();

                mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]     = aT ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_V ]     = aV ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_TAU ]   = mTcrit / aT ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_DELTA ] = mVcrit / aV ;

                mHelmholtzBits.set( BELFEM_HELMHOLTZ_T     );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_V     );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_TAU   );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_DELTA );
            }
        }

//----------------------------------------------------------------------------

        inline void
        Helmholtz::update_Tp( const real aT, const real aP )
        {
            BELFEM_ASSERT( aT   > 0, "Invalid temprerature" );
            BELFEM_ASSERT( aP > 0, "Invalid pressure" );
            if(               aT != mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]
                           || aP != mHelmholtzVals[ BELFEM_HELMHOLTZ_P ] )
            {
                mHelmholtzBits.reset();

                mHelmholtzVals[ BELFEM_HELMHOLTZ_V ]     = this->v( aT, aP );
                mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]     = aT ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_P ]     = aP ;

                mHelmholtzVals[ BELFEM_HELMHOLTZ_TAU ]   = mTcrit / aT ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_DELTA ] = mVcrit / mHelmholtzVals[ BELFEM_HELMHOLTZ_V ] ;


                mHelmholtzBits.set( BELFEM_HELMHOLTZ_T );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_P );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_V );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_TAU );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_DELTA );
            }
        }

//----------------------------------------------------------------------------

        // help function for vapor pressure
        inline real
        Helmholtz::pi_vap( const real aT )
        {
            real tTheta = std::abs( 1.0 - aT / mTcrit ) ;

            return (   mNvap( 0 ) * std::pow( tTheta, mKvap( 0 ) )
                     + mNvap( 1 ) * std::pow( tTheta, mKvap( 1 ) )
                     + mNvap( 2 ) * std::pow( tTheta, mKvap( 2 ) )
                     + mNvap( 3 ) * std::pow( tTheta, mKvap( 3 ) ) )
                     * mTcrit / aT ;
        }

//----------------------------------------------------------------------------

        // help function for vapor pressure
        inline real
        Helmholtz::psi_vap( const real aT )
        {
            real tTheta = std::abs( 1.0 - aT / mTcrit ) ;

            return (   mKvap( 0 ) * mNvap( 0 ) * std::pow( tTheta, mKvap( 0 ) - 1.0 )
                     + mKvap( 1 ) * mNvap( 1 ) * std::pow( tTheta, mKvap( 1 ) - 1.0 )
                     + mKvap( 2 ) * mNvap( 2 ) * std::pow( tTheta, mKvap( 2 ) - 1.0 )
                     + mKvap( 3 ) * mNvap( 3 ) * std::pow( tTheta, mKvap( 3 ) - 1.0 ) ) ;
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HELMHOLTZ_HPP
