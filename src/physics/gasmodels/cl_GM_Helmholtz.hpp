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
         *
         * @ingroup grp_physics_gasmodels
         * @see @ref physics_gasmodels_gasmodels_usage_guide
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
             *
             * Logically const cache: the values are a memo of the state the
             * caller last asked for, not part of the identity of the fluid.
             * The const evaluators therefore write here through mutable.
             * As a consequence, a const Helmholtz is NOT reentrant and not
             * safe to share between threads -- see the note in cl_Gas.hpp.
             */
            mutable real mHelmholtzVals[ BELFEM_HELMHOLTZ_NUMVALS ] = { 0.0 };
            mutable Bitset< BELFEM_HELMHOLTZ_NUMVALS > mHelmholtzBits;

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

            // upper pressure bound of the published validity range, set by
            // the child EoS; enforced with BELFEM_ERROR in Helmholtz::v()
            real mPmax = BELFEM_REAL_MAX ;

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
             * @param T  vapor temperature in K
             * @return    vapor pressure in Pa
             */
            virtual real
            p_vap( const real T ) const;

//----------------------------------------------------------------------------

            /**
             * @param p  vapor pressure in Pa
             * @return    vapor temperature in K
             */
            virtual real
            T_vap( const real p ) const;

//----------------------------------------------------------------------------

            /**
             * vaporization enthalpy in J/kg
             * @param T
             * @param p
             * @return
             */
            real
            hvap( const real T, const real p ) const ;

//----------------------------------------------------------------------------

           /**
           * the main energy function,
           * where phi = f / ( R * T )
           * @param T    temperature in K
           * @param v    mass specific volume in m^3/kg
           * @return
           */
           real
           phi( const real T, const real v ) const;

//----------------------------------------------------------------------------

            /**
             * pressure in Pa
             * @param T
             * @param v
             * @return
             */
           real
           p( const real T, const real v ) const;

//----------------------------------------------------------------------------

           real
           v( const real T, const real p ) const;

//----------------------------------------------------------------------------

           real
           T( const real p, const real v ) const;

//----------------------------------------------------------------------------

           real
           dpdv( const real T, const real v ) const;

//----------------------------------------------------------------------------

           real
           dpdT( const real T, const real v ) const;

//----------------------------------------------------------------------------

           real
           dvdT( const real T, const real v ) const;

//----------------------------------------------------------------------------
// Caloric Functions
//----------------------------------------------------------------------------

           real
           u( const real T, const real p ) const;

           real
           h( const real T, const real p ) const;

           real
           s( const real T, const real p ) const;

           real
           cv( const real T, const real p ) const;

           real
           cp( const real T, const real p ) const;

           real
           w( const real T, const real p ) const;

           real
           dsdT( const real T, const real p ) const;

           real
           dsdp( const real T, const real p ) const;

           real
           alpha( const real T, const real p ) const;

           real
           beta( const real T, const real p ) const;

           real
           kappa( const real T, const real p ) const;

//----------------------------------------------------------------------------

            /**
             * this function remixes the cubic help gas
             */
            void
            remix();

//----------------------------------------------------------------------------

            /**
             * return the critical point data
             * @param T
             * @param p
             * @param v
             */
            void
            eval_critical_point( real & T, real & p, real & v ) const;


//----------------------------------------------------------------------------

            /**
             * reset all state variales
             */
            void
            update_Tv( const real T, const real v ) const ;

//----------------------------------------------------------------------------

            void
            update_Tp( const real T, const real p ) const ;

//----------------------------------------------------------------------------

            /**
             * @name Ideal gas contribution to the reduced Helmholtz energy
             *
             * phi0 and its first and second derivatives with respect to tau.
             * They take no arguments: the state is whatever update_Tp or
             * update_Tv last set, and tau = T_crit / T and delta = v_crit / v
             * are read from the cache. Each is evaluated once per state and
             * then served from that cache.
             * @{
             */
            const real &
            phi0() const;

            const real &
            phi0_t() const;

            const real &
            phi0_tt() const;

            /** @} */

//----------------------------------------------------------------------------

            /**
             * the main energy function, real gas contribution
             */
            const real &
            phir() const;

            const real &
            phir_d() const;

            const real &
            phir_dd() const;

            const real &
            phir_t() const;

            const real &
            phir_tt() const;

            const real &
            phir_dt() const;

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
            set( const index_t aIndex, const real aValue ) const ;

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
            pi_vap( const real T ) const ;

//----------------------------------------------------------------------------

            // help function for vapor pressure
            virtual real
            psi_vap( const real T ) const ;

//----------------------------------------------------------------------------

            // derivative for vapor pressure
            virtual real
            dpvap_dT( const real T, const real aPvap, const real aPiVap ) const;

//----------------------------------------------------------------------------

            // check of this state is liquid or gaseous
            bool
            is_liquid( const real T, const real p ) const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phi0() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phi0_t() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phi0_tt() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_d() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_dd() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_t() const ;

//----------------------------------------------------------------------------

            virtual real
            compute_phir_tt() const ;


//----------------------------------------------------------------------------

            virtual real
            compute_phir_dt() const ;

//----------------------------------------------------------------------------

            /**
             * integer power by binary exponentiation, for the coefficient
             * tables of the residual functions. Measured on the 40 term methane
             * table, building tau^t with hpow() instead of std::pow costs
             * 109 ns against 605 ns. Note that std::pow is still used elsewhere
             * in this class, in particular for the vapor pressure ancillary.
             */
            static inline real
            ipow( real x, uint n )
            {
                real r = 1.0 ;

                while( n )
                {
                    if( n & 1 )
                    {
                        r *= x ;
                    }

                    x *= x ;
                    n >>= 1 ;
                }

                return r ;
            }

//----------------------------------------------------------------------------

            /**
             * x^( n / 2 ) for a signed half integer exponent given as n = 2j.
             * The square root is passed in so that the caller hoists it out of
             * the loop over the coefficient table.
             */
            static inline real
            hpow( const real x, const real sqrt_x, const int n )
            {
                const uint tN = n < 0 ? -n : n ;

                real tValue = ipow( x, tN >> 1 ) ;

                if( tN & 1 )
                {
                    tValue *= sqrt_x ;
                }

                return n < 0 ? 1.0 / tValue : tValue ;
            }

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
        Helmholtz::set( const index_t aIndex, const real aValue ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_HELMHOLTZ_NUMVALS ,
                          "Invalid Helmholtz state index %u", ( unsigned int ) aIndex );

            // set value
            mHelmholtzVals[ aIndex ] = aValue;

            // update flag
            mHelmholtzBits.set( aIndex );
        }

//----------------------------------------------------------------------------

        inline void
        Helmholtz::update_Tv( const real T, const real v ) const
        {
            BELFEM_ASSERT( T   > 0, "Invalid Temprerature" );
            BELFEM_ASSERT( v > 0, "Invalid specific volume" );
            if(               T != mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]
                           || v != mHelmholtzVals[ BELFEM_HELMHOLTZ_V ] )
            {
                mHelmholtzBits.reset();

                mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]     = T ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_V ]     = v ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_TAU ]   = mTcrit / T ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_DELTA ] = mVcrit / v ;

                mHelmholtzBits.set( BELFEM_HELMHOLTZ_T     );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_V     );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_TAU   );
                mHelmholtzBits.set( BELFEM_HELMHOLTZ_DELTA );
            }
        }

//----------------------------------------------------------------------------

        inline void
        Helmholtz::update_Tp( const real T, const real p ) const
        {
            BELFEM_ASSERT( T   > 0, "Invalid temprerature" );
            BELFEM_ASSERT( p > 0, "Invalid pressure" );
            if(               T != mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]
                           || p != mHelmholtzVals[ BELFEM_HELMHOLTZ_P ] )
            {
                // solve for v first: its Newton iterations set phir bits at
                // intermediate states, so the reset must come afterwards
                const real tV = this->v( T, p );

                mHelmholtzBits.reset();

                mHelmholtzVals[ BELFEM_HELMHOLTZ_V ]     = tV;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_T ]     = T ;
                mHelmholtzVals[ BELFEM_HELMHOLTZ_P ]     = p ;

                mHelmholtzVals[ BELFEM_HELMHOLTZ_TAU ]   = mTcrit / T ;
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
        Helmholtz::pi_vap( const real T ) const
        {
            real tTheta = std::abs( 1.0 - T / mTcrit ) ;

            return (   mNvap( 0 ) * std::pow( tTheta, mKvap( 0 ) )
                     + mNvap( 1 ) * std::pow( tTheta, mKvap( 1 ) )
                     + mNvap( 2 ) * std::pow( tTheta, mKvap( 2 ) )
                     + mNvap( 3 ) * std::pow( tTheta, mKvap( 3 ) ) )
                     * mTcrit / T ;
        }

//----------------------------------------------------------------------------

        // help function for vapor pressure
        inline real
        Helmholtz::psi_vap( const real T ) const
        {
            real tTheta = std::abs( 1.0 - T / mTcrit ) ;

            return (   mKvap( 0 ) * mNvap( 0 ) * std::pow( tTheta, mKvap( 0 ) - 1.0 )
                     + mKvap( 1 ) * mNvap( 1 ) * std::pow( tTheta, mKvap( 1 ) - 1.0 )
                     + mKvap( 2 ) * mNvap( 2 ) * std::pow( tTheta, mKvap( 2 ) - 1.0 )
                     + mKvap( 3 ) * mNvap( 3 ) * std::pow( tTheta, mKvap( 3 ) - 1.0 ) ) ;
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HELMHOLTZ_HPP
