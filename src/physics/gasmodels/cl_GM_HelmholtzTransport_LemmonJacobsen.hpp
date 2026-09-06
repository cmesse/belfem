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

#ifndef BELFEM_CL_GM_HELMHOLTZTRANSPORT_LEMMONJACOBSEN_HPP
#define BELFEM_CL_GM_HELMHOLTZTRANSPORT_LEMMONJACOBSEN_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Bitset.hpp"
#include "en_Helmholtz.hpp"
#include "cl_GM_HelmholtzTransport.hpp"

#define BELFEM_LJTRANS_T          0
#define BELFEM_LJTRANS_P          1
#define BELFEM_LJTRANS_V          2
#define BELFEM_LJTRANS_TAU        3
#define BELFEM_LJTRANS_DELTA      4
#define BELFEM_LJTRANS_ETA0       5
#define BELFEM_LJTRANS_ETA        6
#define BELFEM_LJTRANS_LAMBDA     7
#define BELFEM_LJTRANS_NUMVALS    8

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * Viscosity and thermal conductivity after Lemmon and Jacobsen,
         * Int. J. Thermophys. 25:21-69 (2004), doi 10.1023/B:IJOT.0000022327.04529.f3
         *
         * One functional form serves nitrogen, argon, oxygen and air; only the
         * coefficient tables differ. This class therefore carries the fluid as
         * a parameter rather than existing once per species, the way
         * EoS_Hydrogen carries its three spin isomers. Nitrogen and oxygen are
         * implemented, being the two fluids BELFEM has a Helmholtz EoS for.
         *
         * The structure is the same as HelmholtzTransport_Methane -- dilute gas
         * plus residual plus a critical enhancement, served from a cache keyed
         * on the last state -- but the equations are not Friend's; do not carry
         * expressions across between the two.
         *
         *   eta    = eta_0( T ) + eta_r( tau, delta )                   Eq. (1)
         *   lambda = lambda_0( T ) + lambda_r( tau, delta )
         *            + lambda_c( tau, delta )                           Eq. (4)
         *
         * with tau = T_crit / T and delta = rho / rho_crit. There is no
         * critical enhancement on the viscosity: the paper judges it negligible
         * for practical states and fits none.
         *
         * The paper works in micro Pa s and milli W / ( m K ); the accessors
         * return SI, Pa s and W / ( m K ).
         *
         * @ingroup grp_physics_gasmodels
         */
        class HelmholtzTransport_LemmonJacobsen : public HelmholtzTransport
        {
//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            //! which fluid the coefficient tables belong to
            const HelmholtzModel mModel ;

            //! critical temperature in K, Table I
            real mTcrit ;

            //! critical pressure in Pa, Table I
            real mPcrit ;

            //! critical mass density in kg/m^3, Table I converted with M
            real mRhocrit ;

            //! molar mass in g/mol, Table I. The dilute gas equation is
            //! written for that unit, so it is kept as printed
            real mM ;

            //! Lennard-Jones energy parameter epsilon / k in K, Table I
            real mEpsilonKb ;

            //! Lennard-Jones size parameter sigma in nm, Table I
            real mSigma ;

            //! collision integral coefficients b_i, Table II.
            //! shared by every fluid of the correlation
            Vector< real > mB ;

            //! residual viscosity, Table III
            Vector< real > mEtaN ;
            Vector< real > mEtaT ;
            Vector< real > mEtaD ;
            Vector< real > mEtaL ;

            //! dilute gas thermal conductivity, Table IV rows 1 to 3
            real mLambdaN1 ;
            real mLambdaN2 ;
            real mLambdaT2 ;
            real mLambdaN3 ;
            real mLambdaT3 ;

            //! residual thermal conductivity, Table IV from row 4
            Vector< real > mLambdaN ;
            Vector< real > mLambdaT ;
            Vector< real > mLambdaD ;
            Vector< real > mLambdaL ;

            //! critical enhancement, fluid specific terms of Table I
            real mXi0 ;      //!< xi_0 in m ( printed in nm )
            real mGammaC ;   //!< Gamma, dimensionless
            real mQd ;       //!< q_D in m ( printed in nm )
            real mTref ;     //!< reference temperature in K, twice T_crit

            //! constant prefactor of the dilute gas viscosity
            real mConstEta0 ;

            //! state cache. Logically const, written by the const evaluators
            mutable real mVals[ BELFEM_LJTRANS_NUMVALS ] = { 0.0 };
            mutable Bitset< BELFEM_LJTRANS_NUMVALS > mBits;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            /**
             * @param aParent  gas holding the Helmholtz EoS of this fluid
             * @param aModel   HelmholtzModel::Nitrogen or HelmholtzModel::Oxygen
             */
            HelmholtzTransport_LemmonJacobsen(
                    Gas & aParent,
                    const HelmholtzModel aModel );

            ~HelmholtzTransport_LemmonJacobsen() = default ;

//----------------------------------------------------------------------------

            /**
             * dynamic viscosity in Pa*s
             */
            real
            mu( const real T, const real p ) const ;

//----------------------------------------------------------------------------

            /**
             * thermal conductivity in W/(m*K)
             */
            real
            lambda( const real T, const real p ) const ;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            //! coefficient tables of the selected fluid
            void
            init_tables() ;

//----------------------------------------------------------------------------

            //! collision integral Omega( T* ), Table II
            real
            omega( const real T ) const ;

//----------------------------------------------------------------------------

            //! dilute gas viscosity in micro Pa s, Eq. ( 2 )
            real
            eta_0( const real T ) const ;

//----------------------------------------------------------------------------

            //! residual viscosity in micro Pa s, Eq. ( 3 )
            real
            eta_r() const ;

//----------------------------------------------------------------------------

            //! dilute gas thermal conductivity in mW/(m K), Eq. ( 5 )
            real
            lambda_0( const real T ) const ;

//----------------------------------------------------------------------------

            //! residual thermal conductivity in mW/(m K), Eq. ( 6 )
            real
            lambda_r() const ;

//----------------------------------------------------------------------------

            /**
             * critical enhancement of the thermal conductivity in W/(m K),
             * Olchowy and Sengers as repeated in Eqs. ( 7 ) to ( 11 )
             */
            real
            lambda_c( const real T, const real p ) const ;

//----------------------------------------------------------------------------

            //! symmetrised compressibility chi tilde, Eq. ( 11 )
            real
            chi( const real T, const real v ) const ;

//----------------------------------------------------------------------------

            void
            update_Tp( const real T, const real p ) const ;

            bool
            test( const index_t aIndex ) const ;

            void
            set( const index_t aIndex, const real aValue ) const ;

            const real &
            get( const index_t aIndex ) const ;

//----------------------------------------------------------------------------
        } ;

//----------------------------------------------------------------------------

        inline bool
        HelmholtzTransport_LemmonJacobsen::test( const index_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_LJTRANS_NUMVALS,
                          "Invalid Lemmon-Jacobsen state index: %u",
                          ( unsigned int ) aIndex );

            return mBits.test( aIndex );
        }

//----------------------------------------------------------------------------

        inline void
        HelmholtzTransport_LemmonJacobsen::set(
                const index_t aIndex, const real aValue ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_LJTRANS_NUMVALS,
                          "Invalid Lemmon-Jacobsen state index: %u",
                          ( unsigned int ) aIndex );

            mVals[ aIndex ] = aValue;
            mBits.set( aIndex );
        }

//----------------------------------------------------------------------------

        inline const real &
        HelmholtzTransport_LemmonJacobsen::get( const index_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_LJTRANS_NUMVALS,
                          "Invalid Lemmon-Jacobsen state index: %u",
                          ( unsigned int ) aIndex );

            return mVals[ aIndex ];
        }

//----------------------------------------------------------------------------

        inline void
        HelmholtzTransport_LemmonJacobsen::update_Tp(
                const real T, const real p ) const
        {
            BELFEM_ASSERT( T > 0, "Invalid temperature" );
            BELFEM_ASSERT( p > 0, "Invalid pressure" );

            if(    T != mVals[ BELFEM_LJTRANS_T ]
                || p != mVals[ BELFEM_LJTRANS_P ] )
            {
                mBits.reset();

                mVals[ BELFEM_LJTRANS_V ]     = mEoS.v( T, p );

                mVals[ BELFEM_LJTRANS_T ]     = T ;
                mVals[ BELFEM_LJTRANS_P ]     = p ;
                mVals[ BELFEM_LJTRANS_TAU ]   = mTcrit / T ;
                mVals[ BELFEM_LJTRANS_DELTA ] =
                        1.0 / ( mVals[ BELFEM_LJTRANS_V ] * mRhocrit ) ;

                mBits.set( BELFEM_LJTRANS_T );
                mBits.set( BELFEM_LJTRANS_P );
                mBits.set( BELFEM_LJTRANS_V );
                mBits.set( BELFEM_LJTRANS_TAU );
                mBits.set( BELFEM_LJTRANS_DELTA );
            }
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HELMHOLTZTRANSPORT_LEMMONJACOBSEN_HPP
