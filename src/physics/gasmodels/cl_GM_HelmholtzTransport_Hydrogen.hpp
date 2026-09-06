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

#ifndef BELFEM_CL_GM_HELMHOLTZTRANSPORT_HYDROGEN_HPP
#define BELFEM_CL_GM_HELMHOLTZTRANSPORT_HYDROGEN_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "en_Helmholtz.hpp"
#include "cl_GM_HelmholtzTransport.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * Viscosity of normal hydrogen after Muzny, Huber and Kazakov,
         * J. Chem. Eng. Data 58:969-979 (2013), doi 10.1021/je301273j,
         * **with the 2022 erratum applied**, J. Chem. Eng. Data 67:2855,
         * doi 10.1021/acs.jced.2c00523.
         *
         * The erratum matters. Without it the correlation is wrong, and it
         * changes three separate things:
         *
         *   1. Eq. ( 6 ) is missing Avogadro's number. The second viscosity
         *      virial coefficient carries N_A, and the paper as printed does
         *      not say so.
         *   2. In Eq. ( 7 ) the exponent on T* is **-i**, not +i.
         *   3. The density scale of Eq. ( 9 ) is rho_sc = 90.909090909 kg/m^3.
         *      The body text of the paper says 90.5.
         *
         * The erratum also supplies three test values, which this
         * implementation reproduces to better than 0.001 %:
         *
         *   T = 40 K, rho =   0 kg/m^3 : eta =  1.9772 micro Pa s
         *   T = 40 K, rho =  50 kg/m^3 : eta =  5.9905 micro Pa s
         *   T = 40 K, rho = 100 kg/m^3 : eta = 49.034  micro Pa s
         *
         * The correlation is for **normal** hydrogen. BELFEM's Helmholtz EoS
         * distinguishes the para, normal and ortho isomers, but no transport
         * correlation for the spin isomers exists; this one is installed for
         * all three, which is what REFPROP does as well. The density that
         * enters comes from whichever equation of state the gas actually
         * carries, so the isomers do differ in the result, just not in the
         * correlation.
         *
         * The **thermal conductivity** is a second, unrelated correlation:
         * Assael, Assael, Huber, Perkins and Takata, J. Phys. Chem. Ref. Data
         * 40:033101 (2011), doi 10.1063/1.3606499. Unlike the viscosity paper
         * that one does distinguish the isomers and carries separate
         * coefficient tables for normal and parahydrogen, so lambda() is
         * isomer aware where mu() is not. See the note on isomer coverage in
         * the module README.
         *
         * @ingroup grp_physics_gasmodels
         */
        class HelmholtzTransport_Hydrogen : public HelmholtzTransport
        {
//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            //! molar mass in g/mol, as Eq. ( 3 ) is written
            const real mM = 2.01588 ;

            //! Lennard-Jones length scale in nm, after Behnejad and Miralinaghi
            const real mSigma = 0.297 ;

            //! Lennard-Jones energy scale epsilon / k_B in K
            const real mEpsilonKb = 30.41 ;

            //! critical temperature of the density scaling, in K
            const real mTcrit = 33.145 ;

            //! density scale of Eq. ( 9 ) in kg/m^3, **erratum value**
            const real mRhoScale = 90.909090909 ;

            //! effective cross section, Table 2
            Vector< real > mA ;

            //! second viscosity virial, Table 3
            Vector< real > mB ;

            //! symbolic regression term, Table 4
            Vector< real > mC ;

            //! prefactor of Eq. ( 3 )
            real mConstEta0 ;

            //! sigma^3 * N_A / M, the factor of Eq. ( 6 ) after the erratum,
            //! giving the second viscosity virial in m^3/kg
            real mConstBeta ;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  thermal conductivity, Assael et al. 2011
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            //! which spin isomer the coefficient tables belong to
            const HelmholtzModel mModel ;

            //! critical point of the isomer, taken from the equation of state
            //! so that the reduced variables match the density it returns
            real mTcritEoS ;
            real mPcritEoS ;
            real mRhocritEoS ;

            //! dilute gas, Assael Eq. ( 2 ), numerator and denominator
            Vector< real > mLambdaA1 ;
            Vector< real > mLambdaA2 ;

            //! excess, Assael Eq. ( 3 )
            Vector< real > mLambdaB1 ;
            Vector< real > mLambdaB2 ;

            //! critical enhancement, Assael Sec. 3.3.1
            const real mQd    = 1.0 / 4.0e-10 ;  //!< cut-off wavenumber, 1/m
            const real mXi0   = 1.5e-10 ;        //!< amplitude, m
            const real mGammaC = 0.052 ;         //!< amplitude, dimensionless
            real       mTrefC ;                  //!< 3/2 * T_crit

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            /**
             * @param aParent  gas holding the Helmholtz EoS of this fluid
             * @param aModel   which spin isomer the parent carries
             */
            HelmholtzTransport_Hydrogen(
                    Gas & aParent,
                    const HelmholtzModel aModel );

            ~HelmholtzTransport_Hydrogen() = default ;

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

            //! zero density viscosity in micro Pa s, Eqs. ( 3 ) and ( 4 )
            real
            eta_0( const real T ) const ;

//----------------------------------------------------------------------------

            //! initial density coefficient in micro Pa s per kg/m^3,
            //! Eqs. ( 5 ) to ( 7 )
            real
            eta_1( const real T ) const ;

//----------------------------------------------------------------------------

            //! dilute gas thermal conductivity in W/(m K), Assael Eq. ( 2 )
            real
            lambda_0( const real T ) const ;

//----------------------------------------------------------------------------

            //! excess thermal conductivity in W/(m K), Assael Eq. ( 3 )
            real
            lambda_e( const real T, const real rho ) const ;

//----------------------------------------------------------------------------

            //! critical enhancement in W/(m K), Assael Eqs. ( 4 ) to ( 7 )
            real
            lambda_c( const real T, const real p, const real rho ) const ;

//----------------------------------------------------------------------------

            //! symmetrised compressibility of Assael Eq. ( 7 )
            real
            chi( const real T, const real v ) const ;

//----------------------------------------------------------------------------
        } ;

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HELMHOLTZTRANSPORT_HYDROGEN_HPP
