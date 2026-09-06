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

#include "cl_GM_EoS_Nitrogen.hpp"
#include "fn_dot.hpp"

namespace belfem
{
    namespace gasmodels
    {
        EoS_Nitrogen::EoS_Nitrogen( Gas & aParent ) :
                Helmholtz( aParent, "N2" ),
                mA( { 2.5, -12.76952708, -0.00784163, -1.934819e-4, -1.247742e-5,
                              6.678326e-8, 1.012941, 26.65788 }),
                mN( { 0.924803575275, -0.492448489428, 0.661883336938, -1.92902649201,
                    -0.0622469309629, 0.349943957581, 0.564857472498, -1.61720005987,
                    -0.481395031883, 0.421150636384, -0.0161962230825, 0.172100994165,
                    0.00735448924933, 0.0168077305479, -0.00107626664179, -0.0137318088513,
                    0.000635466899859, 0.00304432279419, -0.0435762336045, -0.0723174889316,
                    0.0389644315272, -0.021220136391, 0.00408822981509, -0.0000551990017984,
                    -0.0462016716479, -0.00300311716011, 0.0368825891208, -0.0025585684622,
                    0.00896915264558, -0.0044151337035, 0.00133722924858, 0.000264832491957,
                    19.6688194015, -20.9115600730, 0.0167788306989, 2627.67566274}),
                mI( { 1, 1, 2, 2, 3, 3, 1, 1, 1, 3, 3, 4, 6, 6,
                            7, 7, 8, 8, 1, 2, 3, 4, 5, 8, 4, 5, 5, 8, 3,
                          5, 6, 9, 1, 1, 3, 2 }),
                // Span et al. Table 17, column j_k
                mJ( { 0.25, 0.875, 0.5, 0.875, 0.375, 0.75, 0.5,
                             0.75, 2., 1.25, 3.5, 1., 0.5, 3., 0., 2.75,
                             0.75, 2.5, 4., 6., 6., 3., 3., 6., 16., 11.,
                             15., 12., 12., 7., 4., 16., 0., 1., 2., 3. }),
                mL( { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 2, 2, 2, 2 }),
            mBeta( { 325., 325., 300., 275. } ),
            mGamma( { 1.16, 1.16, 1.13, 1.25 } ),
            mPhi( { 20., 20., 15., 25. } )
        {
            mDeltaPowI.set_size( 37, BELFEM_QUIET_NAN );
            mDeltaPowL.set_size( 36, BELFEM_QUIET_NAN );
            mTauPowJ.set_size( 37, BELFEM_QUIET_NAN );
            mE.set_size( 8, BELFEM_QUIET_NAN );
            mF.set_size( 38, BELFEM_QUIET_NAN );
            mSwap.set_size( 8, BELFEM_QUIET_NAN );

            this->init_tables() ;
            this->set_reference_point() ;
        }

//----------------------------------------------------------------------------

        EoS_Nitrogen::~EoS_Nitrogen()
        {
            this->delete_cubic_eos();
        }

//----------------------------------------------------------------------------

        void
        EoS_Nitrogen::init_tables()
        {
            // triple point
            mTtriple = 63.151 ;

            /* critical data, Span et al. Table 1. The paper gives the
             * critical density in molar form, 11.1839 mol/dm^3; converted
             * with its own molar mass of 28.01348 g/mol that is
             *
             *   11.1839e3 mol/m^3 * 28.01348e-3 kg/mol = 313.300 kg/m^3
             *
             * BELFEM's own N2 molar mass ( thermo.inp, 28.0134 g/mol ) would
             * give 313.299 instead. The 3e-6 difference is far below the
             * uncertainty of the equation of state, and the paper's own value
             * is used so that rho_c and the coefficient tables stay one
             * consistent set. */
            mTcrit   = 126.192 ;
            mPcrit   = 3.3958e6 ;
            mRhocrit = 313.300 ;
            mVcrit   = 1.0 / mRhocrit ;

            // synch data with data object
            this->set_critical_point_in_data_object();

            // validity range of this model
            mTmax = 1000.0 ;
            mPmax = 2200.0e6 ;

            /* vapor pressure coefficients, Span et al. ancillary
             *
             *   ln( p / pc ) = ( Tc / T ) * sum_i N_i * theta^k_i ,
             *   theta = 1 - T / Tc
             *
             * which is the form Helmholtz::pi_vap() already evaluates, so
             * this fluid needs no override of pi_vap / psi_vap.
             *
             * Checked against two points the ancillary does not fit by
             * construction: the normal boiling point ( 101325 Pa at
             * 77.355 K ) comes out to -0.004 % and the triple point
             * ( 12520 Pa at 63.151 K ) to +0.014 %. Both are sensitive to
             * the sign of the third coefficient: with it positive they move
             * to +26 % and +72 %. */
            mNvap = { -6.12445284, 1.26327220, -0.765910082, -1.77570564 } ;
            mKvap = { 1.0, 1.5, 2.5, 5.0 } ;

            /* initial guess for the liquid specific volume, only the starting
             * point of the Newton in Helmholtz::v(). Least squares quadratic
             * on Span et al.'s saturated liquid density ancillary,
             *
             *   rho' / rho_c = exp( sum_i N_i * theta^k_i ) ,  theta = 1 - T / Tc
             *   N = {  1.48654237, -0.280476066, 0.0894143085, -0.119879866 }
             *   k = {  0.3294, 2/3, 8/3, 35/6 }
             *
             * fitted over [ T_triple, 0.95 * T_crit ], where it is within
             * 3.8 % ( rms 1.1 % ). The ancillary's leading theta^0.3294 has
             * infinite slope at the critical point, so no polynomial follows
             * it there; above the fit range the guess degrades, reaching
             * +61 % at T_crit itself. That is the same limitation the oxygen
             * and methane fits carry, and rather milder - the oxygen fit is
             * +93 % at its own critical point. Only the starting point of a
             * Newton that converges to |dp|/p < 1e-8, so it costs iterations,
             * never accuracy. */
            mVliq = { 1.948476e-07, -2.417751e-05, 1.931704e-03 } ;

            // vapor temperature, initial solution
            this->init_Tvap_poly() ;
        }

        real
        EoS_Nitrogen::compute_phi0() const
        {
            this->update_e();

            return std::log( mDelta ) + dot( mA, mE );
        }

        real
        EoS_Nitrogen::compute_phi0_t() const
        {
            this->update_e();

            return
            // ∂/∂τ( a1 * ln τ ) = a1 / τ
            mA( 0 ) * mE( 3 )

            // ∂/∂τ( a3 * τ ) = a3
            + mA( 2 )

            // ∂/∂τ( a4 / τ ) = -a4/τ²
            - mA( 3 ) * mE( 4 )

            // ∂/∂τ( a5 / τ² ) = -2 * a5/τ³
            - 2. * mA( 4 ) * mE( 5 )

            // ∂/∂τ( a6 / τ³ ) = -3 * a6/τ⁴
            - 3. * mA( 5 ) * mE( 4 ) * mE( 4 )

            // ∂/∂τ( a7 * ln( 1 - exp( - a8 * tau ) = ( a7*a8* exp( - a8 * tau )/(1- exp( - a8 * tau ))
            + mA( 6 ) * mA( 7 ) * mChi / ( 1. - mChi );
        }

        real
        EoS_Nitrogen::compute_phi0_tt() const
        {
            this->update_e();

            real x = mChi;

            return
            // ∂²/∂τ²( a1 * ln τ ) = -a1 / τ²
            - mA( 0 ) * mE( 4 )

            // ∂²/∂τ²( a4 / τ ) = 2 * a4 /τ³
            + 2. * mA( 3 ) * mE( 5 )

            // ∂²/∂τ²( a5 / τ² ) = 6 * a5/τ⁴
            + 6. * mA( 4 ) * mE( 4 ) * mE( 4 )

            // ∂²/∂τ²( a6 / τ³ ) = 12 * a6/τ⁵
            + 12. * mA( 5 ) * mE( 4 ) * mE( 5 )


            // ∂²/∂τ²( a7 * ln( 1 - exp( - a8 * tau ) ) = -a7*a8²*exp(...)/(1-exp(...))²
            - mA( 6 ) * mA( 7 ) * mA( 7 ) * ( 1. + x / ( 1. -x ) ) * x / ( 1. - x ) ;
        }

        real
        EoS_Nitrogen::compute_phir() const
        {
            this->update_f();

            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k );
            }
            for ( uint k=6; k<36; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k );
            }
            return alpha_r;
        }

        real
        EoS_Nitrogen::compute_phir_d() const
        {
            this->update_f();


            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mI( k ) ;
            }

            // ∂/∂δ( δ^i * exp( -δ^l ) ) = exp( -δ^l ) * ( i * δ^i - l * δ^(i+l) ) / δ
            for ( uint k=6; k<32; ++k )
            {
                alpha_r += mN( k ) * mTauPowJ( k ) * mF( k ) * mDeltaPowI( k )
                        * ( mI( k ) - mL( k ) * mDeltaPowL( k ) );
            }

            // ∂/∂δ( δ^i * exp( -φ*(δ-1)² - ... ) ) = F * δ^i * ( i - 2*φ*δ*(δ-1) ) / δ
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                alpha_r += mN( k ) * mTauPowJ( k ) * mF( k ) * mDeltaPowI( k )
                        * ( mI( k ) - 2. * mPhi( kk ) * mDelta * ( mDelta - 1. ) );
            }

            alpha_r /= mDelta ;

            return alpha_r;

        }

        real
        EoS_Nitrogen::compute_phir_dd() const
        {
            this->update_f();


            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mI( k ) * ( mI( k ) - 1. );
            }

            for ( uint k=6; k<32; ++k )
            {
                /* delta^2 * d2/ddelta2( delta^i * exp( -delta^l ) )
                  *      = delta^i * exp( -delta^l )
                  *      * [ i*(i-1) - 2*l*i*delta^l
                  *          + l*delta^l * ( l*delta^l - l + 1 ) ]
                  * note delta^l, not delta^i, inside the bracket */
                alpha_r += mN( k ) * mTauPowJ( k ) * mF( k ) * mDeltaPowI( k ) *
                    ( mI( k ) * ( mI( k ) - 1. ) - 2. * mL( k ) * mI( k ) * mDeltaPowL( k )
                        + mL( k ) * mDeltaPowL( k ) * ( mL( k ) * mDeltaPowL( k ) - mL( k ) + 1. ) );
            }

            alpha_r /= mDelta  * mDelta ;

            return alpha_r;
        }

        real
        EoS_Nitrogen::compute_phir_t() const
        {
            this->update_f();

            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mJ( k );
            }
            for ( uint k=6; k<32; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k ) * mJ( k );
            }
            /* tau * d/dtau( tau^j * exp( -beta*(tau-gamma)^2 ) )
              *      = tau^j * exp( ... ) * ( j - 2*beta*tau*( tau - gamma ) ) */
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * ( mJ( k ) - 2. * mBeta( kk ) * mTau * ( mTau - mGamma( kk ) ) );
            }

            alpha_r /= mTau ;
            return alpha_r;
        }

        real
        EoS_Nitrogen::compute_phir_tt() const
        {
            this->update_f();

            real alpha_r = 0. ;

            // tau^2 * d2/dtau2( tau^j ) = tau^j * j * ( j - 1 )
            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mJ( k ) * ( mJ( k ) - 1. );
            }

            // exp( -delta^l ) carries no tau, so the same factor applies
            for ( uint k=6; k<32; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * mJ( k ) * ( mJ( k ) - 1. );
            }
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                /* with  q = j - 2*beta*tau*( tau - gamma ) = tau * dln(f)/dtau,
                  * tau^2 * d2f/dtau2 / f = q^2 - j - 2*beta*tau^2 */
                const real q = mJ( k ) - 2. * mBeta( kk ) * mTau * ( mTau - mGamma( kk ) );

                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * ( q * q - mJ( k ) - 2. * mBeta( kk ) * mTau * mTau );
            }

            alpha_r /= mTau * mTau ;
            return alpha_r;
        }

        real
        EoS_Nitrogen::compute_phir_dt() const
        {
            this->update_f();

            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mI( k ) * mJ( k );
            }
            for ( uint k=6; k<32; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k ) * mJ( k ) * ( mI( k ) - mL( k ) * mDeltaPowL( k ) );
            }
            /* ∂²/∂δ∂τ( δ^i τ^j F ) = δ^i τ^j F
             *      * ( i - 2φδ(δ-1) ) * ( j - 2βτ(τ-γ) ) / ( δτ ) */
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                    * ( mI( k ) - 2. * mPhi( kk ) * mDelta * ( mDelta - 1. ) )
                    * ( mJ( k ) - 2. * mBeta( kk ) * mTau * ( mTau - mGamma( kk ) ) );
            }

            alpha_r /= mTau * mDelta ;

            return alpha_r;
        }
    }
}