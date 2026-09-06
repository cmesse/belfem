//
// Created by christian on 10/23/24.
//
#include "constants.hpp"
#include "globals.hpp"
#include "mt_maxwell_h.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Controller.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "petsctools.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {

            /**
             * adds the rho(|B|,beta) field-derivative tangent row for
             * field-dependent normal conductors ( Kohler metals, alloy
             * tables ): dKdx += ( C^T j ) (x) Row with
             *
             *   Row = drho/d|B| * ( mu + |h| dmu/dh ) * hhat^T E
             *       + drho/dbeta * dbeta/dq
             *
             *   dbeta/dq = -sign(c)/sqrt(1-c^2) *
             *              [ (1/|h|) ( jhat - c bhat )^T E
             *              + (1/|j|) ( bhat - c jhat )^T C ]
             *
             * ( derivation + Codex verification 2026-07-21, exchange thread
             * maxwell_h_kernels_audit ). The beta term exists only in 3D
             * ( 2D: beta = pi/2 by convention ) and is dropped near the
             * beta = 0 / pi/2 kinks of the abs() in bj_angle; the whole row
             * is dropped at the small-field cutoffs of bj_angle.
             */
            static void
            add_rho_field_tangent(
                Calculator * aCalc,
                calculator::MaxwellData * mx,
                TimestepMatrices * aMatrices,
                const uint k,
                const real wdV )
            {
                real drdB    = mx->compute_drhodb( k );
                real drdbeta = mx->compute_drhodbeta( k );

                if ( drdB == 0.0 && drdbeta == 0.0 ) return ;

                const Vector< real > & b = mx->compute_b( k );
                const Vector< real > & j = mx->compute_j( k );

                real B = mx->norm_b( k );
                real J = mx->norm_j( k );

                // small-field guards, consistent with bj_angle's fallback
                if ( B < 1e-6 || J < 1e-6 ) return ;

                // variable-mu chain: |B| = mu(|h|)*|h|
                // -> d|B|/dq = ( mu + |h| dmudh ) hhat^T E ; dmudh = 0 for
                // constant mu, so this covers h_newton_mu0 exactly as well
                real mu    = mx->compute_mu( k );
                real dmudh = mx->compute_dmudh( k );
                real H     = B / mu ;   // |h|, since b is parallel to h

                Vector< real > & vE = aCalc->vector( "vE" );
                Vector< real > & vC = aCalc->vector( "vC" );

                // |B| channel ( bhat = b/B )
                vE = ( drdB * ( mu + H * dmudh ) / B ) * b ;

                bool tHaveBetaC = false ;

                if ( b.length() == j.length() && drdbeta != 0.0 )
                {
                    real c = dot( b, j ) / ( B * J );
                    real absc = std::abs( c );

                    if ( absc > 1e-8 && 1.0 - absc > 1e-8 )
                    {
                        real dbdc = -( c > 0.0 ? 1.0 : -1.0 )
                                  / std::sqrt( 1.0 - c*c );

                        // E channel of beta: (1/|h|) ( jhat - c bhat )
                        vE += ( drdbeta * dbdc / H )
                            * ( ( 1.0/J ) * j - ( c/B ) * b );

                        // C channel of beta: (1/|j|) ( bhat - c jhat )
                        vC = ( drdbeta * dbdc / J )
                           * ( ( 1.0/B ) * b - ( c/J ) * j );

                        tHaveBetaC = true ;
                    }
                }

                Matrix< real > & Ctj = aCalc->matrix( "Ctj" );
                Ctj = trans( aCalc->C( k ) ) * j ;

                Matrix< real > & Rw = aCalc->matrix( "Rw" );
                Rw = trans( aCalc->E( k ) ) * vE ;

                if ( tHaveBetaC )
                {
                    Rw += trans( aCalc->C( k ) ) * vC ;
                }

                // dt scaling happens in assemble_dJdx
                aMatrices->dKdx_times_x() += Ctj * trans( Rw ) * wdV ;
            }

            void
            h_picard( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                const Vector< real > & w = aCalc->integration()->weights();
                calculator::MaxwellData * mx = aCalc->maxwell();

                real rho_m = 0.0;
                real V     = 0.0;
                real dotQ  = 0.0 ;

                const real chi = aCalc->group()->parent()->iwg()->penalty( 2 ) ;
                bool tUseGauging = chi > BELFEM_EPSILON ;

                for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
                {
                    const Matrix< real > & E = aCalc->E( k );
                    const Matrix< real > & C = aCalc->C( k );

                    const Vector< real > & j = mx->compute_j( k );
                    real mu   = mx->compute_mu( k );
                    real rho  = mx->compute_rho( k );
                    real wdV = w( k ) * aCalc->dV( k );

                    rho_m += rho * wdV;
                    V     += wdV;
                    dotQ  += rho * dot( j, j ) * wdV ;

                    aMatrices->M() += trans( E ) * E * ( mu * wdV );
                    aMatrices->K() += trans( C ) * C * ( rho * wdV );

                    // gauging
                    if ( tUseGauging )
                    {
                        const Matrix< real > & G = aCalc->G( k );
                        aMatrices->K() += trans( G ) * G * ( chi * rho * wdV );
                    }

                }
                save_resistivity( aCalc, rho_m / V );
                save_dotQ( aCalc, dotQ );
            }

            void
            h_newton_mu0( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                const Vector< real > & w = aCalc->integration()->weights();
                calculator::MaxwellData * mx = aCalc->maxwell();

                real rho_m = 0.0;
                real V     = 0.0;
                real dotQ  = 0.0 ;

                const real chi = aCalc->group()->parent()->iwg()->penalty( 2 ) ;
                bool tUseGauging = chi > BELFEM_EPSILON ;

                Matrix< real > & Ctj = aCalc->matrix( "Ctj" );

                // q() gathers the dofs on every call — bind once per element
                const Vector< real > & q = aCalc->q();
                Vector< real > & Gq   = aCalc->vector( "Gq" );
                Matrix< real > & GtGq = aCalc->matrix( "GtGq" );

                for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
                {
                    const Matrix< real > & E = aCalc->E( k );
                    const Matrix< real > & C = aCalc->C( k );

                    const Vector< real > & j = mx->compute_j( k );
                    real mu   = mx->compute_mu( k );
                    real rho  = mx->compute_rho( k );
                    real drho = mx->compute_drhodj( k );

                    real wdV = w( k ) * aCalc->dV( k );
                    rho_m += rho * wdV;
                    V     += wdV;
                    dotQ  += rho * dot( j, j ) * wdV ;

                    aMatrices->M() += trans( E ) * E * ( mu * wdV );   // per-term mu placement (T2)
                    aMatrices->K() += trans( C ) * C * ( rho * wdV );

                    if ( mx->norm_j( k ) > BELFEM_EPSILON )
                    {
                        const Vector< real > & j = mx->compute_j( k );
                        Ctj = trans( C ) * j ;
                        aMatrices->dKdx_times_x() += Ctj * trans( Ctj )
                                                   * ( ( drho / mx->norm_j( k ) ) * wdV );
                    }

                    // gauging
                    if ( tUseGauging )
                    {
                        const Matrix< real > & G = aCalc->G( k );
                        aMatrices->K() += trans( G ) * G * ( chi * rho * wdV );

                        // the gauge weight is rho(|J(q)|), so the tangent
                        // carries a drho channel as well; Ctj is valid here
                        // because the power-law block above shares the guard
                        if ( mx->norm_j( k ) > BELFEM_EPSILON )
                        {
                            Gq   = G * q ;
                            GtGq = trans( G ) * Gq ;
                            aMatrices->dKdx_times_x() += GtGq * trans( Ctj )
                                * ( ( chi * drho / mx->norm_j( k ) ) * wdV );
                        }
                    }

                    // rho(|B|,beta) field-derivative channel ( metals/alloys )
                    add_rho_field_tangent( aCalc, mx, aMatrices, k, wdV );
                }
                save_resistivity( aCalc, rho_m / V );
                save_dotQ( aCalc, dotQ );
            }

            void
            h_newton_mu( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                const Vector< real > & w = aCalc->integration()->weights();
                calculator::MaxwellData * mx = aCalc->maxwell();

                const real chi = aCalc->group()->parent()->iwg()->penalty( 2 ) ;
                bool tUseGauging = chi > BELFEM_EPSILON ;

                real rho_m = 0.0;
                real V     = 0.0;
                real dotQ  = 0.0 ;

                Matrix< real > & Ctj = aCalc->matrix( "Ctj" );

                // q() gathers the dofs on every call — bind once per element
                const Vector< real > & q = aCalc->q();
                Vector< real > & Gq   = aCalc->vector( "Gq" );
                Matrix< real > & GtGq = aCalc->matrix( "GtGq" );

                // beta-weighted dof history of the ACTIVE timestepping scheme;
                // the residual contracts M against exactly this combination
                // ( single source of truth: collect_qhist, see phi_ferro )
                const Vector< real > & qhist = static_cast< IWG_Timestep * >(
                    aCalc->group()->parent()->iwg() )->collect_qhist();

                // scratch for the history field at the intpoint ( registered
                // in IWG_Maxwell::create_custom_vectors_and_matrices )
                Vector< real > & hvec_hist = aCalc->vector( "hhist" );

                for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
                {
                    const Matrix< real > & E = aCalc->E( k );
                    const Matrix< real > & C = aCalc->C( k );

                    real dmudh  = mx->compute_dmudh( k );
                    real mu   = mx->compute_mu( k );

                    const Vector< real > & j = mx->compute_j( k );
                    real drho = mx->compute_drhodj( k );
                    real rho  = mx->compute_rho( k );


                    real wdV = w( k ) * aCalc->dV( k );
                    rho_m += rho * wdV;
                    V     += wdV;
                    dotQ  += rho * dot( j, j ) * wdV ;

                    aMatrices->M() += trans( E ) * E * ( mu * wdV );   // per-term mu placement (T2)
                    aMatrices->K() += trans( C ) * C * ( rho * wdV );

                    if ( std::abs( dmudh ) > BELFEM_EPSILON )
                    {
                        // current-iterate field as seen by mu(|h|)
                        const Vector< real > & hvec = mx->compute_h( k );
                        real h = norm( hvec );

                        hvec_hist = E * qhist ;

                        // signed history contraction hhat·hvec_hist ( the
                        // exact isotropized scalar, not a magnitude — the
                        // magnitude flips the tangent sign under field
                        // reversal; cf. mt_thermal_h.cpp ). At h -> 0 the
                        // direction hhat is undefined: take the zero
                        // subgradient and drop the history term — the
                        // dMdx_times_x term already vanishes with h, so the
                        // tangent degrades to the Picard matrix there
                        real h0 = h > BELFEM_EPSILON ?
                                dot( hvec, hvec_hist ) / h : 0.0 ;

                        // isotropized nonlinear-mass tangent blocks; assembled as
                        // dJdx += alpha * dMdx_times_x - dMdx_times_h
                        aMatrices->dMdx_times_x() += trans( E ) * E * ( dmudh * h * wdV );
                        aMatrices->dMdx_times_h() += trans( E ) * E * ( dmudh * h0 * wdV );
                    }

                    if ( mx->norm_j( k ) > BELFEM_EPSILON )
                    {
                        const Vector< real > & j = mx->compute_j( k );
                        Ctj = trans( C ) * j ;
                        aMatrices->dKdx_times_x() += Ctj * trans( Ctj )
                                                   * ( ( drho / mx->norm_j( k ) ) * wdV );
                    }

                    // gauging
                    if ( tUseGauging )
                    {
                        const Matrix< real > & G = aCalc->G( k );
                        aMatrices->K() += trans( G ) * G * ( chi * rho * wdV );

                        // the gauge weight is rho(|J(q)|), so the tangent
                        // carries a drho channel as well; Ctj is valid here
                        // because the power-law block above shares the guard
                        if ( mx->norm_j( k ) > BELFEM_EPSILON )
                        {
                            Gq   = G * q ;
                            GtGq = trans( G ) * Gq ;
                            aMatrices->dKdx_times_x() += GtGq * trans( Ctj )
                                * ( ( chi * drho / mx->norm_j( k ) ) * wdV );
                        }
                    }

                    // rho(|B|,beta) field-derivative channel ( metals/alloys )
                    add_rho_field_tangent( aCalc, mx, aMatrices, k, wdV );
                }
                save_resistivity( aCalc, rho_m / V );
                save_dotQ( aCalc, dotQ );
            }


            void
            h_side_connector( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                const Vector< real > & w = aCalc->integration()->weights();
                calculator::MaxwellData * mx = aCalc->maxwell();

                const real chi = aCalc->group()->parent()->iwg()->penalty( 2 ) ;
                bool tUseGauging = chi > BELFEM_EPSILON ;

                real dotQ  = 0.0 ;
                real rho_m = 0.0;
                real V     = 0.0;

                for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
                {
                    const Matrix< real > & E = aCalc->E( k );
                    const Matrix< real > & C = aCalc->C( k );

                    const Vector< real > & j = mx->compute_j( k );
                    real mu   = mx->compute_mu( k );
                    real rho  = mx->compute_rho( k );

                    real wdV = w( k ) * aCalc->dV( k );
                    rho_m += rho * wdV;
                    V     += wdV;
                    dotQ  += rho * dot( j, j ) * wdV ;

                    aMatrices->M() += trans( E ) * E * ( mu * wdV );
                    aMatrices->K() += trans( C ) * C * ( rho * wdV );

                    // gauging
                    if ( tUseGauging )
                    {
                        const Matrix< real > & G = aCalc->G( k );
                        aMatrices->K() += trans( G ) * G * ( chi * rho * wdV );
                    }
                }
                save_resistivity( aCalc, rho_m / V );
                save_dotQ( aCalc, dotQ );
            }

            void
            h_side_connector_newton( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                const Vector< real > & w = aCalc->integration()->weights();
                calculator::MaxwellData * mx = aCalc->maxwell();

                const real chi = aCalc->group()->parent()->iwg()->penalty( 2 ) ;
                bool tUseGauging = chi > BELFEM_EPSILON ;

                real rho_m = 0.0;
                real V     = 0.0;
                real dotQ  = 0.0 ;

                for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
                {
                    const Matrix< real > & E = aCalc->E( k );
                    const Matrix< real > & C = aCalc->C( k );

                    const Vector< real > & j = mx->compute_j( k );
                    real mu   = mx->compute_mu( k );
                    real rho  = mx->compute_rho( k );

                    real wdV = w( k ) * aCalc->dV( k );
                    rho_m += rho * wdV;
                    V     += wdV;
                    dotQ  += rho * dot( j, j ) * wdV ;

                    aMatrices->M() += trans( E ) * E * ( mu * wdV );
                    aMatrices->K() += trans( C ) * C * ( rho * wdV );

                    // gauging
                    if ( tUseGauging )
                    {
                        const Matrix< real > & G = aCalc->G( k );
                        aMatrices->K() += trans( G ) * G * ( chi * rho * wdV );
                    }

                    // rho(|B|,beta) field-derivative channel; the wall is a
                    // pure metal, so there is no j-channel ( drho/dj = 0 )
                    add_rho_field_tangent( aCalc, mx, aMatrices, k, wdV );

                }
                save_resistivity( aCalc, rho_m / V );
                save_dotQ( aCalc, dotQ );
            }



            void
            h_ghost( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                mesh::Facet * tFacet = aCalc->element()->facet();

                BELFEM_ASSERT( mesh::interpolation_order(  tFacet->element()->type() ) == InterpolationOrder::LINEAR, "second and higher order elements are not supported" );

                // get the blocks
                Block * tMblock = aCalc->group()->parent()->block( tFacet->master()->block_id() ) ;
                Block * tSblock = aCalc->group()->parent()->block( tFacet->slave()->block_id() ) ;

                // the material properties don't need to be exact
                // the order of magnitude is sufficient
                // get the approximate temperature
                real T = gTbulk;

                // check if this is a thermal problem, if so, use the average temperature on the facet
                if ( std::isnan( T ) )
                {
                    T = 0.0 ;
                    const Vector< real > & Tfield = aCalc->group()->parent()->mesh()->field_data( "T" );
                    for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                    {
                        T += Tfield(tFacet->node(k)->index());
                    }
                    T /= tFacet->number_of_corner_nodes();
                }


                // thicknesses
                real hm = tMblock->block()->thickness() ;
                real hs = tSblock->block()->thickness() ;

                // approximate resistivities
                real rho_m  = get_resistivity( aCalc, tFacet->master()->index() );

                real rho_s = get_resistivity( aCalc, tFacet->slave()->index() );

                // harmonic mean of resistivities
                real rho_harm = rho_m + rho_s < BELFEM_EPSILON ? 0.0 : 2.0 * rho_m * rho_s / ( rho_m + rho_s ) ;


                // penalty parameter
                real eta = aCalc->group()->parent()->iwg()->penalty( 0 ) ;
                const real k_reg = aCalc->group()->parent()->iwg()->penalty( 1 ) ;   // Ohm, deck key: nitsche ghost penalty { k_reg }

                // -----------------------------------------------------------
                // Penalty coefficient: regularized harmonic mean of the
                // per-layer stiffnesses k = rho/h.
                //
                // The plain harmonic mean  2 km ks / ( km + ks )  is the
                // standard heterogeneous-DG choice (Burman & Zunino 2006).
                // It is bounded above by  2 * min( km, ks ),  which keeps
                // alpha well behaved when one layer is an insulator
                // (rho -> infinity, k -> infinity). But it collapses to
                // zero when one layer is a superconductor (rho -> 0).
                //
                // Adding a small offset k_reg to each stiffness before the
                // harmonic mean turns both singular limits into bounded
                // smooth limits:
                //   - k_reg -> 0          recovers Burman-Zunino
                //   - both k tiny         alpha ~ 2 k_reg (no collapse)
                //   - one k huge          alpha ~ 2 * min( other_k, k_reg )
                //                                                (no blow-up)
                //
                // k_reg should be chosen as a small fraction of the typical
                // "real" conductor stiffness for the problem class. For HTS
                // tapes at 77 K (Cu / Ag / Hastelloy at the few-micron scale)
                // k_conductor ~ 1e-3 .. 1e-2 Ohm, so k_reg = 1e-3 sits at
                // the low end of the conductor range: invisible for
                // metal-metal interfaces, large enough to keep YBCO from
                // collapsing, small enough to be dominated by metal physics
                // wherever metal physics is real.
                // -----------------------------------------------------------

                real km = ( hm > BELFEM_EPSILON ) ? rho_m / hm : 0.0 ;
                real ks = ( hs > BELFEM_EPSILON ) ? rho_s / hs : 0.0 ;

                real km_reg = km + k_reg ;
                real ks_reg = ks + k_reg ;

                real k_pen = 2.0 * km_reg * ks_reg / ( km_reg + ks_reg );

                real alpha = eta * k_pen ;

                // get the work matrices
                Matrix< real > & Kmm = aCalc->matrix("K++");
                Matrix< real > & Kms = aCalc->matrix("K+-");
                Matrix< real > & Ksm = aCalc->matrix("K-+");
                Matrix< real > & Kss = aCalc->matrix("K--");

                Matrix< real > & Dm = aCalc->matrix("D+");
                Matrix< real > & Ds = aCalc->matrix("D-");

                // reset matrices
                Kmm.fill( 0.0 );
                Kms.fill( 0.0 );
                Ksm.fill( 0.0 );
                Kss.fill( 0.0 );
                Dm.fill( 0.0 );
                Ds.fill( 0.0 );

                const Vector< real > & w = aCalc->integration()->weights();


                uint ng = w.length();

                // if you want other elements, adapt the numbers below and the Dm/Ds assembly
                ElementType tMasterType = tFacet->master()->type() ;
                uint n = mesh::number_of_nedelec_dofs( tMasterType ) ; // number of dofs per ts element
                uint m = mesh::number_of_nedelec_dofs( tFacet->element()->type() ) ;   // number of dofs on facet
                uint d = mesh::dimension( tMasterType ) ;  /// number of dimensions

                for ( uint k=0; k<ng; ++k )
                {

                    const Matrix< real > & Em = aCalc->Em( k );
                    const Matrix< real > & Es = aCalc->Es( k );

                    for ( uint i=0; i<d; ++i )
                    {
                        for ( uint j=0; j<m; ++j )
                        {
                            Dm( i, j )   = -Em( i, j+m );
                            Ds( i, j )   =  Es( i, j );
                        }
                        for ( uint j=0; j<m; ++j )
                        {
                            Dm( i, m+j ) =  Em( i, j+m );
                            Ds( i, m+j ) = -Es( i, j );
                        }
                    }
                    Dm /= hm ;
                    Ds /= hs ;

                    real wdS = w( k ) * aCalc->dS( k );
                    real awdS = alpha * wdS ;
                    real rwdS = rho_harm * wdS ;
                    Kmm +=  trans( Em ) * Em * awdS - trans(Em) * Dm * rwdS - trans( Dm ) * Em * rwdS ;
                    Kms -=  trans( Em ) * Es * awdS + trans(Em) * Ds * rwdS - trans( Dm ) * Es * rwdS ;
                    Ksm -=  trans( Es ) * Em * awdS - trans(Es) * Dm * rwdS + trans( Ds ) * Em * rwdS ;
                    Kss +=  trans( Es ) * Es * awdS + trans(Es) * Ds * rwdS + trans( Ds ) * Es * rwdS ;
                }


                Matrix< real > & K = aMatrices->K();
                for ( uint j=0; j<n; ++j )
                {
                    for ( uint i = 0; i < n; ++i )
                    {
                        K(i,j)     = Kmm( i, j );
                    }
                    for ( uint i = 0; i < n; ++i )
                    {
                        K(n+i,j)    = Ksm( i, j );
                    }
                }
                for ( uint j=0; j<n; ++j )
                {
                    for ( uint i = 0; i < n; ++i )
                    {
                        K(i,n+j) = Kms( i, j );
                    }
                    for ( uint i = 0; i < n; ++i )
                    {
                        K(n+i,n+j) = Kss( i, j );
                    }
                }
            }

        }
    }
}
