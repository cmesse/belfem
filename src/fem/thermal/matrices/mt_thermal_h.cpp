//
// Created by gregorygiard on 10/29/25.
//
#include "mt_thermal_h.hpp"
#include "fn_trans.hpp"
#include "cl_IWG_Timestep.hpp"
namespace belfem
{
    namespace fem
    {
        void
        T_h_picard( Calculator * aCalc, TimestepMatrices * aMatrices )
        {
            const Vector< real > & w =  aCalc->integration()->weights();

            // material helper ( thermal-side instance ): E/C/b/j come from the
            // paired Maxwell calculator, T from this element; the peer element
            // is already linked by Calculator::link. Covers metal, alloy, HTS,
            // defect and piecewise variants via constructor-time dispatch.
            calculator::MaxwellData * mx = aCalc->maxwell();

            for( uint k=0 ; k<aCalc->num_intpoints() ; ++k )
            {
                // gradient operator
                const Matrix< real > & B = aCalc->B( k );

                // shape function matrix
                const Matrix< real > & N = aCalc->N( k ) ;

                // specific heat at the ( clamped ) local temperature
                real cp = mx->compute_cp( k );

                // thermal conductivity ( field-dependence routed on
                // depends( lambda, normB ) inside the helper )
                real lambda = mx->compute_lambda( k );

                // resistivity for the Joule source, clamped [ gRhoMin, gRhoMax ]
                real rho = mx->compute_rho( k );

                // current density magnitude from the Maxwell peer
                real norm_j = mx->norm_j( k );

                // density: reference value of the undeformed mesh
                // ( see the physics-trap note at MaxwellData::mDensity )
                aMatrices->M() += w( k ) * trans( N ) * mx->density() * cp * N * aCalc->dV( k );
                aMatrices->K() += w( k ) * trans( B ) * lambda * B * aCalc->dV( k );

                // Joule heating source plus the artificial volumetric heat
                // load ( W/m³, zero unless the material carries a heating plugin )
                real dotq = rho * norm_j * norm_j + mx->compute_volumetric_heatload( k );

                aMatrices->f() += w( k ) * trans( N ) * dotq * aCalc->dV( k );
            }
        }

        void
        T_h_newton( Calculator * aCalc, TimestepMatrices * aMatrices )
        {
            const Vector< real > & w =  aCalc->integration()->weights();

            // material helper ( thermal-side instance ): E/C/b/j come from the
            // paired Maxwell calculator, T from this element; the peer element
            // is already linked by Calculator::link. Covers metal, alloy, HTS,
            // defect and piecewise variants via constructor-time dispatch.
            calculator::MaxwellData * mx = aCalc->maxwell();

            // beta-weighted dof history of the ACTIVE timestepping scheme;
            // the residual contracts M against exactly this combination
            // ( single source of truth: collect_qhist, see phi_ferro )
            const Vector< real > & qhist = static_cast< IWG_Timestep * >(
                aCalc->group()->parent()->iwg() )->collect_qhist();

            // n x 1 workspace for the conductivity tangent ( bound once,
            // no per-point matrix construction )
            Matrix< real > & Btg = aCalc->matrix( "Btg" );

            for( uint k=0 ; k<aCalc->num_intpoints() ; ++k )
            {
                // gradient operator
                const Matrix< real > & B = aCalc->B( k );

                // shape function matrix
                const Matrix< real > & N = aCalc->N( k ) ;

                // specific heat at the ( clamped ) local temperature
                real cp = mx->compute_cp( k );
                real dcpdT = mx->compute_dcpdT( k );

                // thermal conductivity ( field-dependence routed on
                // depends( lambda, normB ) inside the helper )
                real lambda = mx->compute_lambda( k );
                real dlambdadT = mx->compute_dlambdadT( k );

                // resistivity for the Joule source, clamped [ gRhoMin, gRhoMax ]
                real rho = mx->compute_rho( k );
                real drhodT = mx->compute_drhodT( k );

                // current density magnitude from the Maxwell peer
                real norm_j = mx->norm_j( k );

                real T = mx->compute_T( k );

                // signed history contraction ( exact scalar N·qhist, not a
                // magnitude — cf. thermal_matrices plan §3, rows 1-2 )
                real T0 = dot( aCalc->Nvec( k ), qhist );


                real wdV = w( k ) * aCalc->dV( k );

                // density: reference value of the undeformed mesh
                // ( see the physics-trap note at MaxwellData::mDensity )
                aMatrices->M() += trans( N ) * N * ( mx->density() * cp * wdV );
                aMatrices->K() += trans( B ) * B * ( lambda * wdV );

                aMatrices->dMdx_times_x() += trans( N ) * N * ( mx->density() * dcpdT * T * wdV );
                aMatrices->dMdx_times_h() += trans( N ) * N * ( mx->density() * dcpdT * T0 * wdV );

                // conductivity change: exact mixed-operator tangent
                // ( Bᵀ·∇T ) ⊗ ( dλ/dT · N ) — lambda depends on T = N·q,
                // K is built from B, so the block is unsymmetric
                Btg = trans( B ) * ( B * aCalc->q() );
                aMatrices->dKdx_times_x() += Btg * N * ( dlambdadT * wdV );

                // Joule heating source and the quench-feedback tangent
                // dFdX = ∂f_i/∂x_j = Nᵀ·( dρ/dT·|j|² )·N ( no contraction;
                // sign handled by assemble_dJdx: -Δt·dFdX )
                // The artificial heat load is prescribed in ( x, t ) only,
                // so it enters f and has no tangent
                real dotq = rho * norm_j * norm_j + mx->compute_volumetric_heatload( k );

                aMatrices->f() += trans( N ) * dotq * wdV ;
                aMatrices->dfdx() += trans( N ) * N * ( drhodT * norm_j * norm_j * wdV );
            }
        }

    }
}
