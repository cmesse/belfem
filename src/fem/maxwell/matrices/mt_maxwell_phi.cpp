//
// Created by christian on 9/19/24.
//
#include "constants.hpp"
#include "cl_FEM_Element.hpp"
#include "matrices/mt_maxwell_phi.hpp"

#include "cl_FEM_Controller.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_TimestepMatrices.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
//------------------------------------------------------------------------------

            void
            phi( Calculator * aCalc,
                 TimestepMatrices * aMatrices)
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); ++k )
                {
                    // get the gradient operator
                    const Matrix< real > & B = aCalc->B( k );


                    // add to integration
                    aMatrices->M() += w( k ) * trans( B ) *  B * aCalc->dV( k );
                }


                aMatrices->M() *= constant::mu0 ;
            }

//------------------------------------------------------------------------------

            void
            phi_ferro_picard( Calculator * aCalc,
                 TimestepMatrices * aMatrices)
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                real mu = constant::mu0 ;
                real dmudH = 0.0 ;

                const Vector< real > phi = aCalc->node_data( "phi") ;

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); ++k )
                {

                    const Matrix< real > & B = aCalc->B( k );

                    // compute magnetic field at current timestep
                    real h = norm( B * phi );


                    // compute material properties and derivatives
                    aCalc->material()->dmudH( h, mu, dmudH );

                    // add to integration
                    aMatrices->M() += w( k ) * trans( B ) * mu * B * aCalc->dV( k );
                }
            }

            void
            phi_ferro_newton( Calculator * aCalc,
                 TimestepMatrices * aMatrices)
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                real mu = constant::mu0 ;
                real dmudh = 0.0 ;

                const Vector< real > phi = aCalc->node_data( "phi") ;

                // beta-weighted dof history of the ACTIVE timestepping scheme
                // ( == phi0 for order <= 1 ). The residual contracts M against
                // exactly this combination, so using it here — instead of the
                // old phi0-only shortcut — keeps the dMdX_times_h block a
                // consistent tangent for every BDF order, ramp state and
                // variable step size ( single source of truth: collect_qhist )
                const Vector< real > & phi_hist = static_cast< IWG_Timestep * >(
                    aCalc->group()->parent()->iwg() )->collect_qhist();

                // scratch for the intpoint field vectors ( registered in
                // IWG_Maxwell::create_custom_vectors_and_matrices );
                // the sign of h = -grad(phi) is dropped: h, and the
                // contraction below, are invariant under joint negation
                Vector< real > & hvec      = aCalc->vector( "hcur" );
                Vector< real > & hvec_hist = aCalc->vector( "hhist" );

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); ++k )
                {

                    const Matrix< real > & B = aCalc->B( k );

                    // compute magnetic field at current timestep
                    hvec = B * phi ;
                    real h = norm( hvec );

                    // compute material properties and derivatives
                    aCalc->material()->dmudH( h, mu, dmudh );

                    // add to integration
                    aMatrices->M() += w( k ) * trans( B ) * mu * B * aCalc->dV( k );

                    if ( std::abs( dmudh ) > BELFEM_EPSILON )
                    {
                        hvec_hist = B * phi_hist ;

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
                        aMatrices->dMdx_times_x() += w( k ) * trans( B ) * (dmudh*h) * B * aCalc->dV( k );
                        aMatrices->dMdx_times_h() += w( k ) * trans( B ) * (dmudh*h0) * B * aCalc->dV( k );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            phi_tri3( Calculator * aCalc,
                 TimestepMatrices * aMatrices)
            {
                // get the gradient operator
                const Matrix< real > & B = aCalc->B( 0 );

                aMatrices->M() = trans( B ) * B * ( 0.5 * constant::mu0 * aCalc->dV() );
            }

//------------------------------------------------------------------------------

            void
            phi_tet4( Calculator * aCalc,
                 TimestepMatrices * aMatrices)
            {
                // get the gradient operator
                const Matrix< real > & B = aCalc->B( 0 );

                aMatrices->M() = trans( B ) * B * ( 1./6. * constant::mu0 * aCalc->dV() );
            }

//------------------------------------------------------------------------------

          void
            phi_tri6_tet10( Calculator * aCalc,
                   TimestepMatrices * aMatrices)
          {
                if( aCalc->element()->element()->is_curved() )
                {
                    phi( aCalc, aMatrices );
                }
                else
                {
                    // get integration weights
                    const Vector< real > & w =  aCalc->integration()->weights();

                    // loop over all integration points
                    for( uint k=0; k<aCalc->num_intpoints(); ++k )
                    {
                        // get the gradient operator
                        const Matrix< real > & B = aCalc->B( k );

                        // add to integration
                        aMatrices->M() += w( k ) * trans( B ) *  B ;
                    }

                    aMatrices->M() *= constant::mu0 * aCalc->dV() ;
                }
          }

//------------------------------------------------------------------------------
        }
    }
}
