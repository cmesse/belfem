//
// Created by christian on 9/20/24.
//

#include "matrices/mt_maxwell_symmetry.hpp"

#include "cl_FEM_Block.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_IWG.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"
#include "fn_inv.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            symmetry_phi_2d( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->master_integration()->weights();

                uint tNumNodes = aCalc->element()->master()->element()->number_of_nodes() ;

                Matrix< real > & nxB = aCalc->matrix("nxBm");

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); ++k )
                {

                    // get the gradient operator
                    const Matrix< real > & B = aCalc->Bm( k );

                    // get the normal of the element
                    const Vector< real > & n = aCalc->normal( k );

                    for( uint i=0; i<tNumNodes; ++i )
                    {
                        nxB( 0, i ) = n( 0 ) * B( 1, i ) - n( 1 ) * B( 0, i );
                    }

                    aMatrices->K() += w( k ) * trans( nxB ) * nxB * aCalc->dS( k );
                }
            }

            void
            symmetry_phi_3d( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->master_integration()->weights();

                Matrix< real > & nxB = aCalc->matrix("nxBm");

                uint tNumNodes = aCalc->element()->master()->element()->number_of_nodes() ;

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); ++k )
                {
                    // get the gradient operator
                    const Matrix< real > & B = aCalc->Bm( k );

                    // get the normal of the element
                    const Vector< real > & n = aCalc->normal( k );

                    for( uint i=0; i<tNumNodes; ++i )
                    {
                        nxB( 0, i ) = n( 1 ) * B( 2, i ) - n( 2 ) * B( 1, i );
                        nxB( 1, i ) = n( 2 ) * B( 0, i ) - n( 0 ) * B( 2, i );
                        nxB( 2, i ) = n( 0 ) * B( 1, i ) - n( 1 ) * B( 0, i );
                    }

                    // add to integration using a least squares weighting
                    aMatrices->K() += w( k ) * trans( nxB ) * nxB * aCalc->dS( k );
                }
            }

            void
            h_symmetry_2d( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                Matrix < real  > & nxE = aCalc->matrix("nxEm");

                for( uint k=0 ; k<aCalc->num_intpoints() ; ++k )
                {
                    // edge interpolation operator from master
                    const Matrix< real > & E = aCalc->Em( k );

                    // normal
                    const Vector< real > & n = aCalc->normal( k );

                    nxE.fill( 0.0 );
                    for ( uint i=0; i<aCalc->num_nedelec_dofs(); ++i )
                    {
                        nxE( 0, i ) = n( 0 ) * E( 1, i ) - n( 1 ) * E( 0, i );

                    }
                    aMatrices->K() += trans( nxE ) * nxE * ( w( k ) *  aCalc->dS( k ) );
                }
            }


            void
            h_symmetry_3d( Calculator * aCalc, TimestepMatrices * aMatrices )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                Matrix < real  > & nxE = aCalc->matrix("nxEm");

                for( uint k=0 ; k<aCalc->num_intpoints() ; ++k )
                {
                    // edge interpolation operator from master
                    const Matrix< real > & E = aCalc->Em( k );

                    // normal
                    const Vector< real > & n = aCalc->normal( k );

                    nxE.fill( 0.0 );
                    for ( uint i=0; i<aCalc->num_nedelec_dofs(); ++i )
                    {
                        nxE( 0, i ) =   n(1) * E( 2, i )-n(2) * E( 1, i ) ;
                        nxE( 1, i ) =   n(2) * E( 0, i )-n(0) * E( 2, i ) ;
                        nxE( 2, i ) =   n(0) * E( 1, i )-n(1) * E( 0, i ) ;

                    }
                    aMatrices->K() += trans( nxE ) * nxE * ( w( k ) *  aCalc->dS( k ) );
                }
            }

        }
    }
}