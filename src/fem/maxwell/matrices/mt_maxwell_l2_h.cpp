//
// Created by christian on 10/24/24.
//
#include "matrices/mt_maxwell_l2_phi.hpp"
#include "matrices/fn_maxwell_l2N.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            l2_h( Calculator * aCalc, Matrix< real > & aK,  Vector< real > & aF )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                // grab the phi data for this vector
                const Vector< real > & h =  aCalc->nedelec_data_h() ;

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); k++ )
                {

                    const Matrix< real > & E = aCalc->E( k ) ;
                    const Matrix< real > & N = l2N( aCalc, k );

                    // contribution to striffness
                    aK += w( k ) * trans( N ) * N * aCalc->dV( k );

                    // contribution to right hand side
                    aF  +=  w( k ) * trans( N ) * E * h * aCalc->dV( k );
                }
            }

//------------------------------------------------------------------------------

            void
            l2_ah_2d( Calculator * aCalc, Matrix< real > & aK,  Vector< real > & aF )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                // grab the phi data for this vector
                const Vector< real > & Az =  aCalc->node_data( "a_z");

                Matrix< real > & C = aCalc->matrix("HelpA");
                Vector< real > & b = aCalc->vector("HelpB");

                uint n = aCalc->element()->element()->number_of_nodes() ;

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); k++ )
                {
                    const Matrix< real > & N = l2N( aCalc, k );

                    const Matrix< real > & B = aCalc->B( k ) ;

                    // assemble curl operator
                    for ( uint i=0; i<n; ++i )
                    {
                        C( 0, i ) =  B( 1, i );
                        C( 1, i ) = -B( 0, i );
                    }

                    // contribution to matrix
                    aK += w( k ) * trans( N ) * N * aCalc->dV( k );

                    b = C * Az ;

                    real nu = aCalc->material()->nu( norm( b ) );

                    // contribution to right hand side
                    aF  +=  w( k ) * trans( N ) * b * nu * aCalc->dV( k );
                }
            }

//------------------------------------------------------------------------------
        }
    }
}