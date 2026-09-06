//
// Created by christian on 12/4/24.
//
#include "matrices/mt_maxwell_l2_b.hpp"
#include "matrices/fn_maxwell_l2N.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            l2_b2d( Calculator * aCalc, Matrix< real > & aK,  Vector< real > & aF )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                // grab the phi data for this vector
                const Vector< real > & Az =  aCalc->node_data( "a_z" );

                Matrix< real > & C = aCalc->matrix("HelpA");

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

                    // contribution to right hand side
                    aF  +=  w( k ) * trans( N ) * C * Az * aCalc->dV( k );
                }
            }
        }
    }
}
