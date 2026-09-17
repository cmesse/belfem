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

#include "matrices/mt_maxwell_l2_phi.hpp"
#include "matrices/fn_maxwell_l2N.hpp"
#include "fn_trans.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_IWG.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            l2_phi( Calculator * aCalc, Matrix< real > & aK,  Vector< real > & aF )
            {


                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                // grab the phi data for this vector
                const Vector< real > & phi =  aCalc->node_data( "phi");

                Vector< real > & h = aCalc->vector("h");

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); k++ )
                {
                    const Matrix< real > & N = l2N( aCalc, k );
                    const Matrix< real > & B = aCalc->B( k ) ;

                    // contribution to matrix
                    aK += w( k ) * trans( N ) * N * aCalc->dV( k );

                    // actually, negative h here
                    h = -1.0 * B * phi ;

                    // contribution to right hand side
                    aF  +=  w( k ) * trans( N ) * h * aCalc->dV( k );
                }
            }

            void
            l2_phi_ferro( Calculator * aCalc, Matrix< real > & aK,  Vector< real > & aF )
            {
                // get integration weights
                const Vector< real > & w =  aCalc->integration()->weights();

                // grab the phi data for this vector
                const Vector< real > & phi =  aCalc->node_data( "phi");

                Vector< real > & h = aCalc->vector("h");

                // loop over all integration points
                for( uint k=0; k<aCalc->num_intpoints(); k++ )
                {
                    const Matrix< real > & N = l2N( aCalc, k );
                    const Matrix< real > & B = aCalc->B( k ) ;

                    // contribution to matrix
                    aK += w( k ) * trans( N ) * N * aCalc->dV( k );

                    // actually, negative h here
                    h = -1.0 * B * phi ;

                    real mu = aCalc->material()->mu( norm( h ) );

                    // contribution to right hand side
                    aF  +=  w( k ) * trans( N ) * mu * h * aCalc->dV( k );
                }
            }
        }
    }
}

