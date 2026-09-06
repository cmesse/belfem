//
// Created by gregorygiard on 10/29/25.
//

#include "mt_thermal_phi.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
        void
        T_phi( Calculator * aCalc, TimestepMatrices * aMatrices )
        {
            const Vector< real > & w =  aCalc->integration()->weights();

            // grab the material
            Material * tMaterial = aCalc->group()->material();

            real rho = tMaterial->density( gTroom );

            for( uint k=0 ; k<aCalc->num_intpoints() ; ++k )
            {

                const Matrix< real > & B = aCalc->B( k );

                // get the gradient operator (shape functions)
                const Matrix< real > & N = aCalc->N(k) ;

                real T = std::max(dot( aCalc->Nvec(k), aCalc->q() ), gTmin);

                aMatrices->M() += w( k ) * trans( N ) * rho * tMaterial->cp( T ) *  N * aCalc->dV( k );

                aMatrices->K() += w( k ) * trans( B ) * tMaterial->lambda( T ) * B * aCalc->dV( k );

            }
        }
    }
}