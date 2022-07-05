//
// Created by Christian Messe on 26.08.19.
//
#include "assert.hpp"
#include "fn_GT_idgas_lambda.hpp"
#include "fn_GT_idgas_mu.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        /**
         * calculate the conductivity in case there is no data in the database
         */
        real
        idgas_lambda( RefGas * aGas, const real aT )
        {
            // Method of Chung
            const real & R      =  aGas->data()->R();
            const real & M      =  aGas->data()->M();
            const real & T_crit =  aGas->data()->T_crit();
            const real & omega  =  aGas->data()->acentric();

            // check input
            BELFEM_ASSERT( ! std::isnan( M ),
                        "M not found in database for %s",
                          aGas->label().c_str() );

            if ( aGas->has_thermo() && aGas->data()->has_crit()  )
            {
                // evaluate cp in J/(kg*K)
                real cp = aGas->Cp( aT ) / M;

                // ( 104 b)
                real alpha = cp/R - 2.5;
                real beta  = 0.7862 + omega * ( 1.368 * omega - 0.7109 );
                real gamma = 2.0 + 10.5 * std::pow( aT / T_crit, 2 );

                // ( 104 a )
                real Psi   = 1.0 + alpha * ( 0.215 + 0.28288*alpha
                                             - 1.061 * beta + 0.26665 * gamma ) /
                                   ( 0.6366 + beta * ( gamma + 1.061 * alpha ) );

                real tMu = 0.0;
                if( aGas->has_viscosity() )
                {
                    tMu = aGas->mu( aT );
                }
                if( std::abs( tMu ) < 1e-9 )
                {
                    tMu = idgas_mu( aGas, aT );
                }

                return  3.75 * Psi * R * tMu;

            }
            else
            {
                return 0.0;
            }
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */