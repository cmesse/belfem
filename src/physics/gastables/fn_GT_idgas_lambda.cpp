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
        idgas_lambda( RefGas * aGas, const real T )
        {
            // Method of Chung, as given in Poling, Prausnitz & O'Connell,
            // The Properties of Gases and Liquids, 5th ed., § 10-3.2
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
                real cp = aGas->Cp( T ) / M;

                // Poling Eq. ( 10-3.14 )
                real alpha = cp/R - 2.5;

                /* Poling restricts the omega correlation for beta to nonpolar
                 * species and recommends 1/1.32 as the polar default. Polarity is
                 * gauged by the same reduced dipole moment as in idgas_mu,
                 * Eq. ( 9-4.17 ), with the threshold of Eq. ( 9-4.18 ) */
                real mu_r = 52.46e-5 * std::pow(
                        aGas->data()->dipole() / T_crit, 2 ) * aGas->data()->p_crit();

                real beta  = mu_r < 0.022 ?
                        0.7862 + omega * ( 1.3168 * omega - 0.7109 ) : 1.0 / 1.32;

                real gamma = 2.0 + 10.5 * std::pow( T / T_crit, 2 );

                real Psi   = 1.0 + alpha * ( 0.215 + 0.28288*alpha
                                             - 1.061 * beta + 0.26665 * gamma ) /
                                   ( 0.6366 + beta * ( gamma + 1.061 * alpha ) );

                real tMu = 0.0;
                if( aGas->has_viscosity() )
                {
                    tMu = aGas->mu( T );
                }
                if( std::abs( tMu ) < 1e-9 )
                {
                    tMu = idgas_mu( aGas, T );
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