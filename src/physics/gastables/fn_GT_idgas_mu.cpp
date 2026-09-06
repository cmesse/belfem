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

#include <cmath>
#include <algorithm>
#include "assert.hpp"
#include "fn_GT_idgas_mu.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        real
        idgas_mu( RefGas * aGas, const real T )
        {
            // get refs
            const real & M      =  aGas->data()->M();
            const real & T_crit =  aGas->data()->T_crit();
            const real & p_crit =  aGas->data()->p_crit();
            const real & Z_crit =  aGas->data()->Z_crit();

            // named dipole, not mu: in this module mu is the viscosity
            const real & dipole =  aGas->data()->dipole();

            if( aGas->data()->has_crit() )
            {
                // Method of Lucas, as given in Poling, Prausnitz & O'Connell,
                // The Properties of Gases and Liquids, 5th ed., § 9-4.3

                // reduced temperature
                real T_r = T / T_crit;

                // reduced dipole moment, Eq. ( 9-4.17 )
                real mu_r = 52.46e-5 * std::pow( dipole / T_crit, 2 ) * p_crit;

                // polarity correction factor, Eq. ( 9-4.18 )
                real Fp_id;

                if ( mu_r < 0.022 )
                {
                    Fp_id = 1.0;
                }
                else if ( mu_r < 0.075 )
                {
                    Fp_id = 1.0 + 30.55 * std::pow( std::max( 0.292 - Z_crit, 0.0 ), 1.72 );
                }
                else
                {
                    Fp_id = 1.0 + 30.55 * std::pow( std::max( 0.292 - Z_crit, 0.0 ), 1.72 )
                                  * std::abs( 0.96 + 0.1 * ( T_r - 0.7 ));
                }

                /* quantum gas correction factor, Eq. ( 9-4.19 ); Lucas gives
                 * Q for helium, hydrogen and deuterium only. The exponent 1/M
                 * expects M in g/mol, hence the 1e-3 */
                real Fq_id = 1.0;

                const std::string & tCAS = aGas->data()->cas();

                /* No published Lucas Q exists for helium-3 ( 14762-55-1 );
                 * it deliberately falls through to Fq = 1 rather than
                 * borrowing the helium-4 value */
                real Q = tCAS == "7440-59-7"  ? 1.38 :      // helium
                         tCAS == "1333-74-0"  ? 0.76 :      // hydrogen
                         tCAS == "1333-74-0p" ? 0.76 :      // para-hydrogen, same molecule
                         tCAS == "7782-39-0"  ? 0.52 : 0.0; // deuterium

                if( Q > 0.0 )
                {
                    real sign = T_r < 12.0 ? -1.0 : 1.0;

                    Fq_id = 1.22 * std::pow( Q, 0.15 ) * ( 1.0 + 0.00385 * sign *
                            std::pow( std::pow( T_r - 12.0, 2 ), 1.0e-3 / M ) );
                }

                // inverse of the reduced viscosity, Eq. ( 9-4.15 ), converted to SI
                real inv_zeta = ( 5.0 / 88.0 ) * std::sqrt( M )
                                * std::pow( p_crit, 2.0 / 3.0 )
                                * std::pow( 10.0 / T_crit, 1.0 / 6.0 );

                // Eq. ( 9-4.16 ); 1e-7 converts micropoise to Pa*s
                return Fp_id * Fq_id * inv_zeta * ( 0.807 * std::pow( T_r, 0.618 )
                                            - 0.357 * std::exp( -0.449 * T_r )
                                            + 0.34 * std::exp( -4.058 * T_r )
                                            + 0.018 ) * 1e-7;
            }
            else
            {
                return 0.0;
            }
        }

//------------------------------------------------------------------------------
    }
}