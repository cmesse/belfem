//
// Created by Christian Messe on 26.08.19.
//

#include <cmath>
#include "assert.hpp"
#include "fn_GT_idgas_mu.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        real
        idgas_mu( RefGas * aGas, const real aT )
        {
            // get refs
            const real & M      =  aGas->data()->M();
            const real & T_crit =  aGas->data()->T_crit();
            const real & p_crit =  aGas->data()->p_crit();
            const real & Z_crit =  aGas->data()->Z_crit();
            const real & mu     =  aGas->data()->dipole();

            if( aGas->data()->has_crit() )
            {
                // reduced temperature
                real T_r = aT / T_crit;

                // reduced dipole moment ( 89 )
                real mu_r = 52.46e-5 * std::pow( mu / T_crit, 2 ) * p_crit;

                // correction factor ( 90 )
                real Fp_id;

                if ( mu_r < 0.022 )
                {
                    Fp_id = 1.0;
                }
                else if ( mu_r < 0.075 )
                {
                    Fp_id = 1.0 + 30.55 * std::pow( 0.292 - Z_crit, 1.72 );
                }
                else
                {
                    Fp_id = 1.0 + 30.55 * std::pow( 0.292 - Z_crit, 1.72 )
                                  * std::abs( 0.96 + 0.1 * ( T_r - 0.7 ));
                }

                real inv_zeta = ( 5.0 / 88.0 ) * std::sqrt( M )
                                * std::pow( p_crit, 2.0 / 3.0 )
                                * std::pow( 10.0 / T_crit, 1.0 / 6.0 );

                // ( 88 )
                return Fp_id * inv_zeta * ( 0.807 * std::pow( T_r, 0.618 )
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