//
// Created by Claude (gas module audit) on 05.08.26.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "fn_linspace.hpp"
#include "cl_Gas.hpp"

#include "cl_GM_EoS_Hydrogen.hpp"
#include "cl_GM_EoS_Methane.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

// supercritical caloric consistency dh/dT = cp for the hydrogen
// Helmholtz EoS — control for the methane test below: the hydrogen
// phir_tt carries the correct Gaussian-term form
TEST( GASMODELS, Hydrogen_Cp_FD_Supercritical )
{
    Gas tRef( "H2" );
    EoS_Hydrogen tGas( tRef, HelmholtzModel::NormalHydrogen );

    // ~2.3 p_crit, from above T_crit upward
    const real tP = 3e6;
    const uint tN = 21;
    Vector<real> tT = linspace( 40.0, 200.0, tN );

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tGas.cp( tT( k ), tP );
        tExpect( k ) = (   tGas.h( tT( k ) * 1.001, tP )
                         - tGas.h( tT( k ) * 0.999, tP ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-4 );
}

// same check for methane across the near-critical region ( tau ~ 1,
// delta ~ 1 ), where the Gaussian terms of phir_tt dominate
// ( regression gate for the linearized Gaussian block, audit finding M7 )
TEST( GASMODELS, Methane_Cp_FD_NearCritical )
{
    Gas tRef( "CH4" );
    EoS_Methane tGas( tRef );

    // just above p_crit = 4.599 MPa, sweep across the pseudo-critical line
    const real tP = 5e6;
    const uint tN = 21;
    Vector<real> tT = linspace( 185.0, 240.0, tN );

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    /* the step has to be much tighter than in the hydrogen control: cp spikes
     * to 72 kJ/(kg K) on the pseudo-critical line at 193 K, and a central
     * difference over +-0.1 K truncates by 1.2 % right there — enough to miss
     * an r2 of 1 by 1.7e-4 with an exact equation of state. At +-0.02 K the
     * truncation is 3e-7 in r2, while the difference stays four orders above
     * the noise the density solve leaves in h. */
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tGas.cp( tT( k ), tP );
        tExpect( k ) = (   tGas.h( tT( k ) * 1.0001, tP )
                         - tGas.h( tT( k ) * 0.9999, tP ) )
                       / ( 0.0002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-4 );
}
