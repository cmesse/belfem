//
// Created by Claude (gas module audit) on 05.08.26.
//

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "fn_linspace.hpp"
#include "GT_globals.hpp"
#include "cl_Gas.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

// alpha = p * beta * kappa must hold regardless of query order
TEST( GASMODELS, Cubic_AlphaBetaKappa_Identity )
{
    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    Gas tGas( tSpecies, tMolarFractions, GasModel::SRK );

    const uint tN = 21;
    Vector<real> tT = linspace( 200.0, 800.0, tN );
    const real tP = 50e5;

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tGas.alpha( tT( k ), tP );
        tExpect( k ) = tP * tGas.beta(  tT( k ), tP )
                          * tGas.kappa( tT( k ), tP );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-9 );
}

// calling kappa() before alpha() at the same state must not change alpha
// ( regression gate for the KAPPA statebit mixup, audit finding M3 )
TEST( GASMODELS, Cubic_AlphaAfterKappa )
{
    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    const real tT = 300.0;
    const real tP = 50e5;

    // reference: fresh object, alpha queried first
    Gas tRef( tSpecies, tMolarFractions, GasModel::SRK );
    real tExpect = tRef.alpha( tT, tP );

    // same state, but kappa queried first
    Gas tGas( tSpecies, tMolarFractions, GasModel::SRK );
    tGas.kappa( tT, tP );
    real tValue = tGas.alpha( tT, tP );

    EXPECT_NEAR( tValue, tExpect, std::abs( tExpect ) * 1e-9 );
}

// after remix, properties at an already-queried ( T, p ) must match a
// freshly built gas of the new composition ( regression gate for the
// missing cache invalidation, audit findings M1 + M2 )
TEST( GASMODELS, Gas_Remix_Statevals )
{
    Cell<string> tSpecies = { "N2", "O2" };

    Vector<real> tMixA = { 0.79, 0.21 };
    Vector<real> tMixB = { 0.50, 0.50 };

    const real tT = 300.0;
    const real tP = 50e5;

    // reference gas built directly with composition B
    Gas tRef( tSpecies, tMixB, GasModel::SRK );
    real tCpExpect = tRef.cp( tT, tP );
    real tVExpect  = tRef.v(  tT, tP );

    // gas built with composition A, queried, then remixed to B
    Gas tGas( tSpecies, tMixA, GasModel::SRK );
    tGas.cp( tT, tP );
    tGas.v(  tT, tP );

    tGas.remix( tMixB );

    EXPECT_NEAR( tGas.cp( tT, tP ), tCpExpect,
                 std::abs( tCpExpect ) * 1e-9 );
    EXPECT_NEAR( tGas.v( tT, tP ), tVExpect,
                 std::abs( tVExpect ) * 1e-9 );
}
