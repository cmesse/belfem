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

#include "cl_GM_EoS_Methane.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

// caloric consistency of the real-gas assembly: dh/dT = cp, ds/dT = cp/T
TEST( GASMODELS, Gas_Caloric_FD )
{
    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    Gas tGas( tSpecies, tMolarFractions, GasModel::SRK );

    const uint tN = 31;
    Vector<real> tT = linspace( 250.0, 800.0, tN );
    const real tP = 1e5;

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    // dh/dT == cp
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tGas.cp( tT( k ), tP );
        tExpect( k ) = (   tGas.h( tT( k ) * 1.001, tP )
                         - tGas.h( tT( k ) * 0.999, tP ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    // ds/dT == cp/T
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tGas.cp( tT( k ), tP ) / tT( k );
        tExpect( k ) = (   tGas.s( tT( k ) * 1.001, tP )
                         - tGas.s( tT( k ) * 0.999, tP ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );
}

// in the low-pressure limit the real-gas sound speed must approach the
// ideal value sqrt( gamma R T ) ( regression gate for the missing cp
// factor in realgas_c, audit finding M5 )
TEST( GASMODELS, Gas_SoundSpeed_IdgasLimit )
{
    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    Gas tGas( tSpecies, tMolarFractions, GasModel::SRK );

    const real tT = 300.0;
    const real tP = 1e3;    // low pressure: departures negligible

    real tC     = tGas.c( tT, tP );
    real tGamma = tGas.gamma( tT, tP );
    real tR     = tGas.R( tT, tP );

    real tExpect = std::sqrt( tGamma * tR * tT );

    EXPECT_NEAR( tC, tExpect, tExpect * 1e-2 );
}

// at the reference pressure the cubic models are documented to reduce to
// the ideal gas model ( comment block above Gas::realgas_cp ); serves as
// the numeric probe for the suspected mSref double-count ( cl_Gas.cpp,
// audit finding M8 ). Enable after adjudication.
TEST( GASMODELS, DISABLED_Gas_Entropy_RefState )
{
    Cell<string> tSpecies = { "N2", "O2" };
    Vector<real> tMix = { 0.79, 0.21 };

    Gas tIdgas( tSpecies, tMix, GasModel::IDGAS );
    Gas tSRK(   tSpecies, tMix, GasModel::SRK );

    const real tT = gastables::gTref;
    const real tP = gastables::gPref;

    real tExpect = tIdgas.s( tT, tP );
    real tValue  = tSRK.s( tT, tP );

    EXPECT_NEAR( tValue, tExpect, std::abs( tExpect ) * 1e-4 );
}

// u = h - p v must hold for the Helmholtz reference offsets
// ( regression gate for the double pv subtraction in
// Helmholtz::set_reference_point, audit finding M9 )
TEST( GASMODELS, Helmholtz_u_Consistency )
{
    Gas tRef( "CH4" );
    EoS_Methane tGas( tRef );

    const real tT = 300.0;
    const real tP = 1e5;

    real tExpect = tGas.h( tT, tP ) - tP * tGas.v( tT, tP );
    real tValue  = tGas.u( tT, tP );

    EXPECT_NEAR( tValue, tExpect, std::abs( tExpect ) * 1e-6 + 1.0 );
}
