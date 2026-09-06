//
// Created by Claude (gas module audit) on 05.08.26.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_linspace.hpp"
#include "fn_r2.hpp"

#define protected public
#define private   public
#include "GT_globals.hpp"
#include "cl_GT_RefGasFactory.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_HeatPoly.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;

// finite-difference consistency of the NASA-9 heat polynomial:
// dH/dT = Cp, dS/dT = Cp/T, dCp/dT and d2Cp/dT2 against central differences
TEST( GASTABLES, HeatPoly_Consistency )
{
    gastables::RefGasFactory tFactory;
    gastables::RefGas * tGas = tFactory.create_refgas( "N2" );

    // first tabulated NASA interval of N2
    HeatPoly * tPoly = tGas->mHeatPolys( 0 );

    const uint tN = 101;
    Vector<real> tT = linspace( tPoly->T_min() + 10.0,
                                tPoly->T_max() - 10.0, tN );

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    // dH/dT == Cp
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tPoly->Cp( tT( k ) );
        tExpect( k ) = (   tPoly->H( tT( k ) * 1.001 )
                         - tPoly->H( tT( k ) * 0.999 ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    // dS/dT == Cp/T
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tPoly->Cp( tT( k ) ) / tT( k );
        tExpect( k ) = (   tPoly->S( tT( k ) * 1.001 )
                         - tPoly->S( tT( k ) * 0.999 ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    // dCp/dT
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tPoly->dCpdT( tT( k ) );
        tExpect( k ) = (   tPoly->Cp( tT( k ) * 1.001 )
                         - tPoly->Cp( tT( k ) * 0.999 ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    // d2Cp/dT2
    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tPoly->d2CpdT2( tT( k ) );
        tExpect( k ) = (   tPoly->dCpdT( tT( k ) * 1.001 )
                         - tPoly->dCpdT( tT( k ) * 0.999 ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    delete tGas;
}

// dSdT() must equal Cp/T ( regression gate for the Rm precedence
// defect in HeatPoly::dSdT, audit finding C1 )
TEST( GASTABLES, HeatPoly_dSdT )
{
    gastables::RefGasFactory tFactory;
    gastables::RefGas * tGas = tFactory.create_refgas( "N2" );

    HeatPoly * tPoly = tGas->mHeatPolys( 0 );

    const uint tN = 101;
    Vector<real> tT = linspace( tPoly->T_min() + 10.0,
                                tPoly->T_max() - 10.0, tN );

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tPoly->dSdT( tT( k ) );
        tExpect( k ) = tPoly->Cp( tT( k ) ) / tT( k );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    delete tGas;
}
