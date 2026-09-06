//
// Created by Claude (gas module audit) on 05.08.26.
//

#include <gtest/gtest.h>
#include <cmath>

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
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;

// POLY and SPLINE mode must agree on Cp, H, S, mu and lambda
TEST( GASTABLES, RefGas_SplinePolyCrossCheck )
{
    gastables::RefGasFactory tFactory;
    gastables::RefGas * tGas = tFactory.create_refgas( "N2" );

    const uint tN = 37;
    Vector<real> tT = linspace( 200.0, 2000.0, tN );

    Vector<real> tPoly( tN );
    Vector<real> tSpline( tN );

    // Cp
    tGas->set_mode( RefGasMode::POLY );
    for ( uint k = 0; k < tN; ++k )
    {
        tPoly( k ) = tGas->Cp( tT( k ) );
    }
    tGas->set_mode( RefGasMode::SPLINE );
    for ( uint k = 0; k < tN; ++k )
    {
        tSpline( k ) = tGas->Cp( tT( k ) );
    }
    EXPECT_NEAR( r2( tSpline, tPoly ), 1.0, 1e-5 );

    // H
    tGas->set_mode( RefGasMode::POLY );
    for ( uint k = 0; k < tN; ++k )
    {
        tPoly( k ) = tGas->H( tT( k ) );
    }
    tGas->set_mode( RefGasMode::SPLINE );
    for ( uint k = 0; k < tN; ++k )
    {
        tSpline( k ) = tGas->H( tT( k ) );
    }
    EXPECT_NEAR( r2( tSpline, tPoly ), 1.0, 1e-5 );

    // S
    tGas->set_mode( RefGasMode::POLY );
    for ( uint k = 0; k < tN; ++k )
    {
        tPoly( k ) = tGas->S( tT( k ) );
    }
    tGas->set_mode( RefGasMode::SPLINE );
    for ( uint k = 0; k < tN; ++k )
    {
        tSpline( k ) = tGas->S( tT( k ) );
    }
    EXPECT_NEAR( r2( tSpline, tPoly ), 1.0, 1e-5 );

    // mu
    tGas->set_mode( RefGasMode::POLY );
    for ( uint k = 0; k < tN; ++k )
    {
        tPoly( k ) = tGas->mu( tT( k ) );
    }
    tGas->set_mode( RefGasMode::SPLINE );
    for ( uint k = 0; k < tN; ++k )
    {
        tSpline( k ) = tGas->mu( tT( k ) );
    }
    EXPECT_NEAR( r2( tSpline, tPoly ), 1.0, 1e-5 );

    // lambda
    tGas->set_mode( RefGasMode::POLY );
    for ( uint k = 0; k < tN; ++k )
    {
        tPoly( k ) = tGas->lambda( tT( k ) );
    }
    tGas->set_mode( RefGasMode::SPLINE );
    for ( uint k = 0; k < tN; ++k )
    {
        tSpline( k ) = tGas->lambda( tT( k ) );
    }
    EXPECT_NEAR( r2( tSpline, tPoly ), 1.0, 1e-5 );

    delete tGas;
}

// h_ref must be the mass-specific value H_ref / M ( regression gate
// for the unit inversion in h_ref, audit finding C2 )
TEST( GASTABLES, RefGas_href_units )
{
    gastables::RefGasFactory tFactory;
    gastables::RefGas * tGas = tFactory.create_refgas( "N2" );

    real tExpect = tGas->H_ref() / tGas->mData.M();

    EXPECT_NEAR( tGas->h_ref(), tExpect,
                 std::abs( tExpect ) * 1e-9 + 1e-12 );

    delete tGas;
}

// RefGas::dSdT in POLY mode against a central difference of S
// ( regression gate for audit finding C1 at the RefGas level )
TEST( GASTABLES, RefGas_dSdT )
{
    gastables::RefGasFactory tFactory;
    gastables::RefGas * tGas = tFactory.create_refgas( "N2" );

    tGas->set_mode( RefGasMode::POLY );

    const uint tN = 37;
    Vector<real> tT = linspace( 250.0, 2000.0, tN );

    Vector<real> tValues( tN );
    Vector<real> tExpect( tN );

    for ( uint k = 0; k < tN; ++k )
    {
        tValues( k ) = tGas->dSdT( tT( k ) );
        tExpect( k ) = (   tGas->S( tT( k ) * 1.001 )
                         - tGas->S( tT( k ) * 0.999 ) )
                       / ( 0.002 * tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 1e-6 );

    delete tGas;
}
