//
// Created by Christian Messe on 18.08.20.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Gas.hpp"

#define private public
#define protected public
#include "cl_GM_EoS_Oxygen.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

TEST( GASMODELS, Oxygen_Vapor )
{
    // crate reference mixture
    Gas tRef( "O2" );

    // create EoS
    EoS_Oxygen tO2( tRef );

//----------------------------------------------------------------------------
//  Vapor Pressure
//----------------------------------------------------------------------------

    // VDI heat atlas
    // data in °C
    Vector< real > tTvap = { -218., -216., -214., -212., -210., -208., -206.,
                             -204., -202., -200., -198., -196., -194., -192.,
                             -190., -188., -186., -184., -182., -180., -178.,
                             -176., -174., -172., -170., -168., -166., -164.,
                             -162., -160., -158., -156., -154., -152., -150.,
                             -148., -146., -144., -142., -140., -138., -136.,
                             -134., -132., -130., -128., -126., -124., -122.,
                             -120. };
    // data in K
    tTvap += 273.15;

    // data in bar
    Vector< real > tPvap = { 0.00187, 0.00337, 0.00582, 0.00967, 0.01552,
                             0.02411, 0.0364, 0.05354, 0.0769, 0.10808,
                             0.14892, 0.20148, 0.26807, 0.35123, 0.45372,
                             0.57851, 0.72876, 0.90782, 1.1192, 1.3666,
                             1.6538, 1.9848, 2.3636, 2.7943, 3.2812,
                             3.8286, 4.4408, 5.1223, 5.8776, 6.7111,
                             7.6276, 8.6316, 9.7278, 10.921, 12.216,
                             13.618, 15.131, 16.761, 18.513, 20.393,
                             22.406, 24.558, 26.856, 29.305, 31.915,
                             34.692, 37.646, 40.789, 44.137, 47.71 };
    // convert to Pascal
    tPvap *= 1e5;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function p_vap( T_vap )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    uint tNumPoints = tTvap.length() ;

    Vector< real > tY( tNumPoints );
    Vector< real > tY0( tNumPoints );

    for( uint k=0; k<tNumPoints; ++k )
    {
        tY0( k ) = std::log( tPvap( k ) );
        tY( k ) = std::log( tO2.p_vap( tTvap( k ) ) );
    }

    EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-5 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function dp_vap/dT( T_vap )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // test derivative
    for( uint k=0; k<tNumPoints; ++k )
    {
        tY0( k ) =  ( tO2.p_vap( tTvap( k ) + 0.01 ) - tO2.p_vap( tTvap( k ) - 0.01 ) )/0.02 ;
        tY( k )  = tO2.dpvap_dT( tTvap( k ), tO2.p_vap( tTvap( k ) ), tO2.pi_vap( tTvap( k ) ) );
    }
    EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-6 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function T_vap(p_vap) )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // test derivative
    for( uint k=0; k<tNumPoints; ++k )
    {
        tY( k )  =  tO2.T_vap( tO2.p_vap( tTvap( k ) ) );
    }
    EXPECT_NEAR( r2( tY, tTvap ), 1.0, 1e-6 );

}