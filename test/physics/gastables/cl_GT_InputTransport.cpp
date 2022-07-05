//
// Created by Christian Messe on 25.08.19.
//

//
// Created by Christian Messe on 2019-08-19.
//

#include <gtest/gtest.h>

#define protected public
#define private   public
#include "typedefs.hpp"
#include "fn_GT_data_path.hpp"
#include "cl_GT_InputTransport.hpp"
#include "cl_GT_RefGas.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;

TEST( GASTABLES, input_transport )
{
    // go two steps up
    string tDataPath = gastables::data_path();

    // open the thermo file
    gastables::InputTransport tTransport( tDataPath + "/trans.inp" );

    // create a new reference gas
    gastables::RefGas * tGas = new gastables::RefGas( "CH2Cl2" );

    tTransport.read_data( tGas );

    // test composition
    real tEpsilon = 1e-9;

    // test number of polynomials
    EXPECT_EQ( int( tGas->mViscosityPolys.size() ), 2 );
    EXPECT_EQ( int( tGas->mConductivityPolys.size() ), 2 );

    // test coefficients
    EXPECT_NEAR( tGas->mViscosityPolys(0)->mCoefficients( 0 ), 0.57185884e0 , tEpsilon );
    EXPECT_NEAR( tGas->mViscosityPolys(0)->mCoefficients( 1 ), -0.34599168e3, tEpsilon );
    EXPECT_NEAR( tGas->mViscosityPolys(0)->mCoefficients( 2 ),  0.32975791e5, tEpsilon );
    EXPECT_NEAR( tGas->mViscosityPolys(0)->mCoefficients( 3 ), 0.21786059e1, tEpsilon );

    EXPECT_NEAR( tGas->mConductivityPolys(0)->mCoefficients( 0 ), 0.25979341E00 , tEpsilon );
    EXPECT_NEAR( tGas->mConductivityPolys(0)->mCoefficients( 1 ), -0.10510041e4, tEpsilon );
    EXPECT_NEAR( tGas->mConductivityPolys(0)->mCoefficients( 2 ),  0.11078850e6, tEpsilon );
    EXPECT_NEAR( tGas->mConductivityPolys(0)->mCoefficients( 3 ), 0.51956543e1, tEpsilon );


    // test temperature range
    EXPECT_NEAR( tGas->mViscosityPolys(0)->mTmin, 300.0, tEpsilon );
    EXPECT_NEAR( tGas->mViscosityPolys(0)->mTmax, 1000.0, tEpsilon );
    EXPECT_NEAR( tGas->mConductivityPolys(0)->mTmin, 300.0, tEpsilon );
    EXPECT_NEAR( tGas->mConductivityPolys(0)->mTmax, 1000.0, tEpsilon );
    delete tGas;

}