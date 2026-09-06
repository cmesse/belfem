//
// Created by Christian Messe on 2019-08-19.
//

#include <gtest/gtest.h>

#define protected public
#define private   public
#include "typedefs.hpp"
#include "fn_GT_data_path.hpp"
#include "cl_GT_InputThermo.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_GasData.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;

TEST( GASTABLES, input_thermo )
{
    // go two steps up
    string tDataPath = gastables::data_path();

    // open the thermo file
    gastables::InputThermo tThermo( tDataPath + "/thermo.inp" );

    // create a new reference gas
    gastables::RefGas * tGas = new gastables::RefGas( "CF2ClBr" );

    tThermo.read_data( tGas );


    // test composition

    real tEpsilon = 1e-9;

    EXPECT_NEAR( tGas->component_multiplicity("C"),  1.0, tEpsilon );
    EXPECT_NEAR( tGas->component_multiplicity("F"),  2.0, tEpsilon );
    EXPECT_NEAR( tGas->component_multiplicity("Cl"),  1.0, tEpsilon );
    EXPECT_NEAR( tGas->component_multiplicity("Br"),  1.0, tEpsilon );
    EXPECT_NEAR( tGas->component_multiplicity("Xe"),  0.0, tEpsilon );

    // test molar mass
    EXPECT_NEAR( tGas->M(), 0.1653645064, tEpsilon );

    // test formation enthalpy
    EXPECT_NEAR( tGas->reference_formation_enthalpy(), -435.0e3, tEpsilon );

    // test number of polynomials
    EXPECT_EQ( int( tGas->mHeatPolys.size() ), 2 );

    // test coefficients
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 0 ), 2.6966993650e4 , tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 1 ), -4.67883227e2, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 2 ), 5.24480942e0, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 3 ), 2.530121349e-2, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 4 ), -3.50844352e-5, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 5 ), 2.354353958e-8, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mCoefficients( 6 ), -6.246161130e-12, tEpsilon );

    // test temperature range
    EXPECT_NEAR( tGas->mHeatPolys(0)->mTmin, 200.0, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mTmax, 1000.0, tEpsilon );

    // test reference parameters
    EXPECT_NEAR( tGas->mHeatPolys(0)->mEnthalpyConstant, -5.198382730e4, tEpsilon );
    EXPECT_NEAR( tGas->mHeatPolys(0)->mEntropyConstant, 8.532887803e-1, tEpsilon );

    delete tGas;
}