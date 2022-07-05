//
// Created by Christian Messe on 26.08.19.
//

#include <gtest/gtest.h>

#define protected public
#define private   public
#include "typedefs.hpp"
#include "fn_GT_data_path.hpp"
#include "cl_GT_GasData.hpp"
#include "cl_GT_InputData.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;

TEST( GASTABLES, input_data )
{
    // go two steps up
    string tDataPath = gastables::data_path();

    // open the data file
    InputData tInputData( tDataPath + "/gasdata.inp" );

    // create a new data pointer
    auto tData = new GasData();

    // set the label of the gas
    tData->set_label("BrF3");

    // make sure that entry exists
    EXPECT_TRUE( tInputData.entry_exists( tData ) );

    // read the data from the input file
    tInputData.read_data( tData );

    real tEpsilon = 1e-9;

    EXPECT_NEAR( tData->M(), 0.136899, tEpsilon );
    EXPECT_NEAR( tData->T_crit(), 600.0, tEpsilon );
    EXPECT_NEAR( tData->p_crit(), 49.9e5, tEpsilon );
    EXPECT_NEAR( tData->Z_crit(), 0.115, tEpsilon );

    EXPECT_NEAR( tData->acentric(),  0.413, tEpsilon );
    EXPECT_NEAR( tData->dipole(), 1.1, tEpsilon );

    // delete the data pointer
    delete tData;
}