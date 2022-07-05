//
// Created by Christian Messe on 14.10.19.
//

#ifndef BELFEM_CL_GM_GAS_PRANDTLMEYER_HPP
#define BELFEM_CL_GM_GAS_PRANDTLMEYER_HPP

#include <gtest/gtest.h>


#include "typedefs.hpp"
#include "constants.hpp"

#include "cl_Communicator.hpp"
#include "cl_Gas.hpp"
using namespace belfem;
using namespace belfem::gastables;

TEST( GASMODELS, PrandtlMeyer )
{
    Gas tAir;

    real tT1 = 216.65;
    real tP1 = 19395.0;
    real tMa1 = 2.0;
    real tAlpha = 20 * constant::deg;

    real tU1 = tAir.c( tT1, tP1 ) * tMa1;

    real tNu1 = tAir.prandtl_meyer_angle( tT1, tP1, tU1 ) / constant::deg;

    EXPECT_NEAR( tNu1, 26.4, 0.1 );

    real tT2;
    real tP2;
    real tU2;

    tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha, tT2, tP2, tU2 );

    real tMa2 = tU2 / tAir.c( tT2, tP2 );

    EXPECT_NEAR( tMa2, 2.83, 0.01 );

    EXPECT_NEAR( tP2/tP1, 0.271, 0.01 );
    EXPECT_NEAR( tT2/tT1, 0.689, 0.01 );
}

#endif //BELFEM_CL_GM_GAS_PRANDTLMEYER_HPP
