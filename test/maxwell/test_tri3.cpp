//
// Created by Christian Messe on 18.02.22.
//

#include <gtest/gtest.h>
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "en_IWGs.hpp"
#include "cl_MaxwellJob.hpp"

using namespace belfem ;
using namespace fem ;

TEST( MAXWELL, tri3 )
{
    if( comm_rank() == 0 )
    {
        // create the testcase
        fem::MaxwellJob tJob;

        Vector< real > tResults ;

        tJob.initialize_test( IwgType::MAXWELL_HPHI_TRI3 );
        tJob.run_test( "/tmp/test_maxwell.hdf5", tResults ) ;

        for( real tX : tResults )
        {
            EXPECT_NEAR( tX, 1.0, 1e-12 );
        }
    }
}