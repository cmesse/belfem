/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "banner.hpp"
#include "cl_Communicator.hpp"
#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "cl_Matrix.hpp"
#include "cl_Timer.hpp"
#include "cl_Vector.hpp"


using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );


    print_banner( );

    return gComm.finalize();
}
