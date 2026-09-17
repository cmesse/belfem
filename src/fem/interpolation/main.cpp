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

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "Mesh_Enums.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"


using namespace belfem;

Communicator gComm;

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    fem::InterpolationFunctionFactory tFactory;

    auto tFunction = tFactory.create_lagrange_function( ElementType::LINE2 );

    delete tFunction;
    return gComm.finalize();
}