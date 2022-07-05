//
// Created by Christian Messe on 25.10.19.
//


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
    gComm = Communicator( &argc, &argv );

    fem::InterpolationFunctionFactory tFactory;

    auto tFunction = tFactory.create_lagrange_function( ElementType::LINE2 );

    delete tFunction;
    return gComm.finalize();
}