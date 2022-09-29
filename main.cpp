
#include <iostream>

#include <cholmod.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "fn_rcond.hpp"
using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    print_banner();

    Matrix <uint > tA = { { 0, 1, 1 }, { 1, 1, 1 }, { 1, 2, 3 } } ;

    tA.print("A");

    // close communicator
    return gComm.finalize();
}