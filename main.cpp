
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

    // print_banner();

    Matrix <real> tA( 5, 5, 0.0 );

    tA( 0, 0 ) =  0.0;
    tA( 0, 1 ) = -1.0;
    tA( 0, 3 ) = -3.0;
    tA( 1, 0 ) = -2.0;
    tA( 1, 1 ) =  5.0;
    tA( 2, 2 ) =  4.0;
    tA( 2, 3 ) =  6.0;
    tA( 2, 4 ) =  4.0;
    tA( 3, 0 ) = -4.0;
    tA( 3, 2 ) =  2.0;
    tA( 3, 3 ) =  7.0;
    tA( 4, 1 ) =  8.0;
    tA( 4, 4 ) = -5.0;

    // require suitesparse and openmp!
    tA.print("A");
    SpMatrix aMatrix( tA, SpMatrixType::CSR );

    real tCond = rcond( aMatrix );

    std::cout << "cond " << tCond << std::endl ;

    // close communicator
    return gComm.finalize();
}