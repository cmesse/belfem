//
// Created by christian on 3/1/23.
//
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "FEM_geometrytools.hpp"
using namespace belfem ;
using namespace fem ;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    print_banner();

    std::cout << "eps " << BELFEM_EPS << std::endl ;

    // close communicator
    return gComm.finalize();
}