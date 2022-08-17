//
// Created by Christian Messe on 07.02.22.
//

//
// Created by Christian Messe on 10.01.22.
//

#include <iostream>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Mesh.hpp"
#include "cl_Logger.hpp"
#include "cl_Element_Factory.hpp"
#include "en_IWGs.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_MaxwellJob.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 4 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );


    // close communicator
    return gComm.finalize();
}