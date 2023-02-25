//
// Created by christian on 2/16/23.
//
//
// Created by christian on 2/7/23.
//
#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Timer.hpp"

#include "cl_InputFile.hpp"
#include "cl_MaxwellFactory.hpp"
#include "cl_Profiler.hpp"
#include "fn_FEM_compute_normb.hpp"
#include "fn_sum.hpp"
#include "cl_Maxwell_ThermalMeshExtractor.hpp"

//#include "fn_FEM_compute_element_current_thinshell.hpp"
#include "cl_IWG_Maxwell_2D.hpp"

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

using namespace belfem ;
using namespace fem ;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    if ( comm_rank() == 0 )
    {
        print_banner( "- input file test program -" );
    }


    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
