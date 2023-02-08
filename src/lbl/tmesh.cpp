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
#include "cl_MaxwellJob.hpp" // delete me

#include "cl_Pipette.hpp"
#include "fn_rcond.hpp"
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // initialize job
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // read the input file
    InputFile tInputFile( "input.txt" );


    // create a new factory
    MaxwellFactory * tFactory = new MaxwellFactory( tInputFile );

    // get the mesh
    Mesh * tInMesh = tFactory->mesh();

    // get the kernel
    Kernel * tKernel = tFactory->kernel() ;


    // get the magnetic field
    DofManager * tMagfield = tFactory->magnetic_field() ;

    // get the formulation
    IWG_Maxwell * tFormulation = reinterpret_cast< IWG_Maxwell * >( tMagfield->iwg() );


    // create the mesh extractor
    ThermalMeshExtractor tExtractor( tInMesh, tFormulation->ghost_blocks() );

    // create the output mesh
    Mesh * tOutMesh = tExtractor.create_mesh() ;

    // save the mesh
    tOutMesh->save( "outmesh.vtk");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // tidy up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    delete tKernel ;
    delete tInMesh ;
    delete tOutMesh ;

    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
