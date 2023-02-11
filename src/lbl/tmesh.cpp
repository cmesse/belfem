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


#include "cl_MaxwellMeshSynch.hpp"
#include "cl_IWG_Maxwell_Thermal2D.hpp"

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
    Mesh * tMesh = tFactory->mesh();

    // get the kernel
    Kernel * tKernel = tFactory->kernel() ;


    // get the magnetic field
    DofManager * tMagfield = tFactory->magnetic_field() ;

    // get the formulation
    IWG_Maxwell * tFormulation = reinterpret_cast< IWG_Maxwell * >( tMagfield->iwg() );

    // create the mesh extractor
    ThermalMeshExtractor tExtractor( tMesh, tFormulation->ghost_blocks() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Initialize temperature problem
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    // create the output mesh
    Mesh * tThermalMesh = tExtractor.create_mesh() ;

    comm_barrier() ;

    MaxwellMeshSynch tSynch( tMesh, tThermalMesh );

    KernelParameters * tThermalParams = new KernelParameters( tThermalMesh );
    tThermalParams->set_auto_partition( false );

    Kernel * tThermalKernel = new Kernel( tThermalParams );
    tThermalKernel->claim_parameter_ownership();
    comm_barrier() ;
    tSynch.initialize_tables() ;

    tMagfield->initialize() ;
    IWG * tFourier = new IWG_Maxwell_Thermal2D ;
    tFourier->select_block( 1 );
    DofManager * tThermal = tThermalKernel->create_field( tFourier );


    tFourier->initialize() ;

    // tSynch.maxwell_to_thermal() ;


    // save the mesh
    if( comm_rank() == 0 )
    {
        tThermalMesh->save( "outmesh.vtk");
    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // tidy up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    delete tThermalKernel ;
    // delete tFourier ;

    delete tKernel ;
    delete tMesh ;

    delete tThermalKernel ;
    delete tThermalMesh ;

    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
