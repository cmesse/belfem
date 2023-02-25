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
    Mesh * tMesh = tFactory->magnetic_mesh();

    tMesh->create_field( "elementT").fill( 10.0 );

    // get the kernel
    Kernel * tKernel = tFactory->magnetic_kernel() ;


    // get the magnetic field
    DofManager * tMagfield = tFactory->magnetic_field() ;

    // get the formulation
    IWG_Maxwell * tFormulation = reinterpret_cast< IWG_Maxwell * >( tMagfield->iwg() );


    // create the mesh extractor
    ThermalMeshExtractor tExtractor( tMesh );

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
    tThermalKernel->mesh()->create_field( "T").fill( 10.0 );

    comm_barrier() ;
    tSynch.initialize_tables() ;

    tMagfield->initialize() ;
    IWG_Maxwell_Thermal2D * tFourier = new IWG_Maxwell_Thermal2D( ) ;
    tFourier->select_block( 1 );
    DofManager * tThermal = tThermalKernel->create_field( tFourier );

    tFourier->set_material_table( tFactory->tape_materials() );
    tFourier->compute_geometry_data( tKernel->mesh(),
                                     tThermalKernel->mesh(),
                                     tFactory->tape_thicknesses() );



    tFourier->initialize() ;

    tThermal->set_solver( SolverType::UMFPACK );
    tThermal->initialize() ;


    tSynch.magnetic_to_thermal_b_and_ej()

    tThermal->compute_jacobian_and_rhs();

    tThermal->solve() ;

    tSynch.thermal_to_magnetic_T() ;


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


    delete tThermalMesh ;

    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
