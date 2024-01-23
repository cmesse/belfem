//
// Created by grgia@ge.polymtl.ca on 06/12/23.
//

#include <iostream>
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Homology.hpp"
#include "cl_Cohomology.hpp"
#include "banner.hpp"
#include "cl_Node.hpp"
#include "cl_SideSet.hpp"

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
#include "fn_Smith.hpp"
#include "fn_sum.hpp"
#include "fn_max.hpp"

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
using namespace mesh ;

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

    // Create the simplicial complex
    std::cout << "Creating the simplicial complex ..." << std::endl;
    SimplicialComplex * tSComplex = tFactory->simplicial_complex();
    //SimplicialComplex * tSComplex = new SimplicialComplex(tMesh);

    std::cout << "Number of 0-chains: " << tSComplex->number_of_0simplices() <<std::endl;
    std::cout << "Number of 1-chains: " << tSComplex->number_of_1simplices() <<std::endl;
    std::cout << "Number of 2-chains: " << tSComplex->number_of_2simplices() <<std::endl;


    //------------------------------------------------------------------------------

    // Reduce the simplicial complex
    // start timer
    Timer tTimer;
    std::cout << "Reducing the complex ..." << std::endl;
    tSComplex->reduce_complexCCR();
    std::cout << "Complex reduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

    std::cout << "Number of 0-chains: " << tSComplex->number_of_0simplices() <<std::endl;
    std::cout << "Number of 1-chains: " << tSComplex->number_of_1simplices() <<std::endl;
    std::cout << "Number of 2-chains: " << tSComplex->number_of_2simplices() <<std::endl;

    // Coreduce the simplicial complex
    Timer tTimer2;
    std::cout << "(Co)reducing the complex ..." << std::endl;
    tSComplex->coreduce_complexCCR();
    std::cout << "Complex (co)reduced in: "<< tTimer2.stop()*1e-3 << " s" << std::endl;

    std::cout << "Number of 0-cochains: " << tSComplex->number_of_kcosimplices(0) <<std::endl;
    std::cout << "Number of 1-cochains: " << tSComplex->number_of_kcosimplices(1) <<std::endl;
    std::cout << "Number of 2-cochains: " << tSComplex->number_of_kcosimplices(2) <<std::endl;

    // get the name for the outmesh
    const string tOutFile = tFactory->outfile();

    //Compute the homology and cohomology generators
    Timer tTimer3;
    std::cout << "Computing Homology ..." << std::endl;
    Homology * tHomology = new Homology(tSComplex, tMesh);
    std::cout << "Homology generators computed in: "<< tTimer3.stop() << " ms" << std::endl;
    tHomology->create_kGeneratorsField(1);

    Timer tTimer4;
    std::cout << "Computing Cohomology ..." << std::endl;
    Cohomology * tCohomology = new Cohomology(tSComplex, tMesh);
    std::cout << "Cohomology generators computed in: "<< tTimer4.stop() << " ms" << std::endl;
    tCohomology->create_kGeneratorsField(1);

    //Define chain on conductor edge for homology
    std::cout << "Suggesting an homology from the conductor boundary ..." << std::endl;
    tHomology->suggest_Homology();
    tHomology->create_kGeneratorsField(1);

    std::cout << "Updating Cohomology from suggested homology ..." << std::endl;
    tCohomology->updatekGeneratorsFromHomology(tHomology->get_Generators()(1),1);
    tCohomology->create_kGeneratorsField(1);

    // save the mesh
    tMesh->save(tOutFile);

    delete tFactory;

    delete tHomology;


    delete tCohomology;

    delete tMesh ;
    return 0;

}
