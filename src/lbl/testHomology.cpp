//
// Created by grgia@ge.polymtl.ca on 06/12/23.
//

#include <iostream>
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_BeltedTree.hpp"
#include "cl_Homology.hpp"
#include "cl_Cohomology.hpp"
#include "banner.hpp"
#include "cl_Node.hpp"
#include "cl_SideSet.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Mesh_HDF5Reader.hpp"

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

    bool suggHomology = true; //Do we suggest the homology from conductors' boundaries

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

    // get the name for the outmesh
    const string tOutFile = tFactory->outfile();

    // --------------------------------------------------------------
    // Creating the edge mesh for homology saving and  vizualization
    // --------------------------------------------------------------
    Mesh * tMeshEdge = new Mesh( 2 , 0, false );
    Cell< Node * > & tNodes = tMeshEdge->nodes();
    Map< id_t, Node * > tNodeMap ;
    id_t tNumNodes = tMesh->number_of_nodes();
    id_t tNumEdges = tMesh->number_of_edges();

    // allocate node memory
    tNodes.set_size( tNumNodes, nullptr );
    for ( uint k=0; k<tNumNodes; ++k)
    {
        // grab data
        id_t id = tMesh->nodes()(k)->id();
        real x = tMesh->nodes()(k)->x();
        real y = tMesh->nodes()(k)->y();
        real z = tMesh->nodes()(k)->z();

        // create a new node
        Node * tNode =  new Node( id, x, y, z );

        // add node to container
        tNodes( k ) = tNode ;

        // add node to map
        tNodeMap[ id ] = tNode ;

    }

    mesh::Block * tBlock = new mesh::Block( 1, tNumEdges );
    tBlock->label() = "Edges";
    // create element factory
    ElementFactory tElementFactory ;

    // grab element container
    Cell< mesh::Element * > & tEdges = tBlock->elements();

    //Cell< mesh::Edge * > & tEdges = tMesh->edges();
    tEdges.set_size( tNumEdges, nullptr );
    for ( uint e=0; e<tNumEdges; ++e)
    {
        id_t edge_id = tMesh->edges()(e)->id();
        id_t node_id_1 = tMesh->edges()(e)->node(0)->id();
        id_t node_id_2 = tMesh->edges()(e)->node(1)->id();

        mesh::Element * tEdge = tElementFactory.create_element( ElementType::LINE2, edge_id );
        tEdge->insert_node( tNodeMap(node_id_1), 0  );

        tEdge->insert_node( tNodeMap( node_id_2 ), 1 );
        tEdges( e ) = tEdge ;

    }

    // add the block to the mesh
    tMeshEdge->blocks().push( tBlock );
    tMeshEdge->finalize();

    // ---------------------------------------------------------
    // Computing Homology
    // ---------------------------------------------------------

    Homology * tHomology;
    Cohomology * tCohomology;

    // Check if we already have a homology file
    if (file_exists("Homology.hdf5"))
    {
        std::cout << "Homology already exists, loading data" << std::endl;
        Mesh * tMeshEdgeLoaded = new Mesh( 2 , 0, false );
        HDF5Reader * tReader = new HDF5Reader( "Homology.hdf5", tMeshEdgeLoaded );

        tHomology = new Homology(tMesh, tMeshEdgeLoaded);
        tCohomology = new Cohomology(tMesh, tMeshEdgeLoaded);

        //Create homology and cohomology fields in mesh
        tHomology->create_kGeneratorsField(1, tMeshEdge);
        tCohomology->create_kGeneratorsField(1, tMeshEdge);

        delete tReader;
    }
    else
    {
        // Create the simplicial complex
        std::cout << "Creating the simplicial complex ..." << std::endl;
        SimplicialComplex * tSComplex = tFactory->simplicial_complex();

        std::cout << "Number of 0-chains: " << tSComplex->number_of_kcosimplices(0) <<std::endl;
        std::cout << "Number of 1-chains: " << tSComplex->number_of_kcosimplices(1) <<std::endl;
        std::cout << "Number of 2-chains: " << tSComplex->number_of_kcosimplices(2) <<std::endl;
        std::cout << "Number of 3-chains: " << tSComplex->number_of_kcosimplices(3) <<std::endl;

        // Coreduce the simplicial complex
        //tSComplex->coreduce_complexPellikka();
        tSComplex->coreduce_complexCCR();

        std::cout << "Computing Cohomology ..." << std::endl;
        tCohomology = new Cohomology(tSComplex, tMesh);
        std::cout << tCohomology->get_Generators()(1).size() << " 1-cohomology generators computed" << std::endl;
        std::cout << std::endl;

        // Do we want to suggest a homology?
        if (suggHomology)
        {
            //Define chain on conductor edge for homology
            std::cout << "Suggesting an homology from the conductor boundary ..." << std::endl;
            tHomology = new Homology( tMesh);

            std::cout << "Updating Cohomology from suggested homology ..." << std::endl;
            tCohomology->updatekGeneratorsFromHomology(tHomology->get_Generators()(1),1);
        }
        else // Otherwise, compute the homology as usual
        {
            // Reduce the simplicial complex
            //tSComplex->reduce_complexPellikka();
            tSComplex->reduce_complexCCR();

            //Compute the homology and cohomology generators
            std::cout << "Computing Homology ..." << std::endl;
            tHomology = new Homology(tSComplex, tMesh);
            std::cout << tHomology->get_Generators()(1).size() << " 1-homology generators computed" << std::endl;
            std::cout << std::endl;
        }

        //Create homology and cohomology fields in mesh
        tHomology->create_kGeneratorsField(1, tMeshEdge);
        tCohomology->create_kGeneratorsField(1, tMeshEdge);

        // Save homology information
        tMeshEdge->save("Homology.hdf5");
    }

    // save the mesh
    tMesh->save(tOutFile);


    // Create the Belted Tree (extremely slow for now...)
    std::cout << "Creating new simplicial complex for belted tree" << std::endl ;
    tFactory->flag_nc_simplices();
    SimplicialComplex * tSComplex2 = new SimplicialComplex(tMesh) ;
    tFactory->unflag_nc_simplices();

    Timer tTimer5;
    std::cout << "Creating the belted tree ..." << std::endl ;
    BeltedTree * tBTree = new BeltedTree(tMesh,  tSComplex2, tHomology->get_Generators()(1)) ;
    std::cout << "Belted tree computed in: "<< tTimer5.stop()*1e-3 << " s" << std::endl;
    tBTree->create_TreeField("BeltedTree");

    Timer tTimer6;
    std::cout << "Computing cohomology from belted tree ..." << std::endl ;
    tBTree->compute_cohomology();
    std::cout << "Cohomology computed in: "<< tTimer6.stop() << " ms" << std::endl;
    tBTree->create_cohomologyField( tMeshEdge );

    // save the homology
    tMeshEdge->save("Homology.e-s");

    //delete tBTree ;

    delete tMeshEdge;

    delete tFactory;

    delete tHomology;

    delete tCohomology;

    delete tMesh;


    return 0;

}
