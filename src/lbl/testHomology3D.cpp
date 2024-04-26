//
// Created by Gregory Giard on 2014-03-12.
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

using namespace belfem ;
using namespace mesh ;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );


    // print_banner();

    bool suggHomology = false; //Do we suggest the homology from conductors' boundaries

    // load the mesh
    Mesh * tMesh = new Mesh( "mesh.msh" );

    index_t tCount = tMesh->check() ;
    std::cout << "reoriented " << tCount << " elements." << std::endl ;
    //tMesh->scale_mesh( 0.001 );
    tMesh->create_edges();
    tMesh->create_faces();
    tCount = tMesh->check() ;
    std::cout << "reoriented " << tCount << " elements." << std::endl ;

    //Flag everything for the creation of the simplicial complex
    for (mesh::Element *tElement: tMesh->elements()) {
        tElement->flag();
        tElement->flag_nodes();
        tElement->flag_edges();
        tElement->flag_faces();
    }

    std::cout << "Creating the simplicial complex ..." << std::endl;
    SimplicialComplex * tSComplex = new SimplicialComplex(tMesh) ;
    std::cout << "Number of 0-chains: " << tSComplex->number_of_0simplices() <<std::endl;
    std::cout << "Number of 1-chains: " << tSComplex->number_of_1simplices() <<std::endl;
    std::cout << "Number of 2-chains: " << tSComplex->number_of_2simplices() <<std::endl;
    std::cout << "Number of 3-chains: " << tSComplex->number_of_ksimplices(3) <<std::endl;

    /*tSComplex->get_kchain(2,31)->getBoundary()->print();

    tSComplex->get_kchain(3,1)->getBoundary()->print();
    tSComplex->get_kchain(3,2)->getBoundary()->print();
    tSComplex->get_kchain(3,3)->getBoundary()->print();
    tSComplex->get_kchain(3,4)->getBoundary()->print();
    tSComplex->get_kchain(3,5)->getBoundary()->print();
     */
    tMesh->unflag_everything();

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


        // Coreduce the simplicial complex
        Timer tTimer2;
        std::cout << "(Co)reducing the complex ..." << std::endl;
        //tSComplex->coreduce_complexPellikka();
        tSComplex->coreduce_complexCCR();
        std::cout << "Complex (co)reduced in: "<< tTimer2.stop()*1e-3 << " s" << std::endl;

        std::cout << "Number of 0-cochains: " << tSComplex->number_of_kcosimplices(0) <<std::endl;
        std::cout << "Number of 1-cochains: " << tSComplex->number_of_kcosimplices(1) <<std::endl;
        std::cout << "Number of 2-cochains: " << tSComplex->number_of_kcosimplices(2) <<std::endl;
        std::cout << "Number of 3-cochains: " << tSComplex->number_of_kcosimplices(3) <<std::endl;

        Timer tTimer4;
        std::cout << "Computing Cohomology ..." << std::endl;
        tCohomology = new Cohomology(tSComplex, tMesh);
        std::cout << "Cohomology generators computed in: "<< tTimer4.stop() << " ms" << std::endl;
        //tCohomology->create_kGeneratorsField(1, tMeshEdge);

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
            // start timer
            Timer tTimer;
            std::cout << "Reducing the complex ..." << std::endl;
            //tSComplex->reduce_complexPellikka();
            tSComplex->reduce_complexCCR();
            std::cout << "Complex reduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

            std::cout << "Number of 0-chains: " << tSComplex->number_of_0simplices() <<std::endl;
            std::cout << "Number of 1-chains: " << tSComplex->number_of_1simplices() <<std::endl;
            std::cout << "Number of 2-chains: " << tSComplex->number_of_2simplices() <<std::endl;
            std::cout << "Number of 3-chains: " << tSComplex->number_of_ksimplices(3) <<std::endl;

            //Compute the homology and cohomology generators
            Timer tTimer3;
            std::cout << "Computing Homology ..." << std::endl;
            tHomology = new Homology(tSComplex, tMesh);
            std::cout << "Homology generators computed in: "<< tTimer3.stop() << " ms" << std::endl;
            //tHomology->create_kGeneratorsField(1, tMeshEdge);
        }

        //Create homology and cohomology fields in mesh
        tHomology->create_kGeneratorsField(1, tMeshEdge);
        tCohomology->create_kGeneratorsField(1, tMeshEdge);

        // Save homology information
        tMeshEdge->save("Homology.hdf5");
    }

    // save the mesh
    tMesh->save("Conductors.e-s");

    // save the homology
    tMeshEdge->save("Homology.e-s");


    //delete tBTree ;
    delete tMeshEdge;

    delete tHomology;

    delete tCohomology;

    delete tMesh;

    delete tSComplex ;

    return gComm.finalize();
}
