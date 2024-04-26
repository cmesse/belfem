//
// Created by Christian Messe on 2019-07-25.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_OrientationChecker.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );


    // print_banner();

    // load the mesh
    Mesh * tMesh = new Mesh( "tetsindie.msh" );



    // check orientations (optional)
    index_t tCount = tMesh->check() ;
    std::cout << "reoriented " << tCount << " elements." << std::endl ;


    // convert node coordinates from mm to m
    tMesh->scale_mesh( 0.001 );
    tMesh->create_edges( false );
    tMesh->create_faces( false );

    //index_t tCount = tMesh->check() ;
    //std::cout << "reoriented " << tCount << " elements." << std::endl ;

    std::cout << std::endl ;

    for( mesh::Element * tElement : tMesh->elements() )
    {
        std::cout << tElement->id() << " :" ;
        for( uint k=0; k<tElement->number_of_nodes(); ++k )
        {
            std::cout << " " << tElement->node( k )->id() ;
        }
        std::cout << std::endl ;
    }

    std::cout << std::endl ;

    mesh::Edge * tEdge = tMesh->edges()(0);
    std::cout << tEdge->number_of_faces() << std::endl ;

    for (mesh::Edge * tEdge : tMesh->edges())
    {
        std::cout << "Edge ID: " <<  tEdge->id() << std::endl ;
        std::cout << "Node 1: " <<  tEdge->node(0)->id() << ", Node 2: " << tEdge->node(1)->id() << std::endl ;
        std::cout << "Number of faces: " << tEdge->number_of_faces() << std::endl ;

        for (uint i = 0 ; i < tEdge->number_of_faces() ; i++)
        {
            mesh::Face * tFace = tEdge->face(i);
            std::cout << "Face " << i << " ID: "<<tFace->id() << std::endl;
            for(uint j = 0; j < 3; j++)
            {
                if (tFace->edge(j)->id() == tEdge->id())
                {
                    std::cout << "Edge orientation on face : "<< tFace->edge_orientation(j) << std::endl;
                }
            }
        }
        std::cout << std::endl ;
    }

    for (mesh::Face * tFace : tMesh->faces())
    {
        std::cout << "Face ID: " <<  tFace->id() << std::endl ;
        std::cout << "Edge 1: " <<  tFace->edge(0)->id() << ", Orientation : " << tFace->edge_orientation(0)<< std::endl;
        std::cout << "Edge 2: " <<  tFace->edge(1)->id() << ", Orientation : " << tFace->edge_orientation(1)<< std::endl;
        std::cout << "Edge 3: " <<  tFace->edge(2)->id() << ", Orientation : " << tFace->edge_orientation(2)<< std::endl;
        std::cout << "Number of elements: " << tFace->number_of_elements() << std::endl ;
        std::cout << std::endl ;
    }

    for(mesh::Element * tElement : tMesh->elements())
    {
        std::cout << "Element ID: " <<  tElement->id() << std::endl ;
        std::cout << "Number of faces : " << tElement->number_of_faces() <<std::endl;
        std::cout << "Face 1: " <<  tElement->face(0)->id() << ", Orientation : " << (tElement->face(0)->master()->id() == tElement->id()?1:0)<< std::endl;
        std::cout << "Face 2: " <<  tElement->face(1)->id() << ", Orientation : " << (tElement->face(1)->master()->id() == tElement->id()?1:0)<< std::endl;
        std::cout << "Face 3: " <<  tElement->face(2)->id() << ", Orientation : " << (tElement->face(2)->master()->id() == tElement->id()?1:0)<< std::endl;
        std::cout << "Face 4: " <<  tElement->face(3)->id() << ", Orientation : " << (tElement->face(3)->master()->id() == tElement->id()?1:0)<< std::endl;
        std::cout << std::endl ;
    }


    //mesh::ExodusWriter tWriter( tMesh );
    //tMesh->save( "mesh.vtk" );
    tMesh->save( "mesh.exo" );
    tMesh->save_faces( "faces.exo");

    delete tMesh ;

    return  gComm.finalize();
}