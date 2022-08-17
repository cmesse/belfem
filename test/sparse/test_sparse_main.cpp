//
// Created by Christian Messe on 06.11.19.
//


#include <gtest/gtest.h>
#include "cl_Communicator.hpp"

#include "cl_Logger.hpp"

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 5 );

belfem::Cell< belfem::graph::Vertex * > gGraph;

//------------------------------------------------------------------------------

void
create_graph()
{
    gGraph.set_size( 5, nullptr );

    for( uint k=0; k<5; ++k )
    {
        belfem::graph::Vertex * tVertex = new belfem::graph::Vertex();
            tVertex->set_index( k );
            tVertex->set_id( k+1 );
            gGraph( k ) = tVertex;
    }

    // link graph
    belfem::graph::Vertex * tVertex = nullptr;

    // Vertex 1
    tVertex = gGraph( 0 );
    tVertex->init_vertex_container( 3 );
    tVertex->insert_vertex( gGraph( 0) );
    tVertex->insert_vertex( gGraph( 1 ) );
    tVertex->insert_vertex( gGraph( 3 ) );

    // Vertex 2
    tVertex = gGraph( 1 );
    tVertex->init_vertex_container( 2 );
    tVertex->insert_vertex( gGraph( 0 ) );
    tVertex->insert_vertex( gGraph( 1 ) );


    // Vertex 3
    tVertex = gGraph( 2 );
    tVertex->init_vertex_container( 3 );
    tVertex->insert_vertex( gGraph( 2 ) );
    tVertex->insert_vertex( gGraph( 3 ) );
    tVertex->insert_vertex( gGraph( 4 ) );

    // Vertex 4
    tVertex = gGraph( 3 );
    tVertex->init_vertex_container( 3 );
    tVertex->insert_vertex( gGraph( 0 ) );
    tVertex->insert_vertex( gGraph( 2 ) );
    tVertex->insert_vertex( gGraph( 3 ) );

    // Vertex 2
    tVertex = gGraph( 4 );
    tVertex->init_vertex_container( 2 );
    tVertex->insert_vertex( gGraph( 1 ) );
    tVertex->insert_vertex( gGraph( 4 ) );
}

//------------------------------------------------------------------------------

void
destroy_graph()
{
    for ( auto tVertex : gGraph )
    {
        delete tVertex;
    }
}

//------------------------------------------------------------------------------

int
main( int    argc,
      char * argv[] )
{
    // create communicator
    gComm = belfem::Communicator( argc, argv );

    // start test session
    testing::InitGoogleTest( &argc, argv );

    int aResult = 0;
    if( gComm.rank() == 0 )
    {
        // create the graph
        create_graph();

        // run the tests
        aResult = RUN_ALL_TESTS();

        // tidy up memory
        destroy_graph();
    }

    // close communicator
    gComm.finalize();

    // return the test result
    return aResult;
}

