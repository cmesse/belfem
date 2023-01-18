//
// Created by christian on 9/17/21.
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
#include "cl_Mesh_TapeRoller.hpp"

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

    // read the input file
    InputFile tInputFile( "input.txt" );

    // create a new factory
    MaxwellFactory * tFactory = new MaxwellFactory( tInputFile );

    // get the mesh
    Mesh * tMesh = tFactory->mesh() ;

    if( comm_rank() == 0 )
    {
        //Vector< id_t > tSideSets = { 9, 10, 11, 12 };
        Vector< id_t > tSideSets = { 11, 12, 13, 14 };
        //tMesh->create_edges( false, {}, tSideSets );
        //tMesh->create_faces( false, {}, tSideSets ) ) ;

        // count elements on these sidesets
        tMesh->unflag_all_faces();
        for ( id_t tID: tSideSets )
        {
            tMesh->sideset( tID )->flag_all_facets();
        }
        index_t tCount = 0;
        for ( mesh::Facet * tFacet: tMesh->facets())
        {
            if ( tFacet->is_flagged())
            {
                ++tCount;
            }
        }
        Cell< mesh::Element * > tElements( tCount, nullptr );
        tCount = 0;
        for ( mesh::Facet * tFacet: tMesh->facets())
        {
            if ( tFacet->is_flagged())
            {
                tElements( tCount++ ) = tFacet->element();
            }
        }
        tMesh->finalize_edges( tElements );
        //tMesh->finalize_faces();

        mesh::TapeRoller * tTapeRoller = new mesh::TapeRoller( tMesh, 3 );

        tTapeRoller->add_sidesets( tSideSets );
        tTapeRoller->run();

        delete tTapeRoller ;
    }

    comm_barrier() ;

    Kernel *  tKernel = tFactory->kernel() ;

    if( comm_rank() == 0 )
    {
        tKernel->mesh()->save( "mesh.exo");
    }
    else
    {
        tKernel->mesh()->save( "submesh.exo");
    }



    // create the IWG

    // tMesh->save( "test.exo");

    delete tFactory ;
    delete tKernel ;
    delete tMesh ;

    // close communicator
    return gComm.finalize();
}