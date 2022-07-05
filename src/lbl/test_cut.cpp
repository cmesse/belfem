//
// Created by christian on 9/15/21.
//
#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"

#include "cl_Mesh.hpp"
#include "cl_Mesh_Scissors.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"

using namespace belfem ;


Communicator gComm;
Logger       gLog( 3 );
//------------------------------------------------------------------------------


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );
    Mesh tMesh = Mesh("coilcut.msh") ;
    tMesh.scale_mesh( 0.001 );

    mesh::Scissors tScissors( &tMesh, { 4, 5, 6 } );

    tScissors.cut( 26, 6, 4 );
    tScissors.cut( 27, 6, 5 );

    /*tScissors.cut( { 7, 8 } , 4, 2 );
    tScissors.cut( 22 , 4, 2 );
    tScissors.cut( { 6, 1, 2 } , 4, 1 );
    tScissors.cut( 23 , 4, 2 );
    tScissors.cut( { 11, 12 } , 4, 3 ); */

    // tScissors.cut( 24 , 4, 2 );

    tScissors.finalize();

    // the parametor object is used to init the Kernel
    fem::KernelParameters tParams( tMesh );

    // with the parameters object set, we create the kernel
    fem::Kernel tKernel( &tParams );


    string tString = "boxcut_" + std::to_string( comm_rank() ) + ".vtk";

    tKernel.mesh()->save(tString );

    tMesh.save("mycut.vtk");

    // get facet
    //mesh::Facet * tFacet = tMesh.sideset( 5 )->facet_by_index( 0 );

    // get element
    //mesh::Element * tElement =  tFacet->element() ;
//------------------------------------------------------------------------------

    // close communicator
    return gComm.finalize();
}