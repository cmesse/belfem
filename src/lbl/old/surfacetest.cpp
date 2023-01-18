//
// Created by Christian Messe on 30.06.21.
//

//
// Created by Christian Messe on 29.06.21.
//

//
// Created by christian on 6/9/21.
//
#include <iostream>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"

#include "scratch/cl_IWG_Nedelec_h_static.hpp"
#include "scratch/cl_IWG_Nedelec_L2_Edge2Node_h.hpp"
#include "scratch/cl_IWG_Nedelec_L2_Edge2Node_curlh.hpp"

#include "fn_Mesh_compute_surface_normals.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    print_banner();

    // open mesh
    Mesh tMesh( "square.msh" );

    // scale mesh to meters
    // tMesh.scale_mesh( 0.001 );

    mesh::compute_surface_normals( &tMesh, { 1, 2, 3, 4 }, GroupType::SIDESET );

    KernelParameters tParams( &tMesh );

    Kernel tKernel( &tParams, true );

    Field * tField = tKernel.field( 0 );

    SideSet * tSideSet = tField->sideset( 1 );

    // length
    real tLength = 0 ;


    // close communicator
    return gComm.finalize();
}
