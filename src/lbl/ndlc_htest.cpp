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
    Mesh tMesh( "coil.msh");

    // scale mesh to meters
    tMesh.scale_mesh( 0.001 );

    mesh::compute_surface_normals( &tMesh, { 1, 2, 3, 4, 5, 6, 7, 8 }, GroupType::SIDESET );

    tMesh.create_edges() ;

    KernelParameters tParams( &tMesh );


    tParams.set_field_dimensions( { 0 , 2, 1 }, { 1, 0, 0 } );

    if ( tMesh.max_element_order() == 2 )
    {

        tParams.enforce_linear( { 1, 1, 1 } );
    }


    tParams.set_block_integration_orders({ 5 } );


    Kernel tKernel( &tParams, true );

    IWG_Nedelec_h_static tIWG( tMesh.number_of_dimensions() ) ;


    // now we grab this one field into the pointer
    Field * tField0 = tKernel.field( 0 );

    // we want to solve the plane stress equation on the field ...
    tField0->set_integrated_weak_governing_equation( & tIWG );

    // set the penalty factor
    tIWG.set_psi( 1.0e9 );


    tIWG.set_wetted_sidesets( { 1, 2, 3, 4, 5, 6, 7, 8 } );

    // impose current densities
    tMesh.field("jz")->fill( { 2, 3 },  5.0930e+08 );

    tField0->set_solver( SolverType::UMFPACK );
    // tField0->solver()->set_petsc( Preconditioner::ASM, KrylovMethod::GMRES );

    comm_barrier() ;

    tField0->initialize_jacobian();
    tField0->reset_jacobian();
    tField0->compute_weak_symmetry({ 1, 4, 5, 7, 8 }, 0 );

    tField0->compute_jacobian( false );
    tField0->compute_surface_currents( { 2, 3 } );

    tField0->solve();


    IWG_Nedelec_Edge2Node_h tL2b( tMesh.number_of_dimensions() ) ;

    // now we grab the second field into the pointer
    Field * tField1 = tKernel.field( 1 );

    tField1->set_integrated_weak_governing_equation( & tL2b );
    tField1->set_solver( SolverType::UMFPACK );

    tField1->initialize_jacobian();
    tField1->compute_jacobian();
    tField1->compute_rhs();
    tField1->solve();

    comm_barrier() ;

    Vector< real > & tNorm = tMesh.create_field( "norm_h");

    Vector< real > & aX = tMesh.field_data( "hx");
    Vector< real > & aY = tMesh.field_data( "hy");

    for( index_t k=0; k<tMesh.number_of_nodes(); ++k )
    {
        tNorm( k ) = std::sqrt( aX( k ) * aX( k ) + aY( k ) * aY( k ) );
    }

    IWG_Nedelec_Edge2Node_curlh tL2c( tMesh.number_of_dimensions() ) ;

    // now we grab the second field into the pointer
    Field * tField2 = tKernel.field( 2 );

    tField2->set_integrated_weak_governing_equation( & tL2c );
    tField2->set_solver( SolverType::UMFPACK );

    tField2->initialize_jacobian();
    tField2->compute_jacobian();
    tField2->compute_rhs();
    tField2->solve();

    comm_barrier() ;

    tMesh.save( "mesh.exo");

    // close communicator
    return gComm.finalize();
}