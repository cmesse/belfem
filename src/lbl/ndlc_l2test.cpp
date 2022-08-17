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

#include "scratch/cl_IWG_Nedelec_L2_Node2Edge_h.hpp"
#include "scratch/cl_IWG_Nedelec_L2_Edge2Node_h.hpp"
#include "scratch/cl_IWG_Nedelec_L2_Edge2Node_curlh.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;
Logger       gLog( 5 );

void
create_demo_data( Mesh & aMesh )
{
    // create the fields
    Vector< real > & tAx = aMesh.create_field( "hx");
    Vector< real > & tAy = aMesh.create_field( "hy");

    Vector< real > & tC = aMesh.create_field( "curl0");

    // get the nodes
    Cell< mesh::Node * > & tNodes = aMesh.nodes() ;

    // initialize counter
    index_t tCount = 0 ;

    // loop over all nodes
    for( mesh::Node * tNode : tNodes )
    //for( index_t k=0; k<tNodes.size(); ++k )
    {
        real tX = tNode->x();
        real tY = tNode->y();

        // compute distance to center
        real tR = std::sqrt( tX * tX + tY * tY );

        if ( abs( tR ) < 1e-12 )
        {
            tAx( tCount ) = 0.0;
            tAy( tCount ) = 0.0;
            tC( tCount ) = 0.0;
        }
        else
        {
            tAx( tCount ) = -tY / tR;
            tAy( tCount ) = tX / tR;


            tC( tCount ) = 1 / tR;
        }

        // tAx( tCount ) = 1 ;
       // tAy( tCount ) = 1 ;

        ++tCount ;
    }
}

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    print_banner();

    // open mesh
     Mesh tMesh( "square.msh");
    //Mesh tMesh( "tri6.msh");

    tMesh.create_edges() ;

    if( comm_rank() == 0 )
    {
        create_demo_data( tMesh );

    }
    KernelParameters tParams( &tMesh );


    tParams.set_field_dimensions( { 0, 2, 1 }, { 1, 0, 0 } );

    if ( tMesh.max_element_order() == 2 )
    {

        tParams.enforce_linear( { 1, 1, 1 } );
    }

    tParams.set_block_integration_orders({ 9 });

    // tParams.set_integration_scheme( IntegrationScheme::GAUSSCLASSIC );

    Kernel tKernel( &tParams, true );

    IWG_Nedelec_Node2Edge_h tL2a( tMesh.number_of_dimensions() ) ;

    // now we grab this one field into the pointer
    Field * tField0 = tKernel.field( 0 );

    // we want to solve the plane stress equation on the field ...
    tField0->set_integrated_weak_governing_equation( & tL2a );

    tField0->set_solver( SolverType::UMFPACK );

    comm_barrier() ;

    tField0->initialize_jacobian();
    tField0->distribute_fields( {"hx","hy"} );
    tField0->compute_jacobian();
    tField0->compute_rhs();
    tField0->solve();

    if( comm_rank() == 0 )
    {
        // backup field
        Vector< real > & tAx = tMesh.create_field( "h0x" );
        Vector< real > & tAy = tMesh.create_field( "h0y" );
        tAx = tMesh.field_data( "hx" );
        tAy = tMesh.field_data( "hy" );

        // tMesh.field_data("edge_a").print("edge_a");
    }


    comm_barrier() ;

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

    tMesh.save("mesh.exo");

    // close communicator
    return gComm.finalize();
}