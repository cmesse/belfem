//
// Created by christian on 7/7/21.
//

#include <iostream>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Profiler.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"

#include "cl_IwgFactory.hpp"

#include "fn_Mesh_compute_surface_normals.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_IWG_StationaryHeatConduction.hpp"
#include "cl_Bitset.hpp"
#include "fn_FEM_mises_PlaneStress.hpp"
//#define protected public
//#define private   public
#include "cl_IWG_Maxwell_N2HA.hpp"
//#undef private
//#undef protected
#include "fn_Mesh_compute_volume.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;
Logger       gLog( 4 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    // open mesh
    //Mesh tMesh( "test.msh");
    Mesh tMesh( "coil.msh");

    // scale mesh to meters
    tMesh.scale_mesh( 0.001 );
    tMesh.create_edges() ;

    KernelParameters tParams( &tMesh );

    Kernel tKernel( &tParams, false );

    tMesh.save( "mesh.exo");

    IWG_Maxwell_N2HA * tIWG = reinterpret_cast< IWG_Maxwell_N2HA *  >
            ( tKernel.create_equation( IwgType::MAXWELL_N2HA ) );

    //tIWG->set_block( 2, DomainType::SuperConductor );
    //tIWG->set_block( 1, DomainType::Air );
    //tIWG->set_sideset( 1,  DomainType::Interface );

    tIWG->set_block( 1, DomainType::SuperConductor );
    tIWG->set_block( { 2, 3 }, DomainType::Coil );
    tIWG->set_block( 4, DomainType::Air );
    tIWG->set_sideset( { 1, 2, 3, 4 },  DomainType::Interface );

    tIWG->set_timestep( 1.0 );


    DofManager * tField = tKernel.create_field( tIWG );

    // tField->block( 1 )->set_integration_order( 9 );

    tField->initialize() ;
    tField->zero() ;

    tField->compute_jacobian_and_rhs( true );

    // surface of circle
    real tA = mesh::compute_volume( &tMesh, { 2 } ) ;

    std::cout << "surface " << tA  << std::endl ;

    //tIWG->link_to_group( tField->sideset( 1 ) );

    // grab test element
    //Element * tElement = tField->sideset( 1 )->elements()(0);

    //Matrix< real > tJ( 6, 6 );
    //Vector< real > tX( 6 );

    //tIWG->compute_matrix_ah( tElement, tJ, tX );

    /*tIWG->link_to_group( tField->block( 1 ) );
    Element * tElement = tField->block( 1 )->elements()(0);
    Matrix< real > tJ( 6, 6 );
    Vector< real > tX( 6 );

    tIWG->print_dofs( tElement );

    tIWG->compute_jacobian_and_rhs( tElement, tJ, tX ); */

    // close communicator
    return gComm.finalize();
}