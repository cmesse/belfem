/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"

#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IwgFactory.hpp"

#include "fn_read_node_tmatrix_from_hdf5.hpp"

using namespace belfem;
using namespace fem ;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );


    // load the mesh
    Mesh tMesh( "hanging.msh" );

    if( gComm.rank() == 0 )
    {
        mesh::read_node_tmatrix_from_hdf5( &tMesh, "nodedata.hdf5" );
        tMesh.collect_hanging_basis() ;
    }

    KernelParameters tParams( tMesh );

    Kernel tKernel( &tParams );

    IWG * tIWG = tKernel.create_equation( IwgType::Poisson );

    tIWG->select_block( 5 );
    tIWG->select_sidesets( { 2, 4 } );


    DofManager * tField = tKernel.create_field( tIWG );

    tField->set_solver( SolverType::UMFPACK );



    // tField->block( 5 )->set_material( MaterialType::Aluminum );



    // initialize starting conditions


    tField->sideset( 2 )->impose_dirichlet( 3.0 );
    tField->sideset( 4 )->impose_dirichlet( 1.0 );

    tField->initialize() ;

    tField->compute_jacobian();

    tField->solve();

    tMesh.save( "hanging.exo");

    return gComm.finalize();
}