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

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"

#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IwgFactory.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // load the mesh
    Mesh tMesh( "dipole.msh");

    KernelParameters tParams( tMesh );

    Kernel tKernel( &tParams );

    IWG * tIWG = tKernel.create_equation( IwgType::Poisson );

    tIWG->select_blocks( { 3 } );
    tIWG->select_sidesets( { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } );


    DofManager * tField = tKernel.create_field( tIWG );

    tField->set_solver( SolverType::PETSc );

    // set potential for conductors
    for( id_t k=1; k<=8; ++k )
    {
        tField->sideset( k )->impose_dirichlet( 0.0 );
    }

    // set potential for boundary
    for( id_t k=9; k<=12; ++k )
    {
        tField->sideset( k )->impose_dirichlet( 1.0 );
    }

    // compute and solve system
    tField->compute_jacobian();
    tField->solve() ;


    tMesh.save( "dipole.exo");

    return gComm.finalize();
}
