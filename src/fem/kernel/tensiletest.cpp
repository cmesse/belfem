//
// Created by Christian Messe on 26.10.19.
//


#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "Mesh_Enums.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"

#include "cl_IWG_SimpleDiffusion.hpp"
#include "cl_IWG_PlaneStress.hpp"
#include "cl_IWG_2DGradient.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "fn_FEM_mises_planestress.hpp"

#include "cl_Mesh_HDF5Reader.hpp"

#include "banner.hpp"
#include "commtools.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;
Logger       gLog( 3 );

/**
 * this is an example problem that expains how the finite-element core works
 */
int main( int    argc,
          char * argv[] )
{
    // every program in BELFEM begins with creating a communicator for the parallel
    // communication. If we run without MPI, the communicator does noting
    gComm = Communicator( argc, argv );

    // if we wish, we can display the banner
    print_banner();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // step 1: read the mesh and prepare the Kernel
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // we begin with creating a new mesh that we load from the gmsh file
    Mesh * tMesh = new Mesh( "zugprobe.msh" );

    // the parametor object is used to init the Kernel
    KernelParameters tParams( tMesh );

    // with the parameters object set, we create the kernel
    Kernel tKernel( &tParams );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // step 2: solve the plane stress problem
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // IWG stands for Integrated-Weak-Governing Equation object.
    // In this case, the IWG contains the Linear elastic equation for plane stress
    // now we grab this one field into the pointer
    IWG * tIWG = tKernel.create_equation( IwgType::PlaneStress );

    tIWG->select_sidesets( { 1, 2, 3, 4, 5, 6 } );

    DofManager * tField = tKernel.create_field( tIWG );

    // using the MUMPS solver
    // we can also use SolverType::UMFPACK (default), SolverType::MUMPS (needs MPI), if we have MKL, SolverType::PARDISO
    tField->set_solver( SolverType::UMFPACK );
    tField->solver()->set_petsc( Preconditioner::GAMG, KrylovMethod::GMRES, 1e-8 );

    // before we run the equation, we want to fix a brearing in y-direction
    tField->bearing( 1 )->impose_dirichlet( 0.0, 1 );

    // and the rest of this sidesets
    tField->sideset( 5 )->impose_dirichlet( 0.0, 0 );
    tField->sideset( 3 )->impose_dirichlet( 0.001, 0 );

    // this material only has one block, we use Aluminum
    tField->block( 1 )->set_material( MaterialType::Aluminum );

    // now we wait until all procs have arrived here
    comm_barrier();

    // now we compute the Jacobian
    tField->compute_jacobian();

    // and forward it to the solver
    tField->solve();

    // sava data to hdf file for debugging
    tField->save_system("memdump.hdf5");

    // now we wait until all procs have arrived here
    comm_barrier();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // step 3: compute the derivatives
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // once we have u and v, we need dudx, dudv, dvdx and dvdy

    IWG * tIWG2 = tKernel.create_equation( IwgType::Gradient2D );
    DofManager * tField2 = tKernel.create_field( tIWG2 );
    // the following command creates the adjency graph for the Jacobian
    tField2->set_solver( SolverType::UMFPACK );
    tField2->compute_jacobian();
    tField2->compute_rhs();
    tField2->solve();


    // now we wait until all procs have arrived here
    comm_barrier();

    // with the derivatives computed, we can compute the van mises stress
    mises_planestress( tField2 );

    // only the master proc knows the full mesh.
    if( comm_rank() == 0 )
    {
        // the following command stores the result in a HDF5 file
        // tMesh->save( "Mesh.hdf5" );

        // this command stores the file in a format that can be read by ParaVIEW
        tMesh->save( "Mesh.exo" );
    }
    else
    {
        //tField->mesh()->save( sprint( "submesh_%u.exo", ( unsigned int ) comm_rank() ) );
    }

    // finally, we can delete the mesh
    delete tMesh;

    // and tell the communicator to shut the program down
    return gComm.finalize();
}
