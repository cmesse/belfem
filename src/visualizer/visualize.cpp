/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <cstdio>

#include "banner.hpp"
#include "cl_Arguments.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "cl_Mesh.hpp"
#include "cl_VTK_MeshView.hpp"

#ifdef BELFEM_VTK
#include "vtkSmartPointer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#endif

#include "vtktypes.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    print_banner( );

    // usage: visualize <mesh.msh> [screenshot.png]
    //   <mesh>         -> that mesh, interactive window
    //   <mesh> <png>   -> render one frame off-screen to <png> and exit
    Arguments tArguments( argc, argv );
    const Cell< string > & tArgs = tArguments.data();   // tArgs( 0 ) = program

    if ( tArgs.size() < 2 )
    {
        std::fprintf( stderr, "usage: visualize <mesh> [screenshot.png]\n" );
        return gComm.finalize();
    }

    const string tPath = tArgs( 1 );

    const string tScreenshot = ( tArgs.size() > 2 ) ? tArgs( 2 ) : "";

    // load the mesh serially (full mesh on this process).
    // connectivity computation is off: the visualizer only needs node
    // coordinates and per-element/facet connectivity, not facet->master links.
    Mesh * tMesh = new Mesh( tPath, 0, false, false );
#ifdef BELFEM_VTK
    // build the VTK actors from the mesh
    // (fully qualified: VTK 9 also defines a global ::vtk namespace)
    belfem::vtk::MeshView tVtkMesh( tMesh );


    // rendering is rank-0 only
    if ( gComm.rank() == 0 )
    {
        belfem::vtk::Renderer tRenderer = belfem::vtk::Renderer::New();
        tRenderer->SetBackground( 0.1, 0.1, 0.15 );

        for ( belfem::vtk::Actor tActor : tVtkMesh.actors() )
        {
            tRenderer->AddActor( tActor );
        }

        belfem::vtk::RenderWindow tWindow = belfem::vtk::RenderWindow::New();
        tWindow->AddRenderer( tRenderer );
        tWindow->SetSize( 1280, 800 );
        tWindow->SetWindowName( "BELFEM Visualizer" );

        if ( ! tScreenshot.empty() )
        {
            // off-screen smoke test: render one frame to a PNG and exit
            tWindow->SetOffScreenRendering( 1 );
            tRenderer->ResetCamera();
            tWindow->Render();

            auto tToImage = vtkSmartPointer< vtkWindowToImageFilter >::New();
            tToImage->SetInput( tWindow );
            tToImage->Update();

            auto tWriter = vtkSmartPointer< vtkPNGWriter >::New();
            tWriter->SetFileName( tScreenshot.c_str() );
            tWriter->SetInputConnection( tToImage->GetOutputPort() );
            tWriter->Write();
        }
        else
        {
            belfem::vtk::RenderWindowInteractor tInteractor =
                belfem::vtk::RenderWindowInteractor::New();
            tInteractor->SetRenderWindow( tWindow );

            auto tStyle = vtkSmartPointer< vtkInteractorStyleTrackballCamera >::New();
            tInteractor->SetInteractorStyle( tStyle );

            tRenderer->ResetCamera();
            tWindow->Render();
            tInteractor->Start();
        }
    }
#endif

    delete tMesh ;

    return gComm.finalize();
}
