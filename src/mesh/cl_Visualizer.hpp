//
// Created by christian on 5/7/24.
//

#ifndef BELFEM_CL_VISUALIZER_HPP
#define BELFEM_CL_VISUALIZER_HPP

#ifdef BELFEM_VTK

#if BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#if BELFEM_GCC
#pragma GCC diagnostic pop
#endif


#endif

#include "commtools.hpp"
#include "constants.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    struct VisualizerOptions
    {
        uint WindowSize[3]      = { 800, 600 };

        real BackgroundColor[3] = { 0.1, 0.1, 0.1 } ;
        real EdgeColor[3]       = { 0, 1, 1 };
        real ConeColor[3]       = { 1, 1, 0 };

        real EdgeWidth          = 2 ;

        real ConeHeight         = 0.1 ;
        real ConeRadius         = 0.05 ;
        uint ConeResolution     = 12 ;
        real ConeAngle          = 20 ;

    public:

        VisualizerOptions() = default ;

        ~VisualizerOptions() = default ;

    };

    class Visualizer
    {
        Mesh     *   mMesh;
        const proc_t mRank;
        VisualizerOptions * mOptions = nullptr ;

#ifdef BELFEM_VTK
        vtkSmartPointer<vtkPoints>         mPoints ;
        vtkSmartPointer<vtkCellArray>      mEdges ;
        vtkSmartPointer<vtkPolyData>       mEdgeDirections ;

        vtkSmartPointer<vtkPolyData>       mData ;
        vtkSmartPointer<vtkPolyDataMapper> mMapper ;
        vtkSmartPointer<vtkActor>          mActor ;

        //vtkSmartPointer<vtkCamera>       mCamera ;
        //vtkSmartPointer<vtkLight>        mAmbientLight ;
        vtkSmartPointer<vtkRenderer>       mRenderer ;

        vtkSmartPointer<vtkRenderWindow> mRenderWindow ;
        vtkSmartPointer<vtkRenderWindowInteractor> mRenderWindowInteractor ;

#endif
    public:

        Visualizer( Mesh * aMesh );

        ~Visualizer() ;

    private:

        void
        create_points();

        void
        create_edges();

        void
        create_renderer();

    };
}

#endif //BELFEM_CL_VISUALIZER_HPP
