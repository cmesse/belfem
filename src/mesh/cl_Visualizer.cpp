//
// Created by christian on 5/7/24.
//

#include "cl_Visualizer.hpp"

#ifdef BELFEM_VTK
#if BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

// datatypes
#include <vtkFloatArray.h>
#include <vtkPointData.h>

// element types
#include <vtkLine.h>

// glyphs
#include <vtkGlyph3D.h>
#include <vtkConeSource.h>


#if BELFEM_GCC
#pragma GCC diagnostic pop
#endif

#endif

namespace belfem
{
    Visualizer::Visualizer( Mesh * aMesh ) :
        mMesh( aMesh ),
        mRank( comm_rank() )
    {
        if( mRank == 0 )
        {
            mOptions = new VisualizerOptions();

            this->create_points();
            this->create_edges() ;
            this->create_renderer() ;
        }
    }

//------------------------------------------------------------------------------

    Visualizer::~Visualizer()
    {
        if( mRank == 0 )
        {
            delete mOptions ;
        }
    }


//------------------------------------------------------------------------------

    void
    Visualizer::create_points()
    {
#ifdef BELFEM_VTK

        Cell< mesh::Node * > & tNodes = mMesh->nodes();

        mPoints = vtkSmartPointer<vtkPoints>::New();


        mPoints->SetNumberOfPoints( tNodes.size() );

        vtkIdType tCount = 0 ;

        for( mesh::Node * tNode : tNodes )
        {
            tNode->set_index( tCount );
            mPoints->SetPoint( tCount++, tNode->x(), tNode->y(), tNode->z() );
        }

#endif
    }

//------------------------------------------------------------------------------

    void
    Visualizer::create_edges()
    {
#ifdef BELFEM_VTK
        mEdges          = vtkSmartPointer<vtkCellArray>::New();
        mEdgeDirections = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints>     tEdgeCenters    = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkFloatArray> tEdgeDirections = vtkSmartPointer<vtkFloatArray>::New();

        tEdgeDirections->SetNumberOfComponents(3);
        tEdgeDirections->SetName("Directions");

        real tMeanLength = 0 ;

        for( mesh::Edge * tEdge : mMesh->edges() )
        {
            real tDirection[3];
            real tPoint[3];

            vtkSmartPointer<vtkLine> tLine = vtkSmartPointer<vtkLine>::New() ;

            tLine->GetPointIds()->SetId( 0, tEdge->node( 0 )->index() );
            tLine->GetPointIds()->SetId( 1, tEdge->node( 1 )->index() );

            tDirection[0] = tEdge->node( 1 )->x() - tEdge->node( 0 )->x() ;
            tDirection[1] = tEdge->node( 1 )->y() - tEdge->node( 0 )->y() ;
            tDirection[2] = tEdge->node( 1 )->z() - tEdge->node( 0 )->z() ;

            real tLength = std::sqrt(
                         tDirection[0] * tDirection[0]
                      +  tDirection[1] * tDirection[1]
                      +  tDirection[2] * tDirection[2] );

            tDirection[0] /= tLength ;
            tDirection[1] /= tLength ;
            tDirection[2] /= tLength ;

            tMeanLength += tLength ;

            tPoint[0] = 0.5*( tEdge->node( 1 )->x() + tEdge->node( 0 )->x() ) ;
            tPoint[1] = 0.5*( tEdge->node( 1 )->y() + tEdge->node( 0 )->y() ) ;
            tPoint[2] = 0.5*( tEdge->node( 1 )->z() + tEdge->node( 0 )->z() ) ;

            tEdgeCenters->InsertNextPoint( tPoint );
            tEdgeDirections->InsertNextTuple( tDirection );

            mEdges->InsertNextCell( tLine );
        }

        tMeanLength /= mMesh->edges().size() ;

        mEdgeDirections->SetPoints( tEdgeCenters );
        mEdgeDirections->GetPointData()->SetVectors( tEdgeDirections );

        mOptions->ConeHeight = tMeanLength * 0.1 ;
        mOptions->ConeRadius = tMeanLength * 0.01 ;

#endif
    }

//------------------------------------------------------------------------------

    void
    Visualizer::create_renderer()
    {
#ifdef BELFEM_VTK

        // ------------------------- create the cones

        vtkSmartPointer<vtkConeSource> tConeSource = vtkSmartPointer<vtkConeSource>::New();
        tConeSource->SetHeight(mOptions->ConeHeight );
        tConeSource->SetRadius(mOptions->ConeRadius );
        tConeSource->SetResolution( mOptions->ConeResolution );
        tConeSource->SetAngle( mOptions->ConeAngle  );

        vtkSmartPointer<vtkGlyph3D> tGlyph3D = vtkSmartPointer<vtkGlyph3D>::New();
        tGlyph3D->SetSourceConnection(tConeSource->GetOutputPort());
        tGlyph3D->SetInputData(mEdgeDirections);
        tGlyph3D->SetVectorModeToUseVector();
        tGlyph3D->SetScaleModeToScaleByVector();
        tGlyph3D->SetScaleFactor(1);

        // Create a mapper and actor for the glyphs
        vtkSmartPointer<vtkPolyDataMapper> tGlyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        tGlyphMapper->SetInputConnection(tGlyph3D->GetOutputPort() );

        vtkSmartPointer<vtkActor> tGlyphActor = vtkSmartPointer<vtkActor>::New();
        tGlyphActor->SetMapper( tGlyphMapper );
        tGlyphActor->GetProperty()->SetColor( mOptions->ConeColor[0], mOptions->ConeColor[1], mOptions->ConeColor[2] );

        // ------------------------- end create cones

        mData = vtkSmartPointer<vtkPolyData>::New();
        mData->SetPoints(mPoints);
        mData->SetLines(mEdges); // Ensure lines are being used

        mMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mMapper->SetInputData(mData);

        mActor = vtkSmartPointer<vtkActor>::New();
        mActor->SetMapper(mMapper);
        mActor->GetProperty()->SetColor(
                mOptions->EdgeColor[0],
                mOptions->EdgeColor[1],
                mOptions->EdgeColor[2] );

        mActor->GetProperty()->SetLineWidth(mOptions->EdgeWidth  );

        mRenderer = vtkSmartPointer<vtkRenderer>::New();
        mRenderer->AddActor(mActor );
        mRenderer->AddActor( tGlyphActor );

        mRenderer->SetBackground( mOptions->BackgroundColor[0],
                                  mOptions->BackgroundColor[1],
                                  mOptions->BackgroundColor[2] );

        //mRenderer->AutomaticLightCreationOn(); // Ensure there is light
        mRenderer->ResetCamera(); // Reset the camera to view the objects

        mRenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        mRenderWindow->AddRenderer(mRenderer);
        mRenderWindow->SetSize(mOptions->WindowSize[0], mOptions->WindowSize[1]); // Set window size

        mRenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        mRenderWindowInteractor->SetRenderWindow(mRenderWindow );

        mRenderWindow->Render();
        mRenderWindowInteractor->Start();
#endif
    }

//------------------------------------------------------------------------------
}