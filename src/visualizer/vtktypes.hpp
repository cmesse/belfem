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

#ifndef BELFEM_VTKTYPES_HPP
#define BELFEM_VTKTYPES_HPP

#ifdef BELFEM_VTK
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkTexture.h"
#include "vtkCamera.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkPolyDataMapper.h"
#include "vtkDataSetMapper.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkUnstructuredGrid.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#endif
namespace belfem
{
    namespace vtk
    {
#ifdef BELFEM_VTK
        typedef vtkSmartPointer< vtkActor >                  Actor ;
        typedef vtkSmartPointer< vtkTexture >                Texture ;
        typedef vtkSmartPointer< vtkCamera >                 Camera ;
        typedef vtkSmartPointer< vtkPolyData >               PolyData ;
        typedef vtkSmartPointer< vtkTransform >              Transform ;
        typedef vtkSmartPointer< vtkPolyDataMapper >         Mapper ;
        typedef vtkSmartPointer< vtkDataSetMapper >          DataSetMapper ;
        typedef vtkSmartPointer< vtkPoints >                 Points ;
        typedef vtkSmartPointer< vtkLine >                   Line ;
        typedef vtkSmartPointer< vtkCellArray >              CellArray ;
        typedef vtkSmartPointer< vtkUnstructuredGrid >       UnstructuredGrid ;
        typedef vtkSmartPointer< vtkRenderer >               Renderer ;
        typedef vtkSmartPointer< vtkRenderWindow >           RenderWindow ;
        typedef vtkSmartPointer< vtkRenderWindowInteractor > RenderWindowInteractor ;

#else
        typedef void Actor ;
        typedef void Texture ;
        typedef void Camera ;
        typedef void PolyData ;
        typedef void Transform ;
        typedef void Mapper ;
        typedef void DataSetMapper ;
        typedef void Points ;
        typedef void Line ;
        typedef void CellArray ;
        typedef void UnstructuredGrid ;
        typedef void Renderer ;
        typedef void RenderWindow ;
        typedef void RenderWindowInteractor ;
#endif
    }
}
#endif //BELFEM_VTK_TYPES_HPP
