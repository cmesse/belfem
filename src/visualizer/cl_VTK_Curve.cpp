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

#include "cl_VTK_Curve.hpp"
#include "assert.hpp"

#ifdef BELFEM_VTK
#include "vtkCellArray.h"
#include "vtkProperty.h"
#endif

namespace belfem
{
    namespace vtk
    {
        Curve::Curve(
            const Matrix< real > & aPoints,
            real                   aLineWidth,
            bool                   aClosed ) :
            mPoints( Points::New() ),
            mPolyData( PolyData::New() ),
            mMapper( Mapper::New() ),
            mActor( Actor::New() )
        {
            BELFEM_ERROR( aPoints.n_rows() == 3,
                "vtk::Curve expects a 3 x N point matrix, got %u rows",
                ( unsigned int ) aPoints.n_rows() );

            const index_t tNumPoints = aPoints.n_cols();

            mPoints->Allocate( tNumPoints );
            for ( index_t k = 0; k < tNumPoints; ++k )
            {
                mPoints->InsertNextPoint( aPoints( 0, k ),
                                          aPoints( 1, k ),
                                          aPoints( 2, k ) );
            }

            // a single polyline cell threading all points in order; when
            // closed, the first point is repeated to seal the loop
            CellArray tLines = CellArray::New();
            tLines->InsertNextCell( aClosed ? tNumPoints + 1 : tNumPoints );
            for ( index_t k = 0; k < tNumPoints; ++k )
            {
                tLines->InsertCellPoint( k );
            }
            if ( aClosed )
            {
                tLines->InsertCellPoint( 0 );
            }

            mPolyData->SetPoints( mPoints );
            mPolyData->SetLines( tLines );

            mMapper->SetInputData( mPolyData );

            mActor->SetMapper( mMapper );
            mActor->GetProperty()->SetLineWidth( aLineWidth );
            mActor->GetProperty()->SetLighting( false );
            mActor->GetProperty()->SetColor( 1.0, 1.0, 1.0 );
        }

//------------------------------------------------------------------------------

        void
        Curve::set_color( real aRed, real aGreen, real aBlue )
        {
            mActor->GetProperty()->SetColor( aRed, aGreen, aBlue );
        }
    }
}
