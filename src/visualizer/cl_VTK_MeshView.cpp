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
#include "vtktools.hpp"
#include "cl_VTK_MeshView.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Facet.hpp"


namespace belfem
{
    namespace vtk
    {
        BlockActor::BlockActor(
            Points                   aPoints,
            Map< id_t, vtkIdType > & aNodeMap,
            Map< id_t, vtkIdType > & aElementMap,
            mesh::Block *            aBlock ) :
            mBlock( aBlock ),
            mGrid( UnstructuredGrid::New() ),
            mMapper( DataSetMapper::New() ),
            mActor( Actor::New() )
        {
            // share the parent mesh's point set; this grid only owns cells
            mGrid->SetPoints( aPoints );

            Cell< mesh::Element * > & tElements = mBlock->elements();

            mGrid->Allocate( tElements.size() );

            Cell< vtkIdType > tVtkIDs ;
            Vector< id_t > tNodeIDs ;

            for ( mesh::Element * tElement : tElements )
            {
                uint n = tElement->number_of_nodes() ;
                tNodeIDs.set_size( n );
                tVtkIDs.set_size( n );
                get_node_ids( tElement, tNodeIDs );

                for ( uint k=0; k<n; ++k )
                {
                    tVtkIDs( k ) = aNodeMap( tNodeIDs( k ) );
                }

                aElementMap[ tElement->id() ]
                    = mGrid->InsertNextCell( vtk_type( tElement->type() ), n, tVtkIDs.data() );
            }

            mMapper->SetInputData( mGrid );
            mActor->SetMapper( mMapper );

            if ( mBlock->domain_type() == DomainType::Air )
            {
                mActor->SetVisibility( false );
            }
            else
            {
                mActor->SetVisibility( ! mBlock->is_hidden() );
            }

        }

        SideSetActor::SideSetActor(
            Points                   aPoints,
            Map< id_t, vtkIdType > & aNodeMap,
            Map< id_t, vtkIdType > & aElementMap,
            mesh::SideSet *          aSideSet ) :
            mSideSet( aSideSet ),
            mGrid( UnstructuredGrid::New() ),
            mMapper( DataSetMapper::New() ),
            mActor( Actor::New() )
        {
            // share the parent mesh's point set; this grid only owns cells
            mGrid->SetPoints( aPoints );

            Cell< mesh::Facet * > & tFacets = mSideSet->facets();

            mGrid->Allocate( tFacets.size() );

            Cell< vtkIdType > tVtkIDs ;
            Vector< id_t > tNodeIDs ;

            for ( mesh::Facet * tFacet : tFacets )
            {
                mesh::Element * tElement = tFacet->element();

                uint n = tElement->number_of_nodes() ;
                tNodeIDs.set_size( n );
                tVtkIDs.set_size( n );
                get_node_ids( tElement, tNodeIDs );

                for ( uint k=0; k<n; ++k )
                {
                    tVtkIDs( k ) = aNodeMap( tNodeIDs( k ) );
                }

                aElementMap[ tElement->id() ]
                    = mGrid->InsertNextCell( vtk_type( tElement->type() ), n, tVtkIDs.data() );
            }

            mMapper->SetInputData( mGrid );
            mActor->SetMapper( mMapper );

            mActor->SetVisibility( ! mSideSet->is_hidden() );
        }


//------------------------------------------------------------------------------

        MeshView::MeshView( belfem::Mesh * aMesh ) :
            mMesh( aMesh ),
            mPoints( Points::New() )
        {
            mPoints->Allocate( mMesh->number_of_nodes() );
            Cell< mesh::Node * > & tNodes = mMesh->nodes();
            for( mesh::Node * tNode : tNodes )
            {
                mNodeMap[ tNode->id() ]
                    = mPoints->InsertNextPoint(
                    tNode->x(),
                    tNode->y(),
                    tNode->z() );
            }

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();
            mBlocks.reserve( tBlocks.size() );
            for ( mesh::Block * tBlock : tBlocks )
            {
                mBlocks.push( new BlockActor( mPoints, mNodeMap, mElementMap, tBlock ) );
            }

            Cell< mesh::SideSet * > & tSideSets = mMesh->sidesets();
            mSideSets.reserve( tSideSets.size() );
            for ( mesh::SideSet * tSideSet : tSideSets )
            {
                mSideSets.push( new SideSetActor( mPoints, mNodeMap, mElementMap, tSideSet ) );
            }
        }

//------------------------------------------------------------------------------

        MeshView::~MeshView()
        {
            // delete the lightweight wrappers; the VTK objects they hold are
            // released by their own smart pointers
            for ( BlockActor * tBlock : mBlocks )
            {
                delete tBlock ;
            }
            for ( SideSetActor * tSideSet : mSideSets )
            {
                delete tSideSet ;
            }
        }

//------------------------------------------------------------------------------

        Cell< Actor >
        MeshView::actors()
        {
            Cell< Actor > aActors ;
            aActors.reserve( mBlocks.size() + mSideSets.size() );

            for ( BlockActor * tBlock : mBlocks )
            {
                aActors.push( tBlock->actor() );
            }
            for ( SideSetActor * tSideSet : mSideSets )
            {
                aActors.push( tSideSet->actor() );
            }

            return aActors ;
        }

//------------------------------------------------------------------------------

        void
        MeshView::set_user_transform( const Transform & aTransform )
        {
            for ( BlockActor * tBlock : mBlocks )
            {
                tBlock->actor()->SetUserTransform( aTransform );
            }
            for ( SideSetActor * tSideSet : mSideSets )
            {
                tSideSet->actor()->SetUserTransform( aTransform );
            }
        }
    }
}