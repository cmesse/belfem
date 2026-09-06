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

#include "cl_ProtoMesh.hpp"
#include "cl_Mesh.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Mesh_PeriodicityFactory.hpp"
namespace belfem
{
    namespace mesh
    {

//------------------------------------------------------------------------------

        ProtoMesh::ProtoMesh( Mesh * aMesh ) :
            mMesh( aMesh )
        {
            mMetaData     = new proto::MetaData ;
            mNodeData     = new proto::NodeData ;
            mElementData  = new proto::ElementData ;
            mEdgeData     = new proto::EdgeData ;
            mFaceData     = new proto::FaceData ;
            mFacetData    = new proto::FacetData ;

            mElementExtra = new proto::ElementExtra ;
            mFacetExtra   = new proto::FacetExtra ;

            mVertexData = new proto::ElementData ;

            mControlPointData = new proto::ControlPointData ;

            mPeriodicityData = new proto::PeriodicityData ;

            mTMatrixData = new proto::TMatrixData ;
        }

//------------------------------------------------------------------------------

        ProtoMesh::~ProtoMesh()
        {
            if ( mMetaData != nullptr )
            {
                delete mMetaData ;
            }
            if ( mNodeData != nullptr )
            {
                delete mNodeData ;
            }
            if ( mElementData != nullptr )
            {
                delete mElementData ;
            }
            if ( mEdgeData != nullptr )
            {
                delete mEdgeData ;
            }
            if ( mFaceData != nullptr )
            {
                delete mFaceData ;
            }
            if ( mFacetData != nullptr )
            {
                delete mFacetData ;
            }
            if ( mElementExtra != nullptr )
            {
                delete mElementExtra ;
            }
            if ( mFacetExtra != nullptr )
            {
                delete mFacetExtra ;
            }
            if ( mVertexData != nullptr )
            {
                delete mVertexData ;
            }
            if ( mControlPointData != nullptr )
            {
                delete mControlPointData ;
            }
            if ( mPeriodicityData != nullptr )
            {
                delete mPeriodicityData ;
            }
            if ( mTMatrixData != nullptr )
            {
                delete mTMatrixData ;
            }
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::populate_meta_data()
        {
            proto::MetaData * tMetaData = meta_data();

            tMetaData->mNumberOfDimensions = mMesh->number_of_dimensions();
            Cell< index_t > & tNumEntities = tMetaData->mNumberOfEntities;

            tNumEntities.set_size( static_cast< uint >( EntityType::UNDEFINED ) );
            tNumEntities( static_cast< uint >( EntityType::NODE ) ) = mMesh->number_of_nodes();
            tNumEntities( static_cast< uint >( EntityType::EDGE ) ) = mMesh->number_of_edges();
            tNumEntities( static_cast< uint >( EntityType::FACE ) ) = mMesh->number_of_faces();
            tNumEntities( static_cast< uint >( EntityType::FACET ) ) = mMesh->number_of_facets();
            tNumEntities( static_cast< uint >( EntityType::CONTROLPOINT ) ) = mMesh->number_of_control_points();

            Cell< index_t > & tNumGroups = tMetaData->mNumberOfGroups;

            tNumGroups.set_size( 4 );
            tNumGroups( 0 ) = mMesh->number_of_blocks()  ;
            tNumGroups( 1 ) = mMesh->number_of_sidesets() ;
            tNumGroups( 2 ) = mMesh->curves().size() ;
            tNumGroups( 3 ) = mMesh->thin_shells().size() ;

        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::populate_node_data()
        {
            Cell< Node * > & tNodes = mMesh->nodes();

            uint tNumDim = mMesh->number_of_dimensions();
            index_t tNumNodes = tNodes.size();

            index_t tCount = 0 ;

            // populate node IDs
            Cell< id_t > & tIDs = mNodeData->mIDs ;
            tIDs.set_size( tNumNodes );
            for ( Node * tNode : tNodes )
            {
                tIDs( tCount++ ) = tNode->id();
            }

            // populate node coordinates
            Matrix< real > & tCoords = mNodeData->mCoords ;
            tCoords.set_size( tNumDim, tNumNodes );
            tCount = 0 ;
            if ( tNumDim == 2 )
            {
                for ( Node * tNode : tNodes )
                {
                    tCoords( 0, tCount ) = tNode->x();
                    tCoords( 1, tCount ) = tNode->y();
                    ++tCount ;
                }
            }
            else
            {
                for ( Node * tNode : tNodes )
                {
                    tCoords( 0, tCount ) = tNode->x();
                    tCoords( 1, tCount ) = tNode->y();
                    tCoords( 2, tCount ) = tNode->z();
                    ++tCount ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::create_nodes()
        {
            Cell< Node * > & tNodes = mMesh->nodes();
            BELFEM_ASSERT( tNodes.size() == 0, "Nodes are already allocated" );

            index_t tNumNodes = mNodeData->mIDs.size();

            // allocate container
            tNodes.set_size( tNumNodes, nullptr );

            // create the node entities
            for ( index_t k=0; k<tNumNodes; ++k )
            {
                // create node
                Node * tNode = new Node( mNodeData->mIDs( k ) );

                // set node index
                tNode->set_index( k );

                // add node to map
                mNodeMap[ tNode->id() ] = tNode ;

                // add node to container
                tNodes( k ) = tNode;
            }

            if ( mNodeData->mOwners.size() > 0 )
            {
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    // set node owner
                    tNodes( k )->set_owner( mNodeData->mOwners( k ) );
                }
            }

            // set node coordinates
            if ( mMesh->number_of_dimensions() == 2 )
            {
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    tNodes( k )->set_coords(
                        mNodeData->mCoords( 0, k ),
                        mNodeData->mCoords( 1, k ) );
                }
            }
            else
            {
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    tNodes( k )->set_coords(
                        mNodeData->mCoords( 0, k ),
                        mNodeData->mCoords( 1, k ),
                        mNodeData->mCoords( 2, k ) );
                }
            }

            // assign node duplicates
            index_t tCount = 0 ;

            if ( mNodeData->mDuplicateData.size() > 0 )
            {
                // get the number of originals
                index_t tNumOrignals = mNodeData->mDuplicateData( tCount++ );

                // assign duplicate and original data
                for ( index_t i=0; i<tNumOrignals; ++i )
                {
                    // get the original node
                    Node * tOrg = mNodeMap( mNodeData->mDuplicateData( tCount++ ) ) ;

                    // get the number of duplicates
                    uint tNumDuplicates = mNodeData->mDuplicateData( tCount++ );

                    // allocate duplicate container
                    tOrg->allocate_duplicate_container( tNumDuplicates );

                    for ( uint d=0; d<tNumDuplicates; ++d )
                    {
                        // get the duplicate node
                        Node * tDup = mNodeMap( mNodeData->mDuplicateData( tCount++ ) ) ;
                        tOrg->add_duplicate( tDup );
                        tDup->set_original( tOrg );
                    }
                }
            }

            // delete the node data
            delete mNodeData ;
            mNodeData = nullptr ;
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::reset_node_data()
        {
            // delete the node data
            delete mNodeData ;
            mNodeData = nullptr ;
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::populate_element_data( const bool aPopulateTopology, const bool aPopulateGeo )
        {

            index_t tNumElements = mMesh->elements().size() ;

            Cell< id_t >   & tIDs    = mElementData->mIDs ;
            tIDs.set_size( tNumElements );
            index_t tCount = 0 ;
            bool tHavePhys = false ;

            for ( Element * tElement : mMesh->elements() )
            {
                if ( tElement->physical_tag() != 0 )
                {
                    tHavePhys = true ;
                    break ;
                }
            }

            for ( Block * tBlock : mMesh->blocks() )
            {

                Cell< Element * > & tElements = tBlock->elements();
                for ( Element * tElement : tElements )
                {
                    tIDs( tCount++ )  = tElement->id();
                }
            }

            BELFEM_ASSERT( tCount == tNumElements, "Number of elements does not match (is %lu, expect %lu)",
                ( long unsigned int ) tCount, ( long unsigned int ) tNumElements );

            if ( aPopulateGeo )
            {
                Cell< uint >   & tGeo    = mElementData->mGeometryTags ; //<- block IDs
                tGeo.set_size( tNumElements );

                tCount = 0 ;
                for ( Block * tBlock : mMesh->blocks() )
                {
                    Cell< Element * > & tElements = tBlock->elements();

                    for ( Element * tElement : tElements )
                    {
                        tGeo( tCount++ )  = tElement->geometry_tag();
                    }
                }
            }

            if ( tHavePhys )
            {
                Cell< uint >   & tPhys   = mElementData->mPhysicalTags ;
                tPhys.set_size( tNumElements );

                tCount = 0 ;
                for ( Block * tBlock : mMesh->blocks() )
                {
                    Cell< Element * > & tElements = tBlock->elements();

                    for ( Element * tElement : tElements )
                    {
                        tPhys( tCount++ ) = tElement->physical_tag();
                    }
                }
            }


            if ( aPopulateTopology )
            {
                Cell< id_t >   & tTopo   = mElementData->mTopology ;

                tCount = 0 ;
                for ( Block * tBlock : mMesh->blocks() )
                {
                    tCount += tBlock->number_of_elements() * number_of_nodes( tBlock->element_type() );
                }

                mElementData->mTopology.set_size( tCount );

                tCount = 0 ;
                for ( Block * tBlock : mMesh->blocks() )
                {
                    Cell< Element * > & tElements = tBlock->elements();
                    index_t n = number_of_nodes( tBlock->element_type() );

                    for ( Element * tElement : tElements )
                    {
                        for ( index_t k = 0; k < n; ++k )
                        {
                            tTopo( tCount++ ) = tElement->node( k )->id();
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::populate_facet_data( const bool aPopulateElements, const bool aPopulateTopology, const bool aPopulateGeo )
        {
            bool tHavePhys = false;

            // we need to count the number of facets
            // over the sidesets, because mMesh->facets() also contains thin shell facets.

            index_t tNumFacets = 0 ;
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                Cell< Facet * > & tFacets = tSideSet->facets();
                if ( ! tHavePhys )
                {
                    for ( Facet * tFacet : tFacets )
                    {
                        if ( tFacet->physical_tag() != 0 )
                        {
                            tHavePhys = true;
                            break ;
                        }
                    }
                }
                tNumFacets += tFacets.size();
            }

            Cell< id_t >   & tIDs    = mFacetData->mIDs ;
            tIDs.set_size( tNumFacets );



            for ( Facet * tFacet : mMesh->facets() )
            {
                if ( tFacet->physical_tag() != 0 )
                {
                    tHavePhys = true;
                    break ;
                }

            }


            index_t tCount = 0 ;
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {

                Cell< Facet * > & tFacets = tSideSet->facets();
                for ( Facet * tFacet : tFacets )
                {
                    tIDs( tCount++ )  = tFacet->id();
                }
            }

            BELFEM_ASSERT( tCount == tNumFacets, "Number of facets does not match (is %lu, expect %lu)",
                ( long unsigned int ) tCount, ( long unsigned int ) tNumFacets );

            if ( aPopulateGeo )
            {
                Cell< uint >   & tGeo    = mFacetData->mGeometryTags ; //<- sideset IDs
                tGeo.set_size( tNumFacets );

                tCount = 0 ;
                for ( SideSet * tSideSet : mMesh->sidesets() )
                {
                    Cell< Facet * > & tFacets = tSideSet->facets();

                    for ( Facet * tFacet : tFacets )
                    {
                        tGeo( tCount++ )  = tFacet->element()->geometry_tag();
                    }
                }
            }

            if ( tHavePhys )
            {
                Cell< uint >   & tPhys   = mFacetData->mPhysicalTags ;
                tPhys.set_size( tNumFacets );

                tCount = 0 ;
                for ( SideSet * tSideSet : mMesh->sidesets() )
                {
                    Cell< Facet * > & tFacets = tSideSet->facets();

                    for ( Facet * tFacet : tFacets )
                    {
                        tPhys( tCount++ ) = tFacet->element()->physical_tag();
                    }
                }
            }

            if ( aPopulateTopology )
            {
                Cell< id_t >   & tTopo   = mFacetData->mTopology ;

                tCount = 0 ;
                for ( SideSet * tSideSet : mMesh->sidesets() )
                {
                    tCount += tSideSet->number_of_facets() * number_of_nodes( tSideSet->element_type() );
                }

                mFacetData->mTopology.set_size( tCount );

                tCount = 0 ;
                for ( SideSet * tSideSet : mMesh->sidesets() )
                {
                    Cell< Facet * > & tFacets = tSideSet->facets();


                    index_t n = number_of_nodes( tSideSet->element_type() );

                    for ( Facet * tFacet : tFacets )
                    {
                        for ( index_t k = 0; k < n; ++k )
                        {
                            tTopo( tCount++ ) = tFacet->node( k )->id();
                        }
                    }
                }
            }

            if ( aPopulateElements )
            {
                Cell< id_t > & tMasterIDs = mFacetData->mMasterIDs ;
                Cell< id_t > & tSlaveIDs = mFacetData->mSlaveIDs ;
                Cell< uchar > & tMasterIndices = mFacetData->mIndicesOnMaster ;
                Cell< uchar > & tSlaveIndices = mFacetData->mIndicesOnSlave ;
                Cell< uchar > & tOrientations = mFacetData->mOrientationsOnSlave ;

                tMasterIDs.set_size( tNumFacets, gNoID );
                tSlaveIDs.set_size( tNumFacets, gNoID );
                tMasterIndices.set_size( tNumFacets, BELFEM_UCHAR_MAX );
                tSlaveIndices.set_size( tNumFacets, BELFEM_UCHAR_MAX );
                tOrientations.set_size( tNumFacets, BELFEM_UCHAR_MAX );

                tCount = 0 ;

                for ( SideSet * tSideSet : mMesh->sidesets() )
                {
                    Cell< Facet * > & tFacets = tSideSet->facets();

                    for ( Facet * tFacet : tFacets )
                    {
                        if ( tFacet->has_master() )
                        {
                            tMasterIDs( tCount ) = tFacet->master()->id() ;
                            tMasterIndices( tCount ) = tFacet->index_on_master() ;
                        }
                        if ( tFacet->has_slave() )
                        {
                            tSlaveIDs( tCount ) = tFacet->slave()->id() ;
                            tSlaveIndices( tCount ) = tFacet->index_on_slave() ;
                            tOrientations( tCount ) = tFacet->orientation_on_slave() ;
                        }
                        ++tCount ;
                    }
                }
            }
        }

        void
        ProtoMesh::populate_edge_data( const bool aPopulateTopology )
        {
            if ( ! mMesh->edges_exist() ) return ;

            Cell< Edge * > & tEdges = mMesh->edges();

            // count memory
            index_t tCount = tEdges.size();
            Cell< id_t > & tIDs = mEdgeData->mIDs;
            tIDs.set_size( tCount, 0 );
            index_t e = 0 ;
            for ( Edge * tEdge : tEdges )
            {
                tIDs( e++ ) = tEdge->id();
                tCount += tEdge->number_of_nodes();
            }

            if ( ! aPopulateTopology ) return ;

            Cell< id_t > & tTopology = mEdgeData->mTopology ;
            tTopology.set_size( tCount, 0 );
            tCount = 0 ;

            for ( Edge * tEdge : tEdges )
            {
                uint n = tEdge->number_of_nodes();
                tTopology( tCount++ ) = n ;
                for ( uint k=0; k<n; ++k )
                {
                    tTopology( tCount++ ) = tEdge->node( k )->id();
                }
            }

        }

        void
        ProtoMesh::populate_face_data( const bool aPopulateElements )
        {
            if ( ! mMesh->faces_exist() ) return ;

            index_t tNumFaces = mMesh->faces().size();

            Cell< Face * > & tFaces = mMesh->faces();
            Cell< id_t >   & tIDs    = mFaceData->mIDs ;
            tIDs.set_size( tNumFaces );
            index_t tCount = 0 ;
            for ( Face * tFace : tFaces )
            {
                tIDs( tCount++ ) = tFace->id();
            }
            if ( aPopulateElements )
            {
                Cell< id_t > & tMasterIDs = mFaceData->mMasterIDs ;
                Cell< id_t > & tSlaveIDs = mFaceData->mSlaveIDs ;
                Cell< uchar > & tMasterIndices = mFaceData->mIndicesOnMaster ;
                Cell< uchar > & tSlaveIndices = mFaceData->mIndicesOnSlave ;
                Cell< uchar > & tOrientations = mFaceData->mOrientationsOnSlave ;

                tMasterIDs.set_size( tNumFaces, gNoID );
                tSlaveIDs.set_size( tNumFaces, gNoID );
                tMasterIndices.set_size( tNumFaces, BELFEM_UCHAR_MAX );
                tSlaveIndices.set_size( tNumFaces, BELFEM_UCHAR_MAX );
                tOrientations.set_size( tNumFaces, BELFEM_UCHAR_MAX );

                tCount = 0 ;

                for ( Face * tFace : tFaces )
                {
                    if ( tFace->master() != nullptr )
                    {
                        tMasterIDs( tCount ) = tFace->master()->id();
                        tMasterIndices( tCount ) = tFace->index_on_master();
                    }
                    if ( tFace->slave() != nullptr )
                    {
                        tSlaveIDs( tCount ) = tFace->slave()->id();
                        tSlaveIndices( tCount ) = tFace->index_on_slave();
                        tOrientations( tCount ) = tFace->orientation_on_slave();
                    }
                    tCount++;
                }
            }
        }

        void
        ProtoMesh::populate_control_point_data( const bool aPopulateTopology )
        {
            index_t tNumPoints = mMesh->number_of_control_points();

            if ( tNumPoints == 0 ) return;

            Cell< id_t > & tIDs = mControlPointData->mIDs ;
            tIDs.set_size( tNumPoints );

            uint tNumDim = mMesh->number_of_dimensions();
            Matrix< real > & tCoords = mControlPointData->mCoords ;
            tCoords.set_size( tNumDim, tNumPoints  );

            Cell< ControlPoint * > & tPoints = mMesh->control_points();

            index_t tCount = 0 ;
            for ( ControlPoint * tPoint : tPoints )
            {
                for ( uint d=0; d<tNumDim; ++d )
                {
                    tCoords( d , tCount ) = tPoint->x( d );
                }
                tIDs( tCount++ ) = tPoint->id() ;
            }

            if ( ! aPopulateTopology ) return ;

            tCount = 0 ;
            for ( Block * tBlock : mMesh->blocks() )
            {
                Cell< Element * > & tElements = tBlock->elements() ;

                tCount += 2 * tElements.size() ;
                for ( Element * tElement : tElements )
                {
                    tCount += tElement->number_of_control_points();
                }
            }

            Cell< id_t > & tTopo = mControlPointData->mElementTopology ;
            tTopo.set_size( tCount );
            tCount = 0 ;

            for ( Block * tBlock : mMesh->blocks() )
            {
                Cell< Element * > & tElements = tBlock->elements() ;
                for ( Element * tElement : tElements )
                {
                    tTopo( tCount++ ) = tElement->id() ;
                    uint n = tElement->number_of_control_points() ;
                    tTopo( tCount++ ) = n ;
                    for ( uint k=0 ; k<n; ++k )
                    {
                        tTopo( tCount++ ) = tElement->control_point( k )->id() ;
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::reset_element_data()
        {
            // delete the node data
            delete mElementData ;
            mElementData = nullptr ;
        }

        void
        ProtoMesh::reset_facet_data()
        {
            delete mFacetData ;
            mFacetData = nullptr ;
        }

        void
        ProtoMesh::reset_edge_data()
        {
            delete mEdgeData ;
            mEdgeData = nullptr ;
        }

        void
        ProtoMesh::reset_face_data()
        {
            delete mFaceData ;
            mFaceData = nullptr ;
        }

        void
        ProtoMesh::reset_control_point_data()
        {
            delete mControlPointData ;
            mControlPointData = nullptr ;
        }

        void ProtoMesh::reset_periodicity_data()
        {
            delete mPeriodicityData ;
            mPeriodicityData = nullptr ;
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::create_elements()
        {
            index_t tNumElements = mElementData->mIDs.size();

            Cell< Element * > & tElements = mMesh->elements();
            tElements.set_size( tNumElements, nullptr );

            ElementFactory tFactory;

            Cell< id_t >   & tIDs    = mElementData->mIDs ;
            Cell< proc_t > & tOwners = mElementData->mOwners ;
            Cell< uchar >  & tTypes  = mElementData->mTypes ;
            Cell< uint >   & tGeo    = mElementData->mGeometryTags ; //<- block IDs
            Cell< uint >   & tPhys   = mElementData->mPhysicalTags ;
            Cell< id_t >   & tTopo   = mElementData->mTopology ;


            if ( tTypes.size() > 0 ) // used by mesh distributor
            {
                for ( index_t e=0; e<tNumElements; ++e )
                {

                    // create a new element
                    tElements( e ) = tFactory.create_element(
                        static_cast< ElementType >( tTypes( e ) ), tIDs( e ) );
                }
            }
            else // used by BFM file
            {
                Map< id_t, ElementType > tTypeMap ;

                for ( proto::GroupData & tBlock : mBlockData )
                {
                    tTypeMap[ tBlock.mID ] = tBlock.mElementType ;
                }

                for ( index_t e=0; e<tNumElements; ++e )
                {
                    // create a new element
                    tElements( e ) = tFactory.create_element( tTypeMap( tGeo( e ) ), tIDs( e ) );
                }
            }

            index_t tCount = 0 ;
            for ( index_t e=0; e<tNumElements; ++e )
            {
                Element * tElement = tElements( e ) ;

                // set the local index
                tElement->set_index( e );

                // set block id
                tElement->set_geometry_tag( tGeo( e ) );

                // connect element with nodes
                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    tElement->insert_node( mNodeMap( tTopo( tCount++ ) ), k );
                }

                // add element to map
                mElementMap[ tElement->id() ] = tElement ;
            }

            if ( tOwners.size() > 0 )
            {
                tCount = 0 ;
                for ( Element * tElement : tElements )
                {
                    tElement->set_owner( tOwners( tCount++ ) );
                }
            }

            if ( tPhys.size() > 0 )
            {
                tCount = 0 ;
                for ( Element * tElement : tElements )
                {
                    tElement->set_physical_tag( tPhys( tCount++ ) );
                }
            }

            // clear data container
            delete mElementData ;
            mElementData = nullptr ;
        }

//------------------------------------------------------------------------------
        
        void
        ProtoMesh::create_edges()
        {
            Cell< id_t >   & tIDs    = mEdgeData->mIDs ;
            Cell< proc_t > & tOwners = mEdgeData->mOwners ;
            Cell< id_t >   & tTopo   = mEdgeData->mTopology ;

            index_t tNumEdges = tIDs.size();

            if ( tNumEdges == 0 ) return ;

            Cell< Edge * > & tEdges = mMesh->edges();
            BELFEM_ASSERT( tEdges.size() == 0, "Edges are already allocated" );

            // allocate container
            tEdges.set_size( tNumEdges, nullptr );

            index_t tCount = 0 ;

            for ( index_t e=0; e<tNumEdges; ++e )
            {
                Edge * tEdge = new Edge ;

                // set id and index
                tEdge->set_id( tIDs( e ) );
                tEdge->set_index( e );

                // connect nodes
                uint n = tTopo( tCount++ );
                tEdge->allocate_node_container( n );
                for ( uint k=0; k<n; ++k )
                {
                    tEdge->insert_node( mNodeMap( tTopo( tCount++ ) ) , k );
                }

                // add edge to map
                mEdgeMap[ tEdge->id() ] = tEdge ;

                // add edge to container
                tEdges( e ) = tEdge;
            }

            if ( tOwners.size() > 0 )
            {
                tCount = 0;
                for ( Edge * tEdge : tEdges )
                {
                    tEdge->set_owner( tOwners( tCount++ ) );
                }
            }

            // clear data container
            delete mEdgeData ;
            mEdgeData = nullptr ;
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::create_faces()
        {
            Cell< id_t >   & tIDs     = mFaceData->mIDs ;
            Cell< proc_t > & tOwners  = mFaceData->mOwners ;
            Cell< id_t >   & tMIDs    = mFaceData->mMasterIDs ;
            Cell< id_t >   & tSIDs    = mFaceData->mSlaveIDs ;
            Cell< uchar >  & tMIdx    = mFaceData->mIndicesOnMaster ;
            Cell< uchar >  & tSIdx    = mFaceData->mIndicesOnSlave ;
            Cell< uchar >  & tSorn    = mFaceData->mOrientationsOnSlave ;

            index_t tNumFaces = tIDs.size();
            if ( tNumFaces == 0 ) return ;

            Cell< Face * > & tFaces = mMesh->faces();
            BELFEM_ASSERT( tFaces.size() == 0, "Faces are already allocated" );

            tFaces.set_size( tNumFaces, nullptr );

            if ( mMesh->number_of_dimensions() == 2 )
            {
                for ( index_t f=0; f<tNumFaces; ++f )
                {
                    Face * tFace = new Face( mElementMap( tMIDs( f ) ) );

                    tFace->set_id( tIDs( f ) );
                    tFace->set_index( f );

                    mFaceMap[ tFace->id() ] = tFace ;
                    tFaces( f ) = tFace;
                }
            }
            else
            {
                 for ( index_t f=0; f<tNumFaces; ++f )
                 {
                     Element * tMaster = tMIDs( f ) == gNoID ? nullptr : mElementMap( tMIDs( f ) );
                     Element * tSlave = tSIDs( f ) == gNoID ? nullptr : mElementMap( tSIDs( f ) );

                     Face * tFace = new Face( tMaster, tMIdx( f ), tSlave, tSIdx( f ), tSorn( f ) );

                     tFace->set_id( tIDs( f ) );
                     tFace->set_index( f );

                     mFaceMap[ tFace->id() ] = tFace ;
                     tFaces( f ) = tFace;
                 }
            }

            if ( tOwners.size() > 0 )
            {
                index_t tCount = 0;
                for ( Face * tFace : tFaces )
                {
                    tFace->set_owner( tOwners( tCount++ ) );
                }
            }
            delete mFaceData ;
            mFaceData = nullptr ;
        }
        
//------------------------------------------------------------------------------

        void
        ProtoMesh::create_facets()
        {
             Cell< id_t >   & tIDs     = mFacetData->mIDs ;
             Cell< proc_t > & tOwners  = mFacetData->mOwners ;

             Cell< id_t >   & tMIDs    = mFacetData->mMasterIDs ;
             Cell< id_t >   & tSIDs    = mFacetData->mSlaveIDs ;

             Cell< uchar >   & tMIdx    = mFacetData->mIndicesOnMaster ;
             Cell< uchar >   & tSIdx    = mFacetData->mIndicesOnSlave ;
             Cell< uchar >   & tSors    = mFacetData->mOrientationsOnSlave ;

             Cell< uchar >  & tTypes   = mFacetData->mTypes ;
             Cell< uint >   & tGeo     = mFacetData->mGeometryTags ;
             Cell< uint >   & tPhys    = mFacetData->mPhysicalTags ;

             Cell< id_t >   & tTopo    = mFacetData->mTopology ;

             ElementFactory tFactory;

             index_t tNumFacets = tIDs.size();
             Cell< Facet * > & tFacets = mMesh->facets();
             BELFEM_ASSERT( tFacets.size() == 0, "Facets are already allocated" );

             tFacets.set_size( tNumFacets, nullptr );

             index_t tCount = 0 ;

             if ( tTopo.size() > 0 )
             {
                 for ( index_t f=0; f<tNumFacets; ++f )
                 {
                     Element * tElement = tFactory.create_element( static_cast< ElementType >( tTypes( f ) ), tIDs( f ) );

                     tElement->set_index( f );

                     tElement->set_geometry_tag( tGeo( f ) );

                     // link element with nodes
                     for ( index_t k=0; k<tElement->number_of_nodes(); ++k )
                     {
                         tElement->insert_node( mNodeMap( tTopo( tCount++ ) ), k );
                     }

                     // create facet
                     Facet * tFacet = new Facet( tElement );

                     // set master
                     if ( tMIDs( f ) != gNoID )
                     {
                         tFacet->set_master( mElementMap( tMIDs( f ) ), tMIdx( f ), false );
                     }

                     // set slave
                     if ( tSIDs( f ) != gNoID )
                     {
                         tFacet->set_slave( mElementMap( tSIDs( f ) ), tSIdx( f ),tSors( f ) );
                     }

                     // add facet to map
                     mFacetMap[ tFacet->id() ] = tFacet ;

                     // add facet to container
                     tFacets( f ) = tFacet;
                 }
             }
             else
             {
                 for ( index_t f=0; f<tNumFacets; ++f )
                 {
                     Element * tElement = tFactory.create_element( static_cast< ElementType >( tTypes( f ) ), tIDs( f ) );

                     tElement->set_index( f );
                     tElement->set_geometry_tag( tGeo( f ) );

                     Facet * tFacet = new Facet( tElement );

                     // set master
                     BELFEM_ASSERT( tMIDs( f ) != gNoID, "Facet %lu has neither master nor topology assigned",
                          ( long unsigned int ) tElement->id() );

                     tFacet->set_master( mElementMap( tMIDs( f ) ), tMIdx( f ), true );

                     // set slave
                     if ( tSIDs( f ) != gNoID )
                     {
                         tFacet->set_slave( mElementMap( tSIDs( f ) ), tSIdx( f ),tSors( f ) );
                     }

                     // add facet to map
                     mFacetMap[ tFacet->id() ] = tFacet ;

                     // add facet to container
                     tFacets( f ) = tFacet;
                 }
             }

            if ( tPhys.size() > 0 )
            {
                for ( index_t f=0; f<tNumFacets; ++f )
                {
                    tFacets( f )->set_physical_tag( tPhys( f ) );
                }
            }

            if ( tOwners.size() > 0 )
            {
                for ( index_t f=0; f<tNumFacets; ++f )
                {
                    tFacets( f )->set_owner( tOwners( f ) );
                }
            }
            delete mFacetData ;
            mFacetData = nullptr ;
        }
                
//------------------------------------------------------------------------------

        void
        ProtoMesh::create_vertices()
        {
            ElementFactory tFactory ;
            Cell< Element * > & tVertices = mMesh->vertices();

            index_t tNumVertices = mVertexData->mIDs.size();
            tVertices.set_size( tNumVertices, nullptr );

            for ( index_t v = 0; v<tNumVertices; ++v )
            {
                Element * tVertex = tFactory.create_element( ElementType::VERTEX, mVertexData->mIDs( v ) );

                tVertex->insert_node( mNodeMap( mVertexData->mTopology( v ) ), 0 );

                tVertex->set_index( v );
                tVertex->set_owner( mVertexData->mOwners( v ) );
                tVertex->set_geometry_tag( mVertexData->mGeometryTags( v ) );
                tVertex->set_physical_tag( mVertexData->mPhysicalTags( v ) );

                tVertices( v ) = tVertex;
            }

            delete mVertexData ;
            mVertexData = nullptr;
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::create_control_points()
        {
            index_t tNumControlPoints = mControlPointData->mIDs.size();

            if ( tNumControlPoints == 0 ) return;

            Cell< ControlPoint * > & tControlPoints = mMesh->control_points();
            Matrix< real > & tCoords = mControlPointData->mCoords;

            tControlPoints.set_size( tNumControlPoints, nullptr );

            for ( index_t p=0; p<tNumControlPoints; ++p )
            {
                    tControlPoints( p ) = new ControlPoint(
                        mControlPointData->mIDs( p ),
                        tCoords( 0, p ),
                        tCoords( 1, p ),
                        tCoords( 2, p ));

                mControlPointMap[ mControlPointData->mIDs( p ) ] = tControlPoints( p );
            }

            Cell< id_t > & tTopo = mControlPointData->mElementTopology ;
            index_t tSize = tTopo.size();
            index_t tCount = 0 ;

            while ( tCount < tSize )
            {
                Element * tElement = mElementMap( tTopo( tCount++ ) );
                uint tNumPoints = tTopo( tCount++ );
                if ( tNumPoints > 0 )
                {
                    tElement->allocate_control_points_container( tNumPoints );
                    for ( uint k=0; k<tNumPoints; ++k )
                    {
                        tElement->insert_control_point( mControlPointMap( tTopo( tCount++ ) ) , k );
                    }
                }
            }

            delete mControlPointData ;
            mControlPointData = nullptr ;
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::create_element_extra()
        {
            Cell< id_t > & tNData  = mElementExtra->mNeighborData ;
            Cell< id_t > & tEData  = mElementExtra->mEdgeData ;
            Cell< id_t > & tFData  = mElementExtra->mFaceData ;
            Cell< uchar > & tDData = mElementExtra->mEdgeDirections ;

            index_t tCount = 0 ;

            if ( tNData.size() > 0 )
            {
                index_t tNumElems = tNData( tCount++ );

                for ( index_t i=0; i<tNumElems; ++i )
                {
                    Element * tElement = mElementMap( tNData( tCount++ ) );


                    uint n = tNData( tCount++ ) ;
                    tElement->allocate_element_container( n );

                    for ( uint k=0; k<n; ++k )
                    {
                        tElement->insert_element( mElementMap( tNData( tCount++ ) ) );
                    }
                }
                mMesh->set_connectivity( Connectivity::ElementToElement );
                tCount = 0 ;
            }

            if ( tEData.size() > 0 )
            {
                index_t tNumElems = tEData( tCount++ );
                index_t tDCount = 0 ;

                for ( index_t i=0; i<tNumElems; ++i )
                {
                    Element * tElement = mElementMap( tEData( tCount++ ) );

                    uint n = tElement->number_of_edges();
                    tElement->allocate_edge_container();

                    for ( uint k=0; k<n; ++k )
                    {
                        tElement->insert_edge( mEdgeMap( tEData( tCount++ ) ) , k );
                        tElement->set_edge_direction( k, tDData( tDCount++ ) == 1 );
                    }
                }

                mMesh->set_connectivity( Connectivity::ElementToEdge );
                tCount = 0 ;
            }

            // curved elements, if they exist
            if ( mElementExtra->mCurvedElementIDs.size() > 0 )
            {
                for ( id_t tID : mElementExtra->mCurvedElementIDs )
                {
                    mElementMap( tID )->set_curved_flag();
                }
            }
            if ( tFData.size() > 0 )
            {
                index_t tNumElems = tFData( tCount++ );
                for ( index_t i=0; i<tNumElems; ++i )
                {
                    Element * tElement = mElementMap( tFData( tCount++ ) );
                    uint n = tElement->number_of_faces();
                    tElement->allocate_face_container();

                    for ( uint k=0; k<n; ++k )
                    {
                        tElement->insert_face( mFaceMap( tFData( tCount++ ) ) , k );
                    }
                }

                mMesh->set_connectivity( Connectivity::ElementToFace );
            }

            delete mElementExtra ;
            mElementExtra = nullptr ;
        }
        
//------------------------------------------------------------------------------
        
        void
        ProtoMesh::create_facet_extra()
        {
            Cell< id_t > & tNData = mFacetExtra->mNeighborData ;

            if ( tNData.size() == 0 ) return;

            index_t tCount = 0 ;
            index_t tNumFacets = tNData( tCount++ );

            for ( index_t f=0; f<tNumFacets; ++f )
            {
                Facet * tFacet = mFacetMap( tNData( tCount++ ) );
                uint n = tNData( tCount++ );
                tFacet->set_facet_counter( n );
                tFacet->allocate_facet_container();

                for ( uint i=0; i<n; ++i )
                {
                    tFacet->add_facet( mFacetMap( tNData( tCount++ ) ) );
                }
            }

            // curved facets, if they exist
            for ( id_t tID : mFacetExtra->mCurvedFacetIDs )
            {
                mFacetMap( tID )->element()->set_curved_flag();
            }

            delete mFacetExtra ;
            mFacetExtra = nullptr ;
        }
        
//------------------------------------------------------------------------------

        void
        ProtoMesh::populate_block_data()
        {
            mBlockData.set_size( mMesh->number_of_blocks() );

            index_t tCount = 0 ;
            for ( Block * tBlock : mMesh->blocks() )
            {
                proto::GroupData   & tData = mBlockData( tCount++ ) ;
                tData.mID          = tBlock->id();
                tData.mLabel       = tBlock->label();
                tData.mNumElements = tBlock->number_of_elements();
                tData.mElementType = tBlock->element_type();
                tData.mDomainType  = tBlock->domain_type();
                tData.mHasEdges    = tBlock->has_edges();
                tData.mHasFaces    = tBlock->has_faces();
                tData.mHidden      = tBlock->is_hidden();
                tData.mThickness   = tBlock->thickness();
            }

        }

        void
        ProtoMesh::populate_sideset_data()
        {
            mSideSetData.set_size( mMesh->number_of_sidesets() );

            index_t tCount = 0 ;
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                proto::GroupData & tData = mSideSetData( tCount++ ) ;
                tData.mID          = tSideSet->id();
                tData.mLabel       = tSideSet->label();
                tData.mNumElements = tSideSet->number_of_facets();
                tData.mElementType = tSideSet->element_type();
                tData.mDomainType  = tSideSet->domain_type();
                tData.mHidden      = tSideSet->is_hidden();
            }
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::edge()
        {
             // get the data size
            uint tNumBlocks = mBlockData.size();

            // create a temporary index map
            Map< id_t, index_t > tIndex ;
            for ( uint b=0; b<tNumBlocks; ++b )
            {
                tIndex[ mBlockData( b ).mID ] = b ;
            }

            // allocate the element counter
            Cell< index_t > tCount( tNumBlocks, 0 );

            // count elements per block
            Cell< Element * > & tElements = mMesh->elements();

            // note: using this counter seems overkill when the elements
            // are already sorted by block.
            for ( Element * tElement : tElements )
            {
                ++tCount( tIndex( tElement->block_id() ) ) ;
            }

            // create a temporary block container
            Cell< Block * > tBlocks( tNumBlocks, nullptr );

            // create the blocks
            index_t tElemCount = 0 ;

            for ( uint b=0; b<tNumBlocks; ++b )
            {
                proto::GroupData & tData = mBlockData( b );

                // create a new block
                Block * tBlock = new Block( tData.mID, tCount( b ) );

                if ( tData.mLabel.size() > 0 ) tBlock->label() = tData.mLabel;

                // set the domain type
                tBlock->set_domain_type(  tData.mDomainType );

                // physical thickness ( NaN for volume blocks — same
                // semantics as a freshly created serial block )
                tBlock->set_thickness( tData.mThickness );

                // add block to container
                tBlocks( b ) = tBlock ;

                index_t n = tCount( b );

                for ( index_t e=0; e<n; ++e )
                {
                    BELFEM_ASSERT(  tElements( tElemCount )->block_id() == tBlock->id(),
                        "Corrupted element order in Protomesh" );

                    tBlock->insert_element( tElements( tElemCount++ ) );
                }
            }

            // check for element flags
            for ( uint b=0; b<tNumBlocks; ++b )
            {
                // get the block
                Block * tBlock = tBlocks( b );

                if ( tBlock->number_of_elements() > 0 )
                {
                    // get the first element
                    Element * tElement = tBlock->elements().first();

                    if ( tElement->has_edges() )
                    {
                        tBlock->set_edges_flag();
                    }
                    if ( tElement->has_faces() )
                    {
                        tBlock->set_faces_flag();
                    }

                    // add block to map
                    mBlockMap[ tBlock->id() ] = tBlock ;

                    mMesh->blocks().push( tBlock );
                }
                else
                {
                    delete tBlock ;
                }
            }
        }
        
//------------------------------------------------------------------------------

        void
        ProtoMesh::create_sidesets( const bool aKeepEmpty )
        {
            uint tNumSideSets = mSideSetData.size();
            Cell< SideSet * > tSideSets ;
            tSideSets.set_size( tNumSideSets, nullptr );

            // create a temporary index map
            Map< id_t, index_t > tIndex ;
            for ( uint s=0; s<tNumSideSets; ++s )
            {
                tIndex[ mSideSetData( s ).mID ] = s ;
            }

            // allocate the facet counter
            Cell< index_t > tCount( tNumSideSets, 0 );

            // count elements per block
            Cell< Facet * > & tFacets = mMesh->facets();

            for ( Facet * tFacet : tFacets )
            {
                ++tCount( tIndex( tFacet->sideset_id() ) ) ;
            }

            index_t tFacetCount = 0;

            for ( uint s=0; s<tNumSideSets; ++s )
            {
                proto::GroupData & tData = mSideSetData( s );
                SideSet * tSideSet = new SideSet( tData.mID, tCount( s ) );

                if ( tData.mLabel.size() > 0 ) tSideSet->label() = tData.mLabel ;
                tSideSet->hide( tData.mHidden );

                tSideSet->set_domain_type(  static_cast< DomainType >( tData.mDomainType ) );

                tSideSet->set_index( s );

                tSideSets( s ) = tSideSet;

                index_t n = tCount( s );

                for ( index_t f=0; f<n; ++f )
                {
                    BELFEM_ASSERT( tFacets( tFacetCount )->sideset_id() == tSideSet->id(),
                        "Corrupted facet order in Protomesh" );
                    tSideSet->insert_facet( tFacets( tFacetCount++ ) );
                }
            }

            for ( SideSet * tSideSet : tSideSets )
            {
                if ( tSideSet->number_of_facets() > 0 || aKeepEmpty )
                {
                    mMesh->sidesets().push( tSideSet );

                    // add sideset to mesh
                    mSideSetMap[ tSideSet->id() ] = tSideSet ;

                }
                else
                {
                    delete tSideSet ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::create_thinshells()
        {
            Cell< ThinShell * > & tThinShells = mMesh->thin_shells();

            uint m = mThinShellData.size();

            if ( m==0 ) return ;

            // count thin shells that are visible to this proc
            uint tCount = 0 ;
            for ( uint k=0; k<m; ++k )
            {
                if ( mSideSetMap.key_exists( mThinShellData( k ).mSideSetID ) )
                {
                    uint c=0 ;
                    for ( id_t b : mThinShellData( k ).mBlocksIDs )
                    {
                        if ( mBlockMap.key_exists( b ) ) ++c ;
                    }
                    BELFEM_ASSERT( c == 0 || c == mThinShellData( k ).mBlocksIDs.size(), "Thin shell blocks are not correctly defined" );

                    if ( c > 0 ) ++tCount ;
                }
            }
            tThinShells.reserve( tCount );

            for ( uint k=0; k<m; ++k )
            {
                if ( ! mSideSetMap.key_exists( mThinShellData( k ).mSideSetID ) ) continue ;

                uint c=0 ;
                for ( id_t b : mThinShellData( k ).mBlocksIDs )
                {
                    if ( mBlockMap.key_exists( b ) ) ++c ;
                }
                if ( c == 0 ) continue;

                uint n = mThinShellData( k ).mBlocksIDs.size() ;

                // create the shell object
                SideSet * tGhost = mSideSetMap.key_exists( mThinShellData( k ).mGhostSideSetID ) ? mSideSetMap( mThinShellData( k ).mGhostSideSetID ) : nullptr;
                ThinShell * tShell = new ThinShell( mSideSetMap( mThinShellData( k ).mSideSetID ) , tGhost );

                // set the blocks
                Cell< Block * > & tBlocks = tShell->blocks();
                tBlocks.set_size( n, nullptr );
                for ( uint b=0; b<n; ++b )
                {
                    if ( mBlockMap.key_exists( mThinShellData( k ).mBlocksIDs( b ) ) )
                    {
                        tBlocks( b ) = mBlockMap( mThinShellData( k ).mBlocksIDs( b ) );
                    }
                }

                // convert thicknesses to vector and set

                Vector< real > tThicknesses( n );
                for ( uint b=0; b<n; ++b )
                {
                    tThicknesses( b ) = mThinShellData( k ).mThicknesses( b );
                }
                tShell->set_thicknesses( tThicknesses );

                // set the material labels
                tShell->set_materials( mThinShellData( k ).mMaterials );

                // add shell to map
                mThinShellMap[ tShell->id() ] = tShell ;

                // add shell to container
                tThinShells.push( tShell );
            }
        }

        void
        ProtoMesh::create_periodicitiy( const bool aCrosslinkEntities )
        {
            if ( mPeriodicityData->mMasterPlane.size() == 0 && mPeriodicityData->mSlavePlane.size() == 0 ) return;

            PeriodicityFactory tFactory( mMesh, this );

            mMesh->set_periodicity( tFactory.from_proto( mPeriodicityData , aCrosslinkEntities ) );
        }

        void
        ProtoMesh::create_t_matrices()
        {
            index_t tNumNodes = mTMatrixData->mNumTargets( static_cast< uint >( EntityType::NODE ) );
            index_t tNumEdges = mTMatrixData->mNumTargets( static_cast< uint >( EntityType::EDGE ) );
            index_t tNumFaces = mTMatrixData->mNumTargets( static_cast< uint >( EntityType::FACE ) );
            index_t tNumElements = mTMatrixData->mNumTargets( static_cast< uint >( EntityType::ELEMENT ) );
            index_t tNumFacets = mTMatrixData->mNumTargets( static_cast< uint >( EntityType::FACET ) );
            index_t tNumControlPoints = mTMatrixData->mNumTargets( static_cast< uint >( EntityType::CONTROLPOINT ) );

            index_t tPivot = 0 ;
            index_t tCount = 0 ;
            const Cell< id_t >  & tTargets    = mTMatrixData->mTargetIDs ;
            const Cell< uint >  & tNumSources = mTMatrixData->mCounters ;
            const Cell< id_t >  & tSourceIDs  = mTMatrixData->mSourceIDs ;
            const Cell< uchar > & tTypes      = mTMatrixData->mTypes ;
            const Cell< real >  & tCoeffs     = mTMatrixData->mWeights ;
            Cell< Basis * > tSources ;
            Vector< real >  tWeights ;

            for ( index_t k=0; k<tNumNodes; ++k )
            {
                Node * tNode = mMesh->node( tTargets( tPivot ) );
                uint s = tNumSources( tPivot++ );

                tSources.set_size( s, nullptr );
                tWeights.set_size( s, 0 );

                for ( uint i=0; i<s; ++i )
                {
                    tSources( i ) = mMesh->basis( static_cast< EntityType >( tTypes( tCount ) ), tSourceIDs( tCount ) );
                    tWeights( i ) = tCoeffs( tCount++ );
                }
                tNode->set_sources( tSources, tWeights );
            }
            for ( index_t e=0; e<tNumEdges; ++e )
            {
                Edge * tEdge = mMesh->edge( tTargets( tPivot ) );
                uint s = tNumSources( tPivot++ );
                tSources.set_size( s, nullptr );
                tWeights.set_size( s, 0 );
                for ( uint i=0; i<s; ++i )
                {
                    tSources( i ) = mMesh->basis( static_cast< EntityType >( tTypes( tCount ) ), tSourceIDs( tCount ) );
                    tWeights( i ) = tCoeffs( tCount++ );
                }
                tEdge->set_sources( tSources, tWeights );
            }
            for ( index_t f=0; f<tNumFaces; ++f )
            {
                Face * tFace = mMesh->face( tTargets( tPivot ) );
                uint s = tNumSources( tPivot++ );
                tSources.set_size( s, nullptr );
                tWeights.set_size( s, 0 );
                for ( uint i=0; i<s; ++i )
                {
                    tSources( i ) = mMesh->basis( static_cast< EntityType >( tTypes( tCount ) ), tSourceIDs( tCount ) );
                    tWeights( i ) = tCoeffs( tCount++ );
                }
                tFace->set_sources( tSources, tWeights );
            }

            for ( index_t e=0; e<tNumElements; ++e )
            {
                Element * tElement = mMesh->element( tTargets( tPivot ) );
                uint s = tNumSources( tPivot++ );
                tSources.set_size( s, nullptr );
                tWeights.set_size( s, 0 );
                for ( uint i=0; i<s; ++i )
                {
                    tSources( i ) = mMesh->basis( static_cast< EntityType >( tTypes( tCount ) ), tSourceIDs( tCount ) );
                    tWeights( i ) = tCoeffs( tCount++ );
                }
                tElement->set_sources( tSources, tWeights );
            }

            for ( index_t f=0; f<tNumFacets; ++f )
            {
                Facet * tFacet = mMesh->facet( tTargets( tPivot ) );
                uint s = tNumSources( tPivot++ );
                tSources.set_size( s, nullptr );
                tWeights.set_size( s, 0 );
                for ( uint i=0; i<s; ++i )
                {
                    tSources( i ) = mMesh->basis( static_cast< EntityType >( tTypes( tCount ) ), tSourceIDs( tCount ) );
                    tWeights( i ) = tCoeffs( tCount++ );
                }
                tFacet->set_sources( tSources, tWeights );
            }

            for ( index_t c=0; c<tNumControlPoints; ++c )
            {
                ControlPoint * tControlPoint = mMesh->control_point( tTargets( tPivot ) );
                uint s = tNumSources( tPivot++ );
                tSources.set_size( s, nullptr );
                tWeights.set_size( s, 0 );
                for ( uint i=0; i<s; ++i )
                {
                    tSources( i ) = mMesh->basis( static_cast< EntityType >( tTypes( tCount ) ), tSourceIDs( tCount ) );
                    tWeights( i ) = tCoeffs( tCount++ );
                }
                tControlPoint->set_sources( tSources, tWeights );
            }

            mMesh->collect_hanging_basis();

        }

        void
        ProtoMesh::reconstruct_edge_connectivity()
        {
            if ( ! mMesh->edges_exist() ) return ;

            Cell< Edge * > & tEdges = mMesh->edges() ;
            key_t N = mMesh->number_of_nodes() ;

            // node-pair key of an edge. NOT unique on cut-enriched meshes:
            // twin edges along the cut boundary curves share both end nodes
            // when neither node was duplicated
            auto tEdgeKey = [ N ]( Edge * aEdge ) -> key_t
            {
                key_t A = aEdge->node( 0 )->index();
                key_t B = aEdge->node( 1 )->index();
                return std::min( A, B ) * N + std::max( A, B );
            };

            // edge pointers sorted by key plus parallel key array. The sort
            // must be stable so that edges sharing a key keep container order
            Cell< Edge * > tSortedEdges( tEdges );
            std::stable_sort( tSortedEdges.begin(), tSortedEdges.end(),
                [ &tEdgeKey ]( Edge * aA, Edge * aB )
                { return tEdgeKey( aA ) < tEdgeKey( aB ); } );

            index_t tNumEdges = tSortedEdges.size() ;
            Cell< key_t > tKeys( tNumEdges, 0 );
            for ( index_t k=0; k<tNumEdges; ++k )
            {
                tKeys( k ) = tEdgeKey( tSortedEdges( k ) );
            }

            // returns the edge for a node-pair key. For twin edges the LAST
            // one in container order wins, matching the last-write-wins
            // behavior of the map this replaces. Primary edges are appended
            // to the mesh before their material-interface duplicates, so
            // last = duplicate sheet
            auto tFindEdge = [ &tKeys, &tSortedEdges ]( key_t aKey ) -> Edge *
            {
                index_t tIndex = std::upper_bound( tKeys.begin(), tKeys.end(), aKey )
                    - tKeys.begin() ;

                BELFEM_ERROR( tIndex > 0 && tKeys( tIndex - 1 ) == aKey,
                    "no edge found for node pair key (corrupted file?)" );

                return tSortedEdges( tIndex - 1 );
            };

            // returns the FIRST edge for a node-pair key: with primaries
            // appended before duplicates, first = primary sheet
            auto tFindFirstEdge = [ &tKeys, &tSortedEdges ]( key_t aKey ) -> Edge *
            {
                index_t tIndex = std::lower_bound( tKeys.begin(), tKeys.end(), aKey )
                    - tKeys.begin() ;

                BELFEM_ERROR( tIndex < tKeys.size() && tKeys( tIndex ) == aKey,
                    "no edge found for node pair key (corrupted file?)" );

                return tSortedEdges( tIndex );
            };

            Cell< Node * > tNodes ;

            for ( proto::GroupData & tData : mBlockData )
            {
                if ( tData.mHasEdges )
                {
                    Block * tBlock = mBlockMap( tData.mID );

                    Cell< Element * > & tElements = tBlock->elements();
                    uint n = number_of_edges( tBlock->element_type() );

                    // side connector walls mix the edge sheets by construction:
                    // the bottom slots 0 and 1 take the duplicate sheet when it
                    // exists ( = last ), the top slots 2 and 3 always the
                    // primary sheet ( = first ), see
                    // ThinShellFactory::create_side_connectors. The node-pair
                    // key cannot tell the twins apart, so the slot decides
                    const bool tIsWall = tBlock->element_type() == ElementType::HEX8TB ;

                    for ( Element * tElement : tElements )
                    {
                        tElement->allocate_edge_container();
                        for ( uint e=0; e<n; ++e )
                        {
                            tElement->get_nodes_of_edge( e, tNodes );

                            key_t A = tNodes( 0 )->index();
                            key_t B = tNodes( 1 )->index();

                            key_t tKey = std::min( A, B ) * N + std::max( A, B );

                            tElement->insert_edge( tIsWall && e >= 2 ?
                                tFindFirstEdge( tKey ) : tFindEdge( tKey ), e );
                        }
                    }
                }
            }

#if !defined( NDEBUG ) || defined( DEBUG )

            // sanity check
            for ( proto::GroupData & tData : mBlockData )
            {
                if ( tData.mHasEdges )
                {
                    Block * tBlock = mBlockMap( tData.mID );

                    Cell< Element * > & tElements = tBlock->elements();

                    uint n = number_of_edges( tBlock->element_type() );

                    for ( Element * tElement : tElements )
                    {
                        for ( uint e=0; e<n; ++e )
                        {
                            BELFEM_ERROR( tElement->edge( e ) != nullptr,
                                "Element %lu has no edge %u", ( long unsigned int ) tElement->id(), ( unsigned int ) e );
                        }
                    }
                }
            }
#endif

            mMesh->set_connectivity( Connectivity::ElementToEdge );
        }

//------------------------------------------------------------------------------

        void
        ProtoMesh::normalize_edge_hangs()
        {
            // must run AFTER the hanging entities are loaded ( the sources
            // do not exist before load_hanging_entities() )
            if ( ! mMesh->edges_exist() ) return ;

            Cell< Edge * > & tEdges = mMesh->edges() ;
            key_t N = mMesh->number_of_nodes() ;

            // same node-pair key as reconstruct_edge_connectivity
            auto tEdgeKey = [ N ]( Edge * aEdge ) -> key_t
            {
                key_t A = aEdge->node( 0 )->index();
                key_t B = aEdge->node( 1 )->index();
                return std::min( A, B ) * N + std::max( A, B );
            };

            Cell< Edge * > tSortedEdges( tEdges );
            std::stable_sort( tSortedEdges.begin(), tSortedEdges.end(),
                [ &tEdgeKey ]( Edge * aA, Edge * aB )
                { return tEdgeKey( aA ) < tEdgeKey( aB ); } );

            index_t tNumEdges = tSortedEdges.size() ;
            Cell< key_t > tKeys( tNumEdges, 0 );
            for ( index_t k=0; k<tNumEdges; ++k )
            {
                tKeys( k ) = tEdgeKey( tSortedEdges( k ) );
            }

            // LAST twin for a key = the sheet the slot reconstruction hands
            // to the elements ( cf. reconstruct_edge_connectivity )
            auto tFindEdge = [ &tKeys, &tSortedEdges ]( key_t aKey ) -> Edge *
            {
                index_t tIndex = std::upper_bound( tKeys.begin(), tKeys.end(), aKey )
                    - tKeys.begin() ;

                BELFEM_ERROR( tIndex > 0 && tKeys( tIndex - 1 ) == aKey,
                    "no edge found for node pair key (corrupted file?)" );

                return tSortedEdges( tIndex - 1 );
            };

            // normalize 1:1 edge-on-edge hangs onto the tie representative:
            // the slot reconstruction hands tied element slots to the LAST
            // twin, so a hang that was created against the other sheet would
            // reference an edge no element flags for dofs — the side
            // connector inner copies are the first such consumers ( the dof
            // machinery would abort with "invalid number of dofs on source
            // edge : 0" ). At a tied station the surviving twin carries the
            // trace, so the hang follows it. Layer twins are same-oriented
            // by construction, so the weight transfers unchanged; sources
            // with a unique key resolve to themselves ( no-op )
            for ( Edge * tEdge : tEdges )
            {
                if ( ! tEdge->is_hanging() ) continue ;
                if ( tEdge->number_of_sources() != 1 ) continue ;
                if ( tEdge->source( 0 )->entity_type() != EntityType::EDGE ) continue ;

                Edge * tSource = reinterpret_cast< Edge * >( tEdge->source( 0 ) );
                Edge * tSurvivor = tFindEdge( tEdgeKey( tSource ) );

                if ( tSurvivor != tSource )
                {
                    const real tWeight = tEdge->weight( 0 );
                    tEdge->reset_source_container();
                    tEdge->allocate_source_container( 1 );
                    tEdge->add_source( tSurvivor, tWeight );
                }
            }
        }

        void
        ProtoMesh::reconstruct_face_connectivity()
        {
            if ( ! mMesh->faces_exist() ) return ;

            for ( proto::GroupData & tData : mBlockData )
            {
                if ( tData.mHasFaces )
                {
                    Block * tBlock = mBlockMap( tData.mID );

                    Cell< Element * > & tElements = tBlock->elements();

                    for ( Element * tElement : tElements )
                    {
                        tElement->allocate_face_container();
                    }
                }
            }

            Cell< Face * > & tFaces = mMesh->faces();

            for ( Face * tFace : tFaces )
            {
                if ( tFace->master() != nullptr )
                {
                    tFace->master()->insert_face( tFace, tFace->index_on_master() );
                }

                if ( tFace->slave() != nullptr )
                {
                    tFace->slave()->insert_face( tFace, tFace->index_on_slave() );
                }
            }

#if !defined( NDEBUG ) || defined( DEBUG )

            // sanity check
            for ( proto::GroupData & tData : mBlockData )
            {
                if ( tData.mHasFaces )
                {
                    Block * tBlock = mBlockMap( tData.mID );

                    Cell< Element * > & tElements = tBlock->elements();

                    uint n = number_of_faces( tBlock->element_type() );

                    for ( Element * tElement : tElements )
                    {
                        for ( uint f=0; f<n; ++f )
                        {
                            BELFEM_ERROR( tElement->face( f ) != nullptr,
                                "Element %lu has no face %u", ( long unsigned int ) tElement->id(), ( unsigned int ) f );
                        }
                    }
                }
            }
#endif

            mMesh->set_connectivity( Connectivity::ElementToFace );
        }

    }
}