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

#include "cl_Mesh_OrderConverter.hpp"
#include "fn_unique.hpp"
#include "meshtools.hpp"
#include "assert.hpp"
#include "fn_max.hpp"
#include "fn_sum.hpp"
#include "cl_Element_Factory.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    OrderConverter::OrderConverter( belfem::Mesh * aInMesh, belfem::Mesh * aOutMesh ) :
        mInMesh( aInMesh )
    {
        bool tMode3D = this->check_input_mesh();

        this->create_shared_facets();
        this->create_unique_facets();

        this->get_max_node_id();

        this->copy_nodes();

        if ( tMode3D )
        {
            this->create_edges();
            this->create_edge_nodes();
            this->create_facet_nodes_3d();
            this->create_center_nodes_3d();
        }
        else
        {
            this->create_facet_nodes_2d();
            this->create_center_nodes_2d();
        }

        if( aOutMesh == NULL )
        {
            mMesh = new Mesh( aInMesh->number_of_dimensions(), aInMesh->master() );
        }
        else
        {
            mMesh = aOutMesh;
        }

        this->create_nodes();
        this->create_elements();

        if ( tMode3D )
        {
            this->link_edge_nodes();
            this->link_facet_nodes();
            this->link_center_nodes_3d();
        }
        else
        {
            this->link_facet_nodes();
            this->link_center_nodes_2d();
        }

        this->create_blocks();
        this->create_sidesets();

        mMesh->finalize();
    }

//------------------------------------------------------------------------------

    Mesh *
    OrderConverter::mesh()
    {
        return mMesh;
    }

//------------------------------------------------------------------------------
    bool
    OrderConverter::check_input_mesh()
    {
        uint tNumBlocks = mInMesh->number_of_blocks();

        Vector< uint > tOrder( tNumBlocks, 0 );
        Vector< uint > tDimension( tNumBlocks, 0 );

        Cell< mesh::Block * > & tBlocks = mInMesh->blocks();

        for( uint b=0; b<tNumBlocks; ++b )
        {
            mesh::Element * tElement = tBlocks( b )->element( 0 );

            if( mesh::interpolation_order( tElement->type() )
                == InterpolationOrder::LINEAR )
            {
                tOrder( b ) = 1;
            }

            if( mesh::dimension( tElement->type() ) == 3 )
            {
                tDimension( b ) = 1;
            }
        }

        BELFEM_ERROR( sum( tOrder ) == tNumBlocks,
                "All input blocks must contain linear elements only" );

        uint tNumberOf3DBlocks = sum( tDimension );

        if( tNumberOf3DBlocks == tNumBlocks )
        {
            return true;
        }
        else if( tNumberOf3DBlocks == 0 )
        {
            return false;
        }
        else
        {
            BELFEM_ERROR( false, "can't upgrade a mesh with mixed 2D and 3D elements" );
            return false;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_shared_facets()
    {
        Cell<mesh::Element *> & tElements = mInMesh->elements();

        index_t tCount = 0;

        for ( mesh::Element * tElement : tElements )
        {
            tCount += tElement->number_of_elements();
        }

        Vector<index_t> tIDs( tCount, 0 );

        tCount = 0;

        index_t tMasterIndex;
        index_t tSlaveIndex;
        index_t tSwap;

        index_t tNumberOfElements = tElements.size();

        for ( mesh::Element * tElement : tElements )
        {
            for ( uint e = 0; e < tElement->number_of_elements(); ++e )
            {
                mesh::Element * tNeighbor = tElement->element( e );

                uint tDim = std::min( mesh::dimension( tElement->type() ),
                                      mesh::dimension( tNeighbor->type() ) );

                tMasterIndex = tElement->index();
                tSlaveIndex  = tNeighbor->index();

                if ( tSlaveIndex < tMasterIndex )
                {
                    tSwap = tMasterIndex;
                    tMasterIndex = tSlaveIndex;
                    tSlaveIndex = tSwap;
                }

                if ( this->count_common_nodes( tElement, tNeighbor ) >= tDim )
                {
                    tIDs( tCount++ ) = tSlaveIndex * tNumberOfElements + tMasterIndex;
                }
            }
        }

        unique( tIDs );

        tCount = 0;

        if( tIDs.length() > 0 )
        {
            mSharedFacetTable.set_size( 4, tIDs.length() - 1 );

            Cell<mesh::Node *> tFacetNodes;

            for ( uint k = 1; k < tIDs.length(); ++k )
            {
                tSlaveIndex = tIDs( k ) / tNumberOfElements;
                tMasterIndex = tIDs( k ) - tSlaveIndex * tNumberOfElements;

                this->get_facet_nodes( tElements( tMasterIndex ),
                                       tElements( tSlaveIndex ), tFacetNodes );

                mSharedFacetTable( 0, tCount ) = tMasterIndex;
                mSharedFacetTable( 1, tCount ) = tSlaveIndex;
                mSharedFacetTable( 2, tCount ) =
                        this->get_facet_index( tElements( tMasterIndex ), tFacetNodes );
                mSharedFacetTable( 3, tCount ) =
                        this->get_facet_index( tElements( tSlaveIndex ), tFacetNodes );

                ++tCount;
            }
        }

        mInMesh->unflag_all_nodes();
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_unique_facets()
    {
        Cell<mesh::Element *> & tElements = mInMesh->elements();

        Vector< uint > tNumFacetsPerElement( tElements.size(), 0 );

        index_t tNumSharedFacets = mSharedFacetTable.n_cols();

        index_t tCount = 0;
        for( mesh::Element * tElement : tElements )
        {
            tNumFacetsPerElement( tCount++ ) = tElement->number_of_facets();
        }
        index_t tMaxFacetsPerElement = max( tNumFacetsPerElement );

        Matrix< uint > tFacetsPerElement(
                tMaxFacetsPerElement,
                tElements.size(), tNumSharedFacets );

        for( index_t f=0; f<tNumSharedFacets; ++f )
        {
            tFacetsPerElement(
                    mSharedFacetTable( 2, f ),
                    mSharedFacetTable( 0, f ) ) = f;

            tFacetsPerElement(
                    mSharedFacetTable( 3, f ),
                    mSharedFacetTable( 1, f ) ) = f;
        }

        tCount = 0;
        index_t e = 0;
        for( mesh::Element * tElement : tElements )
        {
            for ( uint k = 0; k < tElement->number_of_facets(); ++k )
            {
                if ( tFacetsPerElement( k, e ) == tNumSharedFacets )
                {
                    ++tCount;
                }
            }

            ++e;
        }

        mUniqueFacetTable.set_size( 2, tCount );

        tCount = 0;
        e = 0;
        for( mesh::Element * tElement : tElements )
        {
            for ( uint k = 0; k < tElement->number_of_facets(); ++k )
            {
                if ( tFacetsPerElement( k, e ) == tNumSharedFacets )
                {
                    mUniqueFacetTable( 0, tCount ) = e;
                    mUniqueFacetTable( 1, tCount++ ) = k;
                }
            }

            ++e;
        }
    }

//------------------------------------------------------------------------------

    uint
    OrderConverter::count_common_nodes(
            mesh::Element * aElementA,
            mesh::Element * aElementB )
    {
        uint tNumNodesA = mesh::number_of_corner_nodes( aElementA->type() );
        uint tNumNodesB = mesh::number_of_corner_nodes( aElementB->type() );

        for( uint k=0; k<tNumNodesA; ++k )
        {
            aElementA->node( k )->unflag();
        }

        for( uint k=0; k<tNumNodesB; ++k )
        {
            aElementB->node( k )->flag();
        }

        uint aCount = 0;

        for( uint k=0; k<tNumNodesA; ++k )
        {
            if( aElementA->node( k )->is_flagged() )
            {
                ++aCount;
            }
        }

        return aCount;
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::get_facet_nodes(
            mesh::Element * aElementA,
            mesh::Element * aElementB,
            Cell< mesh::Node * > & aNodes )
    {
        uint tCount = this->count_common_nodes( aElementA, aElementB );

        uint tNumNodesA = mesh::number_of_corner_nodes( aElementA->type() );
        uint tNumNodesB = mesh::number_of_corner_nodes( aElementB->type() );

        if( tCount > 0 )
        {
            aNodes.set_size( tCount, nullptr );

            for( uint k=0; k<tNumNodesA; ++k )
            {
                aElementA->node( k )->unflag();
            }

            for( uint k=0; k<tNumNodesB; ++k )
            {
                aElementB->node( k )->flag();
            }

            tCount = 0;

            for( uint k=0; k<tNumNodesA; ++k )
            {
                if( aElementA->node( k )->is_flagged() )
                {
                    aNodes( tCount ++ ) = aElementA->node( k );
                }
            }
        }
        else
        {
            aNodes.clear();
        }
    }

//------------------------------------------------------------------------------

    uint
    OrderConverter::get_facet_index(
               mesh::Element * aElement,
            Cell< mesh::Node * > & aNodes )
    {

        uint tNumFacets = aElement->number_of_facets();
        Cell< mesh::Node * > tNodes;

        for( uint f=0; f<tNumFacets; ++f )
        {
            for( mesh::Node * tNode : aNodes )
            {
                tNode->unflag();
            }

            aElement->get_nodes_of_facet( f, tNodes );

            for( mesh::Node * tNode : tNodes )
            {
                tNode->flag();
            }

            uint tCount = 0;
            for( mesh::Node * tNode : aNodes )
            {
                if ( tNode->is_flagged() )
                {
                    ++tCount;
                }
            }

            if( tCount == aNodes.size() )
            {
                return f;
            }
        }

        BELFEM_ERROR( false, "Could not find facet index for element %lu",
                     ( long unsigned int ) aElement->id() );

        return BELFEM_UINT_MAX;
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_edges()
    {
        Cell< mesh::Element * > & tElements = mInMesh->elements();

        index_t tCount = 0;

        for( mesh::Element * tElement : tElements )
        {
            tCount += tElement->number_of_edges();
        }

        Vector< index_t > tIDs( tCount );

        index_t tNumNodes =  mInMesh->number_of_nodes();

        index_t tMasterIndex;
        index_t tSlaveIndex;
        index_t tSwap;

        Cell< mesh::Node * > tEdge;

        tCount = 0;
        for( mesh::Element * tElement : tElements )
        {
            for( uint e=0; e<tElement->number_of_edges(); ++e )
            {
                tElement->get_nodes_of_edge( e, tEdge );

                tMasterIndex = tEdge( 0 )->index();
                tSlaveIndex = tEdge( 1 )->index() ;

                if( tMasterIndex > tSlaveIndex )
                {
                    tSwap = tMasterIndex;
                    tMasterIndex = tSlaveIndex;
                    tSlaveIndex = tSwap;
                }

                tIDs( tCount++ ) = tSlaveIndex * tNumNodes + tMasterIndex;
            }
        }

        unique( tIDs );

        mEdges.set_size( 2, tIDs.length() );

        for( index_t k=0; k<tIDs.length(); ++k )
        {
            tSlaveIndex = tIDs( k ) / tNumNodes;
            tMasterIndex = tIDs( k ) - tSlaveIndex * tNumNodes;

            mEdges( 0, k ) = tMasterIndex;
            mEdges( 1, k ) = tSlaveIndex;

            mEdgeMap[ tIDs( k ) ] = k;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::get_max_node_id()
    {
        Cell< mesh::Node * > & tNodes = mInMesh->nodes();

        Vector< id_t > tIDs( tNodes.size() );

        index_t tCount = 0;

        for(  mesh::Node * tNode : tNodes )
        {
            tIDs( tCount++ ) = tNode->id();
        }

        mNodeID = max( tIDs ) + 1;
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::copy_nodes()
    {
        Cell< mesh::Node * > & tNodes = mInMesh->nodes();

        mOriginalNodes.set_size( tNodes.size(), nullptr );

        index_t tCount = 0;

        for(  mesh::Node * tNode : tNodes )
        {
            mOriginalNodes( tCount ) = new mesh::Node(
                    tNode->id(),
                    tNode->x(),
                    tNode->y(),
                    tNode->z() );

            mOriginalNodes( tCount )->set_index( tCount );

            ++tCount;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_edge_nodes()
    {
        Cell< mesh::Node * > & tNodes = mInMesh->nodes();

        index_t tNumberOfEdges = mEdges.n_cols();

        mEdgeNodes.set_size( tNumberOfEdges, nullptr );

        for( index_t e = 0; e<tNumberOfEdges; ++e )
        {
            mesh::Node * tMaster = tNodes( mEdges( 0, e ) );
            mesh::Node * tSlave  = tNodes( mEdges( 1, e ) );

            mEdgeNodes( e ) = new mesh::Node(
                    mNodeID++,
                    0.5 * ( tMaster->x() + tSlave->x() ),
                    0.5 * ( tMaster->y() + tSlave->y() ),
                    0.5 * ( tMaster->z() + tSlave->z() ) );
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_facet_nodes_2d()
    {
        Cell< mesh::Element * > & tElements = mInMesh->elements();

        Cell< mesh::Node * > tNodes;

        index_t tNumUniqueFacets = mUniqueFacetTable.n_cols();
        index_t tNumSharedFacets = mSharedFacetTable.n_cols();

        mFacetNodes.set_size( tNumUniqueFacets + tNumSharedFacets, nullptr );

        index_t tCount = 0;

        mFacesPerElement.set_size( 4, tElements.size(), tCount );

        // create nodes that sit on unique facets
        for( uint f=0; f<tNumUniqueFacets; ++f )
        {
            mesh::Element * tElement = tElements( mUniqueFacetTable( 0, f ) );

            tElement->get_nodes_of_facet( mUniqueFacetTable( 1, f ), tNodes );

            real tX = ( tNodes( 0 )->x() + tNodes( 1 )->x() ) * 0.5;
            real tY = ( tNodes( 0 )->y() + tNodes( 1 )->y() ) * 0.5;
            real tZ = ( tNodes( 0 )->z() + tNodes( 1 )->z() ) * 0.5;

            mFacesPerElement(
                    mUniqueFacetTable( 1, f ),
                    mUniqueFacetTable( 0, f ) ) = tCount;

            mFacetNodes( tCount ++ ) = new mesh::Node( mNodeID++, tX, tY, tZ );
        }

        // create nodes that sit on shared facets
        for( uint f=0; f<tNumSharedFacets; ++f )
        {
            mesh::Element * tElement = tElements( mSharedFacetTable( 0, f ) );

            tElement->get_nodes_of_facet( mSharedFacetTable( 2, f ), tNodes );

            real tX = ( tNodes( 0 )->x() + tNodes( 1 )->x() ) * 0.5;
            real tY = ( tNodes( 0 )->y() + tNodes( 1 )->y() ) * 0.5;
            real tZ = ( tNodes( 0 )->z() + tNodes( 1 )->z() ) * 0.5;

            mFacesPerElement( mSharedFacetTable( 2, f ),
                              mSharedFacetTable( 0, f ) ) = tCount;

            mFacesPerElement( mSharedFacetTable( 3, f ),
                             mSharedFacetTable( 1, f ) ) = tCount;

           mFacetNodes( tCount++ ) =
                new mesh::Node( mNodeID++, tX, tY, tZ );
        }

    }
//------------------------------------------------------------------------------

    void
    OrderConverter::create_facet_nodes_3d()
    {
        Cell< mesh::Element * > & tElements = mInMesh->elements();

        index_t tCount = 0;

        Cell< mesh::Node * > tNodes;

        index_t tNumUniqueFacets = mUniqueFacetTable.n_cols();
        for( uint f=0; f<tNumUniqueFacets; ++f )
        {
            mesh::Element * tElement = tElements( mUniqueFacetTable( 0, f ) );

            tElement->get_nodes_of_facet( mUniqueFacetTable( 1, f ), tNodes );

            if( tNodes.size() == 4 )
            {
                ++tCount;
            }
        }

        index_t tNumSharedFacets = mSharedFacetTable.n_cols();
        for( uint f=0; f<tNumSharedFacets; ++f )
        {
            mesh::Element * tElement = tElements( mSharedFacetTable( 0, f ) );

            tElement->get_nodes_of_facet( mSharedFacetTable( 2, f ), tNodes );

            if( tNodes.size() == 4 )
            {
                ++tCount;
            }
        }

        mFacesPerElement.set_size( 6, tElements.size(), tCount );
        mFacetNodes.set_size( tCount, nullptr );

        tCount = 0;
        for( uint f=0; f<tNumUniqueFacets; ++f )
        {
            mesh::Element * tElement = tElements( mUniqueFacetTable( 0, f ) );

            tElement->get_nodes_of_facet( mUniqueFacetTable( 1, f ), tNodes );

            if( tNodes.size() == 4 )
            {
                real tX = 0.0;
                real tY = 0.0;
                real tZ = 0.0;
                for ( uint k = 0; k < 4; ++k )
                {
                    tX += tNodes( k )->x();
                    tY += tNodes( k )->y();
                    tZ += tNodes( k )->z();
                }
                tX  *= 0.25;
                tY  *= 0.25;
                tZ  *= 0.25;

                mFacesPerElement( mUniqueFacetTable( 1, f ),
                        mUniqueFacetTable( 0, f ) ) = tCount;

                mFacetNodes( tCount++ ) =
                        new mesh::Node( mNodeID++, tX, tY, tZ );
            }
        }

        // create shared faces
        for( uint f=0; f<tNumSharedFacets; ++f )
        {
            mesh::Element * tElement = tElements( mSharedFacetTable( 0, f ) );

            tElement->get_nodes_of_facet( mSharedFacetTable( 2, f ), tNodes );

            if( tNodes.size() == 4 )
            {
                real tX = 0.0;
                real tY = 0.0;
                real tZ = 0.0;
                for ( uint k = 0; k < 4; ++k )
                {
                    tX += tNodes( k )->x();
                    tY += tNodes( k )->y();
                    tZ += tNodes( k )->z();
                }
                tX *= 0.25;
                tY *= 0.25;
                tZ *= 0.25;

                mFacesPerElement( mSharedFacetTable( 2, f ),
                                  mSharedFacetTable( 0, f ) ) = tCount;

                mFacesPerElement( mSharedFacetTable( 3, f ),
                                  mSharedFacetTable( 1, f ) ) = tCount;

                mFacetNodes( tCount++ ) =
                        new mesh::Node( mNodeID++, tX, tY, tZ );
            }
        }

    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_center_nodes_2d()
    {
        Cell< mesh::Element * > & tElements = mInMesh->elements();

        index_t tCount = 0;

        for ( mesh::Element * tElement : tElements )
        {
            if ( mesh::geometry_type( tElement->type() )
                == GeometryType::QUAD )
            {
                ++tCount;
            }
        }

        mCenterNodes.set_size( tCount, nullptr );

        tCount = 0;
        for ( mesh::Element * tElement : tElements )
        {
            if ( mesh::geometry_type( tElement->type() )
                 == GeometryType::QUAD )
            {
                real tX = 0.0;
                real tY = 0.0;
                real tZ = 0.0;

                for( uint k=0; k<4; ++k )
                {
                    tX += tElement->node( k )->x();
                    tY += tElement->node( k )->y();
                    tZ += tElement->node( k )->z();
                }

                tX *= 0.25;
                tY *= 0.25;
                tZ *= 0.25;

                mCenterNodes( tCount ++ ) = new mesh::Node( mNodeID++, tX, tY, tZ );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_center_nodes_3d()
    {
        Cell< mesh::Element * > & tElements = mInMesh->elements();

        index_t tCount = 0;

        for ( mesh::Element * tElement : tElements )
        {
            if (   tElement->type() == ElementType::HEX8 )

            {
                ++tCount;
            }
        }

        mCenterNodes.set_size( tCount, nullptr );

        tCount = 0;

        for ( mesh::Element * tElement : tElements )
        {
            if ( tElement->type() == ElementType::HEX8 )
            {
                real tX = 0.0;
                real tY = 0.0;
                real tZ = 0.0;

                for( uint k=0; k<8; ++k )
                {
                    tX += tElement->node( k )->x();
                    tY += tElement->node( k )->y();
                    tZ += tElement->node( k )->z();
                }

                tX *= 0.125;
                tY *= 0.125;
                tZ *= 0.125;

                mCenterNodes( tCount ++ ) = new mesh::Node( mNodeID++, tX, tY, tZ );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_nodes()
    {
        index_t tCount =
                  mOriginalNodes.size()
                + mEdgeNodes.size()
                + mFacetNodes.size()
                + mCenterNodes.size();

        Cell< mesh::Node * > & tNodes = mMesh->nodes();

        tNodes.set_size( tCount, nullptr );

        tCount = 0;

        for( mesh::Node * tNode : mOriginalNodes )
        {
            tNodes( tCount ++ ) = tNode;
        }
        for( mesh::Node * tNode : mEdgeNodes )
        {
            tNodes( tCount ++ ) = tNode;
        }
        for( mesh::Node * tNode : mFacetNodes )
        {
            tNodes( tCount ++ ) = tNode;
        }
        for( mesh::Node * tNode : mCenterNodes )
        {
            tNodes( tCount ++ ) = tNode;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_elements()
    {
        Cell< mesh::Node * >    & tNodes    = mMesh->nodes();
        Cell< mesh::Element * > & tElements = mMesh->elements();
        Cell< mesh::Element * > & tInElements = mInMesh->elements();

        mesh::ElementFactory tFactory;

        tElements.set_size( tInElements.size(), nullptr );

        index_t tCount = 0;

        for( mesh::Element * tInElement : tInElements )
        {
            mesh::Element * tElement = tFactory.create_element(
                    this->upgrade_type( tInElement->type() ),
                    tInElement->id() );

            tElement->set_geometry_tag( tInElement->geometry_tag() );
            tElement->set_physical_tag( tInElement->physical_tag() );
            tElement->set_owner( tInElement->owner() );

            for( uint k=0; k<tInElement->number_of_nodes(); ++k )
            {
                tElement->insert_node( tNodes( tInElement->node( k )->index() ), k );
            }

            tElements( tCount++ ) = tElement;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::link_edge_nodes()
    {
        Cell< mesh::Element *> & tElements = mMesh->elements();

        Cell<mesh::Node *> tEdge;

        index_t tMasterIndex;
        index_t tSlaveIndex;
        index_t tSwap;
        index_t tNumberOfNodes = mInMesh->number_of_nodes();

        for ( mesh::Element * tElement : tElements )
        {
            uint tOff = mesh::number_of_corner_nodes( tElement->type() );

            for( uint e=0; e<tElement->number_of_edges(); ++e )
            {
                tElement->get_nodes_of_edge( e, tEdge );

                tMasterIndex = tEdge( 0 )->index();
                tSlaveIndex = tEdge( 1 )->index();

                if ( tMasterIndex > tSlaveIndex )
                {
                    tSwap = tMasterIndex;
                    tMasterIndex = tSlaveIndex;
                    tSlaveIndex = tSwap;
                }

                mesh::Node * tNode = mEdgeNodes(
                        mEdgeMap( tSlaveIndex * tNumberOfNodes + tMasterIndex ));

                tElement->insert_node( tNode, e + tOff );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::link_facet_nodes()
    {
        Cell< mesh::Element* > & tElements = mMesh->elements();

        Vector< uint > tTable;

        index_t e = 0;
        for( mesh::Element* tElement : tElements )
        {
            this->facet_table( tElement->type(), tTable );

            for( uint f=0; f<tTable.length(); ++f )
            {
                tElement->insert_node(
                        mFacetNodes( mFacesPerElement( f, e ) ),
                        tTable( f ) );
            }

            ++e;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::facet_table( const ElementType & aType, Vector< uint > & aTable )
    {
        switch ( aType )
        {
            case ( ElementType::TRI6 ) :
            {
                aTable = { 3, 4, 5 };
                break;
            }
            case ( ElementType::QUAD9 ) :
            {
                aTable = { 4, 5, 6, 7 };
                break;
            }
            case ( ElementType::TET10 ) :
            {
                aTable = {};
                break;
            }
            case ( ElementType::PENTA18 ) :
            {
                aTable = { 15, 16, 17 };
                break;
            }
            case ( ElementType::HEX27 ) :
            {
                aTable = { 25, 24, 26, 23, 21, 22 };
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "invalid type " );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::link_center_nodes_2d()
    {
        Cell< mesh::Element* > & tElements = mMesh->elements();

        index_t tCount = 0;

        for( mesh::Element* tElement : tElements )
        {
            if ( tElement->type() == ElementType::QUAD9 )
            {
                tElement->insert_node( mCenterNodes( tCount++ ), 8 );
            }
        }

    }

//------------------------------------------------------------------------------

    void
    OrderConverter::link_center_nodes_3d()
    {
        Cell< mesh::Element* > & tElements = mMesh->elements();

        index_t tCount = 0;

        for( mesh::Element* tElement : tElements )
        {
            if ( tElement->type() == ElementType::HEX27 )
            {
                tElement->insert_node( mCenterNodes( tCount++ ), 20 );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_blocks()
    {
        Cell< mesh::Block * > & tInBlocks = mInMesh->blocks();
        Cell< mesh::Block * > & tBlocks = mMesh->blocks();
        Cell< mesh::Element * > & tElements = mMesh->elements();

        uint tNumBlocks = mInMesh->number_of_blocks();
        tBlocks.set_size( tNumBlocks, nullptr );

        uint tCount = 0;

        for( mesh::Block * tInBlock :  tInBlocks )
        {
            index_t tNumElements = tInBlock->number_of_elements();

            mesh::Block * tBlock = new mesh::Block(
                    tInBlock->id(),
                    tNumElements );

            tBlock->label() = tInBlock ->label();

            for( index_t e=0; e<tNumElements; ++e )
            {
                tBlock->insert_element( tElements( tInBlock->element( e )->index() ) );
            }

            tBlocks( tCount++ ) = tBlock;
        }
    }

//------------------------------------------------------------------------------

    void
    OrderConverter::create_sidesets()
    {
        Cell< mesh::Element * > & tElements = mMesh->elements();
        Cell< mesh::SideSet * > & tInSideSets = mInMesh->sidesets();
        Cell< mesh::SideSet * > & tSideSets = mMesh->sidesets();

        mesh::ElementFactory tFactory;

        index_t tCount = 0;

        if( tInSideSets.size() > 0 )
        {
            tSideSets.set_size( tInSideSets.size(), nullptr );

            for( mesh::SideSet * tInSideset : tInSideSets )
            {
                index_t tNumFacets = tInSideset->number_of_facets();

                mesh::SideSet * tSideSet = new mesh::SideSet(
                        tInSideset->id(),
                        tNumFacets );

                tSideSet->label() = tInSideset->label();

                for( uint f=0; f<tNumFacets; ++f )
                {
                    mesh::Facet * tInFacet = tInSideset->facet_by_index( f );

                    mesh::Element * tInElement = tInFacet->element();

                    mesh::Element * tElement = tFactory.create_element(
                            this->upgrade_type( tInElement->type() ),
                            tInElement->id() );

                    tElement->set_geometry_tag( tInElement->geometry_tag() );
                    tElement->set_physical_tag( tInElement->physical_tag() );
                    tElement->set_owner( tInElement->owner() );

                    mesh::Element * tMaster = tElements(
                            tInFacet->master()->index() );

                    Cell< mesh::Node * > tNodes;
                    tMaster->get_nodes_of_facet(
                            tInFacet->index_on_master(), tNodes );

                    for( uint k=0; k<tNodes.size(); ++k )
                    {
                        tElement->insert_node( tNodes( k ), k );
                    }

                    mesh::Facet * tFacet = new mesh::Facet( tElement );

                    tFacet->set_master( tMaster, tInFacet->index_on_master() );

                    if( tInFacet->has_slave() )
                    {
                        mesh::Element * tSlave = tElements( tInFacet->slave()->index() );
                        tFacet->set_slave( tSlave, tInFacet->index_on_slave() );
                    }

                    tSideSet->insert_facet( tFacet );
                }

                tSideSets( tCount++ ) = tSideSet;
            }
        }
    }

//------------------------------------------------------------------------------

    ElementType
    OrderConverter::upgrade_type( const ElementType & aType )
    {
        switch( aType )
        {
            case( ElementType::LINE2 ):
            {
                return ElementType::LINE3;
            }
            case( ElementType::TRI3 ):
            {
                return ElementType::TRI6;
            }
            case( ElementType::QUAD4 ) :
            {
                return ElementType::QUAD9;
            }
            case( ElementType::TET4 ) :
            {
                return ElementType::TET10;
            }
            case( ElementType::PENTA6 ) :
            {
                return ElementType::PENTA18;
            }
            case( ElementType::HEX8 ) :
            {
                return ElementType::HEX27;
            }
            default:
            {
                BELFEM_ERROR( false, "UNSUPPORTED ELEMENT TYPE" );
                return ElementType::UNDEFINED;
            }
        }
    }

//------------------------------------------------------------------------------
}
