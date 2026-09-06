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

#include "commtools.hpp"
#include "cl_ThinShellFactory.hpp"

#include "assert.hpp"
#include "cl_Logger.hpp"
#include "stringtools.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_unique.hpp"
#include "cl_FaceFactory.hpp"
#include "cl_Element_Factory.hpp"
#include "op_Graph_Vertex_ID.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace mesh
    {
        ThinShellFactory::ThinShellFactory(
        Mesh * aMesh,
            Cell< Node * > & aMasterNodes,
            Cell< Node * > & aSlaveNodes,
            Map< string, Material * > * aMaterialMap,
            const bool aCreateGhostFacets ):
            mCommRank( comm_rank() ),
            mNumDimensions( aMesh->number_of_dimensions() ),
            mMesh( aMesh ),
            mMaterialMap( aMaterialMap ),
            mMaxGroupID( aMesh->max_block_and_sideset_id() ),
            mMaxElementID( aMesh->max_element_id() ),
            mMaxNodeID( this->max_node_id() ),
            mMasterNodes( aMasterNodes ),
            mSlaveNodes( aSlaveNodes )
        {
            // assigned here rather than in the initializer list so the
            // member's position in the header cannot reorder the list
            mCreateGhostFacets = aCreateGhostFacets ;

            for ( Block * tBlock : aMesh->blocks() )
            {
                tBlock->unflag_elements();
                if ( tBlock->domain_type() == DomainType::Air || tBlock->domain_type() == DomainType::Ferro )
                {
                    Cell < Element * > & tElements = tBlock->elements();
                    for ( Element * tElement : tElements )
                    {
                        tElement->flag( 1 );
                    }
                }
                else
                {
                    tBlock->unflag_elements( 1 );
                }
            }
        }

        ThinShell *
        ThinShellFactory::create( Protoshell * aProtoShell )
        {

            this->reset_node_indices();

            // get the sidesets that are connected to this shell
            Cell< SideSet * > tSideSets;
            this->collect_sidesets( aProtoShell->sidesets(), tSideSets );

            // these are all facets that form a sideset, also determines the element type

            SideSet * tSideSet = new SideSet( ++mMaxGroupID, 0 );
            tSideSet->set_domain_type( DomainType::GeometryOnly );

            Cell< Facet * > & tFacets = tSideSet->facets();

            ElementType tElementType = this->collect_facets( tSideSets, tFacets );

            // the nodes that we'll need to use for the extrusion
            Cell< Node * > tNodes ;
            this->collect_nodes( tFacets, tNodes );

            Matrix< real > tNodeNormals;
            switch ( tElementType )
            {
                case ElementType::LINE2 :
                {
                    this->process_nodes_line2( tFacets, tNodes, tNodeNormals );
                    break ;
                }
                case ElementType::LINE3 :
                {
                    this->process_nodes_line3( tFacets, tNodes, tNodeNormals );
                    break ;
                }
                case ElementType::TRI3 :
                {
                    this->process_nodes_tri3( tFacets, tNodes, tNodeNormals );
                    break ;
                }
                case ElementType::TRI6 :
                {
                    this->process_nodes_tri6( tFacets, tNodes, tNodeNormals );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Element type not supported in %s", "ThinShellFactory::create" );
                }
            }

            // compute the offsets of each layer to the center
            Vector< real > tDistances ;

            uint tOrder = interpolation_order_numeric( tElementType );
            this->compute_distances( tOrder, aProtoShell->thicknesses(), tDistances );

            // create the temporary layer objects
            uint tNumLayers = tDistances.length() ;
            Cell< Layer * > tLayers( tNumLayers, nullptr );

            for ( uint l=0; l<tNumLayers; ++l )
            {
                tLayers( l ) = new Layer();
            }
            uint tCount = tOrder ;

            // a layer names its material by the section label of the input
            // file. The lookups below are raw find() calls: an unknown label
            // dereferences the end iterator and takes the run down in setup
            const Cell< string > & tMaterials = aProtoShell->materials() ;

            for ( uint b=0; b<tMaterials.size(); ++b )
            {
                BELFEM_ERROR( mMaterialMap->key_exists( string_to_lower( tMaterials( b ) ) ),
                    "Layer %u of thin shell '%s' uses the material '%s', which is not defined in the materials section of the input file.",
                    ( unsigned int ) b+1,
                    aProtoShell->label().c_str(),
                    tMaterials( b ).c_str() );
            }

            if ( mCreateGhostFacets )
            {
                // this sequence makes sure that we only create duplicates
                // if the materials differ between layers
                uint tNumMat = aProtoShell->materials().size() ; // equal to number of blocks
                for ( uint b=1; b<tNumMat; ++b )
                {
                    Material * tA = mMaterialMap->find( string_to_lower( aProtoShell->materials()( b-1 ) ) )->second ;
                    Material * tB = mMaterialMap->find( string_to_lower( aProtoShell->materials()( b ) ) )->second ;

                    tLayers( tCount )->hasDuplicates = ( tA != tB ) && tA->have( MaterialProperty::rho ) && tB->have( MaterialProperty::rho );
                    tCount += tOrder ;
                }
            }

            // create the nodes on the layers
            mMesh->unflag_all_nodes( 0 );
            mMesh->unflag_all_nodes( 1 );
            mMesh->unflag_all_nodes( 2 );
            this->create_nodes_on_layers( tNodes, tNodeNormals, tDistances, tLayers );

            Cell< Edge * > tEdges;

            // create the temporary edges that we need to construct the layer edges
            this->create_temporary_edges( tNodes, tFacets, tEdges );
            this->update_node_edge_tables( tNodes, tEdges );
            this->update_node_facet_tables( tNodes, tFacets );

            Cell< Block * > tBlocks ;
            switch ( tElementType )
            {
                case ElementType::LINE2 :
                {
                    this->create_elements_on_blocks_line2( mMaxGroupID, mMaxElementID, tFacets, tLayers, tBlocks );
                    this->create_edges_on_layers( tEdges, mMaxElementID, tLayers );
                    this->link_elements_with_edges( tFacets, tLayers, tBlocks );
                    break ;
                }
                case ElementType::LINE3 :
                {
                    this->create_elements_on_blocks_line3( mMaxGroupID, mMaxElementID, tFacets, tLayers, tBlocks );
                    this->create_edges_on_layers(  tEdges, mMaxElementID, tLayers );
                    this->link_elements_with_edges( tFacets, tLayers, tBlocks );
                    break ;
                }
                case ElementType::TRI3 :
                {
                    this->create_elements_on_blocks_tri3( mMaxGroupID, mMaxElementID, tFacets, tLayers, tBlocks );
                    this->create_edges_on_layers( tEdges, mMaxElementID, tLayers );
                    this->link_elements_with_edges( tFacets, tLayers, tBlocks );
                    break ;
                }
                case ElementType::TRI6 :
                {
                    this->create_elements_on_blocks_tri6( mMaxGroupID, mMaxElementID, tFacets, tLayers, tBlocks );
                    this->create_edges_on_layers( tEdges, mMaxElementID, tLayers );
                    this->create_faces_on_layers( tFacets, mMaxElementID, tLayers );
                    this->link_elements_with_edges( tFacets, tLayers, tBlocks );
                    this->link_elements_with_faces( tLayers, tBlocks );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Element type not supported" );
                }
            }



            for ( Facet * tFacet : tFacets )
            {
                // set facet id to thin shell id
                tFacet->set_sideset_id( tSideSet->id() );
            }


            // this->write_normals_to_mesh( tNodes, tNodeNormals );

            tCount = 0 ;

            Cell< id_t > tBlockIDs( tBlocks.size() );
            for ( Block * tBlock : tBlocks )
            {
                // set block type
                tBlock->set_domain_type( DomainType::ThinShell );

                // add block elements to mesh
                append( mMesh->elements(), tBlock->elements() );

                // set the thickness of the block
                tBlock->set_thickness( aProtoShell->thicknesses()( tCount++ ) );

                tBlockIDs.push( tBlock->id() );
                mMesh->add_block( tBlock );
            }

            this->create_buffers( aProtoShell->materials(), tLayers, tBlocks );

            index_t g = 0 ;
            if ( mCreateGhostFacets )
            {
                this->create_ghost_facets( tOrder, tFacets, tLayers, tBlocks );


                for ( Layer * tLayer : tLayers )
                {
                    g += tLayer->GhostFacets.size() ;
                }
            }

            // ghost
            SideSet * tGhostSideset = nullptr ;
            if ( g > 0 )
            {
                // make a new sideset
                tGhostSideset = new SideSet(  ++mMaxGroupID, g );
                tGhostSideset->set_domain_type( DomainType::Ghost );

                // populate the sideset
                for ( Layer * tLayer : tLayers )
                {
                    for ( Facet * tFacet : tLayer->GhostFacets )
                    {
                        tFacet->set_sideset_id( mMaxGroupID );
                        tGhostSideset->insert_facet( tFacet );
                    }
                }
            }

            if ( mMesh->has_periodicity() )
            {
                this->flag_layer_nodes( tNodes, tLayers );

                SideSet * tPeriodicMaster = this->create_periodic_sideset( tBlocks, true );
                SideSet * tPeriodicSlave = this->create_periodic_sideset( tBlocks, false );

                if ( tPeriodicMaster != nullptr )
                {
                    mMesh->sidesets().push( tPeriodicMaster );
                }

                if ( tPeriodicSlave  != nullptr )
                {
                    mMesh->sidesets().push( tPeriodicSlave );
                }
            }

            // create the thin shell object
            ThinShell * aThinShell = new ThinShell( tSideSet, tGhostSideset );

            // move blocks into thin shell object
            aThinShell->blocks().vector_data() = std::move( tBlocks.vector_data() );

            // set the thicknesses
            aThinShell->set_thicknesses( aProtoShell->thicknesses() );

            // add thin shell to mesh
            mMesh->thin_shells().push( aThinShell );

            // also add sideset with facets to mesh container
            mMesh->sidesets().push( tSideSet );

            if ( tGhostSideset != nullptr )
            {
                mMesh->sidesets().push( tGhostSideset );
            }


            // #BEGIN SIDE CONNECTOR MOD
            if ( mMesh->number_of_dimensions() == 3 )
            {
                // per-tape opt-in from the input file ( edge coating : on );
                // width 0 = derive from the outer stabilizer layer
                mCreateSideConnectors = aProtoShell->edge_coating() ;
                mConnectorWidth       = aProtoShell->edge_coating_width() ;

                if ( mCreateSideConnectors )
                {
                    this->create_side_connectors(
                        tOrder,
                        aProtoShell,
                        tBlockIDs,
                        tLayers,
                        tNodeNormals,
                        tEdges,
                        tFacets,
                        aThinShell );
                }
                else if ( mFuseEdges )
                {
                    index_t tNumNodes = tLayers.first()->Nodes.size() ;

                    Map< key_t, index_t > tEdgeMap ;
                    this->create_edge_map(
                        tNumNodes,
                        tEdges, tEdgeMap );


                    if ( mFuseEdges )
                    {
                        for ( Curve * tCurve : aProtoShell->side_curves() )
                        {
                            // curve pairs of higher order elements are not in the edge map
                            BELFEM_ERROR( tOrder == 1,
                                "side edge fusing is only implemented for first order thin shells" );

                            Cell< index_t > tIndices ;
                            this->compute_side_edge_indices(
                                tCurve, tEdgeMap, tNumNodes, tOrder, tIndices );

                            // single authority per rim station: all interior levels
                            // and the node fusing tie to the same cut-composed
                            // branch the outer-interface anchors resolve to
                            Cell< Node * > tEdgeSources ;
                            Cell< Node * > tNodeSources ;
                            this->compute_side_authority(
                                tCurve, tFacets, tEdges, tIndices,
                                tEdgeSources, tNodeSources );

                            for ( uint l=1; l<tNumLayers-1; l++ )
                            {
                                this->connect_side_edges( tLayers( l ), tEdges, tIndices, tEdgeSources );
                                this->connect_side_nodes( tCurve->nodes(), tNodeSources, tLayers( l )->Nodes );
                            }
                        }
                    }
                }
            }
            
            // #END SIDE CONNECTOR MOD


            for ( Layer * tLayer : tLayers )
            {
                append_move( mMesh->nodes(), tLayer->Nodes );
                if ( tLayer->Edges.size() > 0 )
                {
                    append_move( mMesh->edges(), tLayer->Edges );
                }
                if ( tLayer->EdgeDuplicates.size() > 0 )
                {
                    append_move( mMesh->edges(), tLayer->EdgeDuplicates );
                }
                if ( tLayer->Faces.size() > 0 )
                {
                    append_move( mMesh->faces(), tLayer->Faces );
                }
                if ( tLayer->FaceDuplicates.size() > 0 )
                {
                    append_move( mMesh->faces(), tLayer->FaceDuplicates );
                }
            }

            // delete the temporary layer containers
            for ( Layer * tLayer : tLayers )
            {
                delete tLayer;
            }

            // remove the temporary edges from the facets before deleting
            // them. This must happen after the side connector branch, which
            // reads the facet edge containers to build the edge-to-facet map
            for ( Facet * tFacet : tFacets )
            {
                tFacet->element()->reset_edge_container() ;
            }

            for ( Edge * tEdge : tEdges )
            {
                delete tEdge;
            }

            for ( Node * tNode : tNodes )
            {
                tNode->original()->reset_edge_container();
            }

            // note : we don't delete facets because they belong to the sideset on the mesh
            aThinShell->move_node_indices( mNodeIndices );
            aThinShell->set_materials( aProtoShell->materials() );
            return aThinShell;
        }

        ThinShellFactory::~ThinShellFactory()
        {
        }

        void
        ThinShellFactory::create_side_connectors(
            const uint aOrder,
            Protoshell * aProtoShell,
            const Cell< id_t > & aBlockIDs,
            Cell< Layer * > & aLayers,
            Matrix< real > & aNodeNormals,
            Cell< Edge* > & aEdges,
            Cell< Facet * > & aFacets,
            ThinShell * aThinShell )
        {
            // selective layers are not wired into the new algorithm; fail
            // loudly instead of silently building connectors on every layer
            BELFEM_ERROR( mConnectorsForAllLayers,
                "selective side connector layers are not supported by the new algorithm" );

            BELFEM_ERROR( aOrder == 1,
                "side connectors are only implemented for first order thin shells" );

            // edge-coating policy guards: the side coating deposits in the
            // same electroplating step as the outer stabilizer layers, so top
            // and bottom must agree in material and thickness, and the wall
            // inherits both. Loud errors, never silent skips: with the
            // edge coating active, a missing wall would falsify the physics
            const Cell< string > & tMaterials   = aProtoShell->materials() ;
            const Vector< real > & tThicknesses = aProtoShell->thicknesses() ;

            BELFEM_ERROR( tMaterials.size() >= 3,
                "edge coating requires at least three tape layers, have %u",
                ( unsigned int ) tMaterials.size() );

            BELFEM_ERROR( tThicknesses.length() == tMaterials.size(),
                "edge coating: thickness count %u does not match material count %u",
                ( unsigned int ) tThicknesses.length(),
                ( unsigned int ) tMaterials.size() );

            BELFEM_ERROR( string_to_lower( tMaterials.first() )
                       == string_to_lower( tMaterials.last() ),
                "edge coating requires identical top and bottom layer materials, have %s and %s",
                tMaterials.first().c_str(),
                tMaterials.last().c_str() );

            const real tBottom = tThicknesses( 0 );
            const real tTop    = tThicknesses( tThicknesses.length() - 1 );

            // relative tolerance: the values come from the input file parser
            BELFEM_ERROR( std::abs( tTop - tBottom )
                    <= 1e-6 * std::max( tTop, tBottom ),
                "edge coating requires identical top and bottom layer thicknesses, have %g and %g m",
                tBottom, tTop );

            // zero thickness passes the relative test above ( 0 <= 0 ) and
            // would build a zero-width wall
            BELFEM_ERROR( tBottom > 0.0,
                "edge coating: outer layer thickness must be positive, have %g m",
                tBottom );

            // wall width: derived from the outer layer unless explicitly set
            const real tWidth = mConnectorWidth > 0.0 ? mConnectorWidth : tBottom ;

            int_t tNumNodes = aNodeNormals.n_cols() ;

            Map< key_t, index_t > tEdgeMap ;
            this->create_edge_map(
                tNumNodes,
                aEdges, tEdgeMap );

            Matrix< real > tBinomials ;

            ElementFactory tFactory ;

            uint tNumLayers = aLayers.size();

            Cell< Element * > tLeftCoatings ;
            Cell< Element * > tRightCoatings ;
            id_t tL = ++mMaxGroupID ;
            id_t tR = ++mMaxGroupID ;

            for ( Curve * tCurve : aProtoShell->side_curves() )
            {
                real tSign = this->compute_binomial_vectors( aProtoShell,
                    tCurve,
                    aNodeNormals,
                    tBinomials );

                Cell< index_t > tIndices ;

                Cell< SideLayer * > tSideLayers( aLayers.size() );

                this->compute_side_edge_indices(
                    tCurve, tEdgeMap, tNumNodes, aOrder, tIndices );

                // facet lookup for this curve's side edges only
                Map< index_t, std::pair< index_t, int8_t > > tEdgeToFacetMap ;
                this->create_edge_to_face_map(
                    aEdges, aFacets, tIndices, tEdgeToFacetMap );

                for ( uint l=0; l<tNumLayers; ++l )
                {
                    tSideLayers.push( new SideLayer(
                        tSign,
                        tCurve,
                        aLayers( l ),
                        aEdges,
                        tIndices,
                        tBinomials,
                        tWidth,
                        mMaxNodeID,
                        mMaxElementID,
                        mFuseEdgesWhenHavingSideConnectors ) );
                }

                uint n = tIndices.size() ;
                uint m = aLayers.size() - 1 ;

                Cell< Element * > & tElements = tSign == 1.0 ? tRightCoatings : tLeftCoatings ;
                tElements.reserve( tElements.size() + n*m );


                // @AI: this is a bit risky because we are linking the same sideset to multiple
                //      thin shell blocks. We don't store this sideset in the exodus file to avoid issues
                SideSet * tSideSet = new SideSet( ++mMaxGroupID, n * m );
                tSideSet->hide();

                for ( uint j=0; j<m; ++j )
                {
                    Block * tOtherBlock = mMesh->block( aBlockIDs( j ) );

                    // bottom layer
                    SideLayer * tBm = tSideLayers( j );

                    // top layer
                    SideLayer * tTp = tSideLayers( j + 1 );

                    Cell< Edge * > & tInnerEdgesBottom = tBm->InnerEdgeDuplicates.size() == 0 ? tBm->InnerEdges : tBm->InnerEdgeDuplicates ;
                    Cell< Edge * > & tOuterEdgesBottom = tBm->OuterEdgeDuplicates.size() == 0 ? tBm->OuterEdges : tBm->OuterEdgeDuplicates ;

                    Cell< Edge * > & tA = tSign == 1.0 ? tOuterEdgesBottom : tInnerEdgesBottom ;

                    Cell< Edge * > & tB = tSign == 1.0 ? tInnerEdgesBottom  : tOuterEdgesBottom ;

                    Cell< Edge * > & tC = tSign == 1.0 ? tTp->OuterEdges : tTp->InnerEdges ;
                    Cell< Edge * > & tD = tSign == 1.0 ? tTp->InnerEdges : tTp->OuterEdges ;



                    for ( uint i=0; i<n; ++i )
                    {
                        Element * tElement = tFactory.create_element( ElementType::HEX8TB, ++mMaxElementID );

                        Element * tFace    = tFactory.create_element( ElementType::QUAD4, ++mMaxElementID );

                        // the physical tag stays free for the material machinery
                        // ( set_physical_tags_for_elements ); the layer block is
                        // recovered through the facet's master instead

                        Edge * a = tA( i );
                        Edge * b = tB( i );
                        Edge * c = tC( i );
                        Edge * d = tD( i );
                        Edge * e = aEdges( tIndices( i ) );

                        Edge * f =  tInnerEdgesBottom( i );
                        Edge * g =  tTp->InnerEdges( i );

                        bool tIsReversed = e->node( 1 )->original() == tCurve->nodes()( i )->original();

                        BELFEM_ASSERT( tIsReversed || e->node( 0 )->original() == tCurve->nodes()( i )->original(),
                            "Can't determine edge orientation along curve %lu",
                            ( long unsigned int ) tCurve->id() );

                        if ( tIsReversed )
                        {
                            tElement->insert_node( a->node( 1 ), 0 );
                            tElement->insert_node( a->node( 0 ), 1 );
                            tElement->insert_node( b->node( 0 ), 2 );
                            tElement->insert_node( b->node( 1 ), 3 );
                            tElement->insert_node( c->node( 1 ), 4 );
                            tElement->insert_node( c->node( 0 ), 5 );
                            tElement->insert_node( d->node( 0 ), 6 );
                            tElement->insert_node( d->node( 1 ), 7 );
                        }
                        else
                        {
                            tElement->insert_node( a->node( 0 ), 0 );
                            tElement->insert_node( a->node( 1 ), 1 );
                            tElement->insert_node( b->node( 1 ), 2 );
                            tElement->insert_node( b->node( 0 ), 3 );
                            tElement->insert_node( c->node( 0 ), 4 );
                            tElement->insert_node( c->node( 1 ), 5 );
                            tElement->insert_node( d->node( 1 ), 6 );
                            tElement->insert_node( d->node( 0 ), 7 );
                        }

                        tElement->allocate_edge_container();

                        tElement->insert_edge( a, 0 );
                        tElement->insert_edge( b, 1 );
                        tElement->insert_edge( c, 2 );
                        tElement->insert_edge( d, 3 );

                        tElements.push( tElement );

                        // Get the master element. The key must be the temp edge:
                        // layer edge copies carry no index of their own
                        auto tPair = tEdgeToFacetMap( e->index() );
                        Element * tMaster = tOtherBlock->elements()( tPair.first );


                        Facet * tFacet = new Facet( tFace );

                        // linking writes the master face nodes onto the QUAD4 in
                        // canonical master order, so the facet normal points out of
                        // the shell element toward the wall for both connector signs
                        tFacet->set_master( tMaster, tPair.second, true );

                        tFacet->set_slave( tElement, tSign > 0. ? 2 : 0 );

                        // the offset between the master and slave face cycles depends
                        // on which lateral face of the shell element the station hits
                        // ( PENTA6TS faces 0/1 wind bottom-edge first, face 2 winds
                        // up-leg first ), so it is not a constant per sign; the
                        // node-identity search resolves it exactly, once, here
                        tFacet->compute_orientation();

                        // the recovery facet carries only the two longitudinal dof
                        // edges; the container lives on the wrapped element because
                        // Facet::edge() delegates there. Slots 2 and 3 stay empty
                        // on purpose
                        tFace->allocate_edge_container();
                        tFace->insert_edge( f, 0 );
                        tFace->insert_edge( g, 1 );

                        tSideSet->insert_facet( tFacet );
                    }
                }



                // register once per curve: the sideset spans all layer gaps,
                // and the mesh would delete a multiply-pushed pointer m times
                mMesh->add_sideset( tSideSet );

                for ( SideLayer * tSideLayer : tSideLayers )
                {
                    append_move( mMesh->nodes(), tSideLayer->OuterNodes );
                    append_move( mMesh->edges(), tSideLayer->OuterEdges );
                    append_move( mMesh->nodes(), tSideLayer->InnerNodes );
                    append_move( mMesh->edges(), tSideLayer->InnerEdges );
                    if ( tSideLayer->OuterEdgeDuplicates.size() > 0 )
                    {
                        append_move( mMesh->edges(), tSideLayer->OuterEdgeDuplicates );
                    }
                    if ( tSideLayer->InnerEdgeDuplicates.size() > 0 )
                    {
                        append_move( mMesh->edges(), tSideLayer->InnerEdgeDuplicates );
                    }
                    delete tSideLayer;
                }

                // one block and one sideset per curve, pushed together so that
                // both containers stay index aligned ( the bfm file pairs them
                // by position, see BfmFile::save_thinshell_data )
                aThinShell->side_connector_sidesets().push( tSideSet );
            }

            if ( ! tLeftCoatings.empty() )
            {
                auto * tBlock = new Block( tL, 0 );
                tBlock->set_thickness( tWidth );
                tBlock->set_domain_type( DomainType::LeftCoating );
                tBlock->set_element_counter( tLeftCoatings.size() );
                tBlock->elements() = std::move( tLeftCoatings );
                tBlock->set_thickness( tWidth );
                aThinShell->side_connector_blocks().push( tBlock );
                mMesh->add_block( tBlock );
            }

            if ( ! tRightCoatings.empty() )
            {
                auto * tBlock = new Block( tR, 0 );
                tBlock->set_thickness( tWidth );
                tBlock->set_domain_type( DomainType::RightCoating );
                tBlock->set_element_counter( tRightCoatings.size() );
                tBlock->elements() = std::move( tRightCoatings );
                tBlock->set_thickness( tWidth );
                aThinShell->side_connector_blocks().push( tBlock );
                mMesh->add_block( tBlock );
            }
        }
        void
        ThinShellFactory::reset_node_indices()
        {
            index_t tCount = 0 ;
            for ( Node * tNode : mMasterNodes )
            {
                tNode->set_index( tCount++ );
                BELFEM_ASSERT( ! tNode->is_duplicate(), "Expect only non duplicates in master nodes" );
            }
            tCount = 0 ;
            for ( Node * tNode : mSlaveNodes )
            {
                tNode->set_index( tCount++ );
            }

        }

        void
        ThinShellFactory::write_normals_to_mesh( Cell< Node * > & aNodes, const Matrix< real > & aNodeNormals )
        {
            mMesh->update_node_indices() ;

            Vector< real > & tNx = mMesh->field_exists( "Nx" ) ?
                    mMesh->field_data( "Nx" ) : mMesh->create_field( "Nx" );
            Vector< real > & tNy =  mMesh->field_exists( "Ny" ) ?
                    mMesh->field_data( "Ny" ) : mMesh->create_field( "Ny" );

            Vector< real > & tNz =  mMesh->field_exists( "Nz" ) ?
                    mMesh->field_data( "Nz" ) : mMesh->create_field( "Nz" );

            index_t i = 0 ;
            for ( Node * tNode : aNodes )
            {
                index_t k = tNode->index();

                tNx( k ) = aNodeNormals( 0, i );
                tNy( k ) = aNodeNormals( 1, i );
                tNz( k ) = aNodeNormals( 2, i );
                ++i ;
            }
        }

        void
        ThinShellFactory::collect_sidesets( const Vector< id_t > & aSideSetIDs, Cell< SideSet * > & aSideSets )
        {
            index_t tCount = 0 ;
            for ( id_t tID : aSideSetIDs )
            {
                if ( mMesh->sideset( tID )->number_of_facets() > 0 )
                {
                    ++tCount;
                }
            }
            aSideSets.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( id_t tID : aSideSetIDs )
            {
                if ( mMesh->sideset( tID )->number_of_facets() > 0 )
                {
                    aSideSets( tCount++ ) = mMesh->sideset( tID );
                }
            }
        }


        ElementType
        ThinShellFactory::collect_facets( Cell< SideSet * > & aSideSets, Cell< Facet * > & aFacets )
        {
            index_t tCount = 0 ;

            BELFEM_ASSERT( aSideSets.size() > 0, "No sidesets defined" );
            ElementType aType = ElementType::UNDEFINED ;

            // note: every facet of the selected sidesets is taken; the sidesets' own
            //       containers are emptied below ( reset_facet_container )

            for( SideSet * tSideSet : aSideSets )
            {
                Cell< Facet * > & tFacets = tSideSet->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    tFacet->set_index( tCount++ );
                }

                if( tSideSet->number_of_facets() > 0 )
                {
                    ElementType tType = tSideSet->facets()(0)->element()->type();
                    if( aType == ElementType::UNDEFINED )
                    {
                        aType = tType ;
                    }
                    else if( aType != tType )
                    {
                        BELFEM_ERROR( false, "Sidesets must have the same element type" );
                    }
                }
            }

            aFacets.set_size( tCount, nullptr );
            tCount = 0 ;
            for( SideSet * tSideSet : aSideSets )
            {
                Cell< Facet * > & tFacets = tSideSet->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    aFacets( tCount++ ) = tFacet ;
                }
                tSideSet->reset_facet_container();
            }

            // finally, we make sure that all facets are flagged,
            // which is needed later for the normal computation
            mMesh->unflag_all_facets() ;
            for ( Facet * tFacet : aFacets )
            {
                tFacet->flag();
            }
            return aType ;
        }

        void
        ThinShellFactory::collect_nodes(
            Cell< Facet * > & aFacets, Cell< Node * > & aNodes )
        {
            // select nodes
            DynamicBitset tBitset( mMasterNodes.size() );

            // flag nodes and count facets
            for( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    Node * tNode = tFacet->node( k );

                    // a facet node whose original is outside the master list means
                    // the tape facets were rewired after the cut factory restored
                    // them, e.g. by a cut duplicate or abstract node copied in from
                    // a jump-side master element
                    BELFEM_ERROR( tNode->original()->index() < mMasterNodes.size(),
                                  "Facet %lu carries node %lu whose original has index %lu (expect < %lu master nodes)",
                                  ( long unsigned int ) tFacet->id(),
                                  ( long unsigned int ) tNode->id(),
                                  ( long unsigned int ) tNode->original()->index(),
                                  ( long unsigned int ) mMasterNodes.size() );

                    tBitset.set( tNode->original()->index() );
                }
            }

            tBitset.where( mNodeIndices );

            // populate node container
            aNodes.set_size( mNodeIndices.size(), nullptr );
            index_t tCount = 0 ;
            for ( Node * tNode : mMasterNodes )
            {
                tNode->set_index( gNoIndex );
            }

            for ( index_t k : mNodeIndices )
            {
                aNodes( tCount++ ) = mMasterNodes( k ) ;
            }

            // ordering the nodes in this container should make sure
            // that the edge orientations are preserved
            // during thin shell creation
            sort( aNodes, opVertexID );

            tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tCount++ );
            }
        }

        void
        ThinShellFactory::process_nodes_line2(
            const Cell< Facet * > & aFacets,
            const Cell< Node * >  & aNodes,
                Matrix< real >    & aNodeNormals )
        {
            // the node coordinates
            Vector< real > tA( 2 );
            Vector< real > tB( 2 );

            // the normal
            Vector< real > tN( 3, 0.0 );

            // the center
            Vector< real > tM( 2 );

            index_t tNumFacets = aFacets.size();

            Matrix< real > tFacetNormals;
            Matrix< real > tFacetCenters;

            tFacetNormals.set_size( 3 , tNumFacets );
            tFacetCenters.set_size( 2 , tNumFacets );

            index_t tCount = 0 ;

            for ( Facet * tFacet : aFacets )
            {
                tA( 0 ) = tFacet->node( 0 )->x();
                tA( 1 ) = tFacet->node( 0 )->y();

                tB( 0 ) = tFacet->node( 1 )->x();
                tB( 1 ) = tFacet->node( 1 )->y();

                // center of line
                tM = tA + tB ;
                tM *= 0.5 ;

                tN( 0 ) = tB( 1 ) - tA( 1 );
                tN( 1 ) = tA( 0 ) - tB( 0 );

                tN /= norm( tN );

                tFacetCenters.set_col( tCount, tM );
                tFacetNormals.set_col( tCount, tN );

                ++tCount;
            }

            aNodeNormals.set_size( 3, aNodes.size() );

            Vector< real > tX( 2 );

            for ( Node * tNode : aNodes )
            {
                if ( tNode->number_of_facets() == 1 && tNode->number_of_duplicates() == 0 )
                {
                    aNodeNormals.set_col( tNode->index(), tFacetNormals.col( tNode->facet( 0 )->index() ) );
                    continue;
                }

                tX( 0 ) = tNode->x();
                tX( 1 ) = tNode->y();

                tN.fill( 0.0 );

                for ( uint f=0; f<tNode->number_of_facets(); ++f )
                {
                    if ( tNode->facet( f )->is_flagged() )
                    {
                        index_t k = tNode->facet( f )->index();
                        real tR = norm( tX - tFacetCenters.col( k ) );
                        tN += tFacetNormals.col( k ) / ( tR * tR );
                    }
                }
                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    Node * tDup = tNode->duplicate( d );
                    for ( uint f=0; f<tDup->number_of_facets(); ++f )
                    {
                        if ( tDup->facet( f )->is_flagged() )
                        {
                            index_t k = tDup->facet( f )->index();
                            real tR = norm( tX - tFacetCenters.col( k ) );
                            tN += tFacetNormals.col( k ) / ( tR * tR );
                        }
                    }
                }
                tN /= norm( tN );
                aNodeNormals.set_col( tNode->index(), tN );
            }
        }

        void
        ThinShellFactory::process_nodes_line3(
            const Cell< Facet * > & aFacets,
            const Cell< Node * >  & aNodes,
            Matrix< real >        & aNodeNormals )
        {
            // node 0
            Vector< real > tA( 2 );

            // node 1
            Vector< real > tB( 2 );

            // node 2
            Vector< real > tC( 2 );

            // normal
            Vector< real > tN( 2 );

            aNodeNormals.set_size( 3, aNodes.size(), 0.0 );

            real tR ;
            index_t i ;

            // process the facets
            for ( Facet * tFacet : aFacets )
            {
                // grab the coordinates
                tA( 0 ) = tFacet->node( 0 )->x();
                tA( 1 ) = tFacet->node( 0 )->y();

                tB( 0 ) = tFacet->node( 1 )->x();
                tB( 1 ) = tFacet->node( 1 )->y();

                tC( 0 ) = tFacet->node( 2 )->x();
                tC( 1 ) = tFacet->node( 2 )->y();

                // node 0
                tR = norm( tC-tA );
                tN( 0 ) = 4.*tC(1) - tB(1) - 3.*tA(1) ; //  dy/dxi
                tN( 1 ) =  tB(0) + 3.*tA(0) - 4.*tC(0) ; // -dx/dxi
                tN *= ( tR * tR ) / norm( tN );
                i = tFacet->node( 0 )->original()->index();
                aNodeNormals( i, 0 ) += tN( 0 );
                aNodeNormals( i, 1 ) += tN( 1 );

                // node 1 :
                tR = norm( tC-tB );
                tN( 0 ) =  tA(1) + 3.*tB(1) - 4.*tC(1) ;
                tN( 1 ) =  4.*tC(0)-tA(0) - 3.*tB(0) ;
                tN *= ( tR * tR ) / norm( tN );
                i = tFacet->node( 1 )->original()->index() ;
                aNodeNormals( i, 0 ) += tN( 0 );
                aNodeNormals( i, 1 ) += tN( 1 );

                // node 2
                tN( 0 ) =  tB(1) - tA(1) ;
                tN( 1 ) =  tA(0) - tB(0) ;
                tN /= norm( tN );
                i = tFacet->node( 2 )->original()->index() ;
                aNodeNormals( i, 0 ) += tN( 0 );
                aNodeNormals( i, 1 ) += tN( 1 );
            }

            // normalize the normals
            for ( Node * tNode : aNodes )
            {
                tN = aNodeNormals.col( tNode->index() );
                tN /= norm( tN );
                aNodeNormals.set_col( tNode->index(), tN );
            }
        }

        void
        ThinShellFactory::process_nodes_tri3(
            const Cell< Facet * > & aFacets,
            const Cell< Node * >  & aNodes,
                Matrix< real > & aNodeNormals )
        {
            // the node coordinates
            Vector< real > tA( 3 );
            Vector< real > tB( 3 );
            Vector< real > tC( 3 );

            // the first and the second edge
            Vector< real > tP( 3 );
            Vector< real > tQ( 3 );

            // the normal
            Vector< real > tN( 3 );

            // the center
            Vector< real > tM( 3 );

            index_t tNumFacets = aFacets.size();

            Matrix< real > tFacetNormals;
            Matrix< real > tFacetCenters;

            tFacetNormals.set_size( 3 , tNumFacets );
            tFacetCenters.set_size( 3 , tNumFacets );

            index_t tCount = 0 ;

            for ( Facet * tFacet : aFacets )
            {
                tA( 0 ) = tFacet->node( 0 )->x();
                tA( 1 ) = tFacet->node( 0 )->y();
                tA( 2 ) = tFacet->node( 0 )->z();

                tB( 0 ) = tFacet->node( 1 )->x();
                tB( 1 ) = tFacet->node( 1 )->y();
                tB( 2 ) = tFacet->node( 1 )->z();

                tC( 0 ) = tFacet->node( 2 )->x();
                tC( 1 ) = tFacet->node( 2 )->y();
                tC( 2 ) = tFacet->node( 2 )->z();

                // center of triangle
                tM = tA + tB + tC ;
                tM /= 3.0;

                tP = tB - tA;
                tQ = tC - tA;
                tN = cross( tP, tQ );

                tN /= norm( tN );

                tFacetCenters.set_col( tCount, tM );
                tFacetNormals.set_col( tCount, tN );

                ++tCount;
            }

            aNodeNormals.set_size( 3, aNodes.size() );

            Vector< real > tX( 3 );

            for ( Node * tNode : aNodes )
            {
                if ( tNode->number_of_facets() == 1 && tNode->number_of_duplicates() == 0 )
                {
                    aNodeNormals.set_col( tNode->index(), tFacetNormals.col( tNode->facet( 0 )->index() ) );
                    continue;
                }

                tX( 0 ) = tNode->x();
                tX( 1 ) = tNode->y();
                tX( 2 ) = tNode->z();

                tN.fill( 0.0 );

                for ( uint f=0; f<tNode->number_of_facets(); ++f )
                {
                    if ( tNode->facet( f )->is_flagged() )
                    {
                        index_t k = tNode->facet( f )->index();
                        real tR = norm( tX - tFacetCenters.col( k ) );
                        tN += tFacetNormals.col( k ) / ( tR * tR );
                    }
                }

                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    Node * tDup = tNode->duplicate( d );
                    for ( uint f=0; f<tDup->number_of_facets(); ++f )
                    {
                        if ( tDup->facet( f )->is_flagged() )
                        {
                            index_t k = tDup->facet( f )->index();
                            real tR = norm( tX - tFacetCenters.col( k ) );
                            tN += tFacetNormals.col( k ) / ( tR * tR );
                        }
                    }
                }
                tN /= norm( tN );
                aNodeNormals.set_col( tNode->index(), tN );
            }
        }

        void
        ThinShellFactory::process_nodes_tri6(
            const Cell< Facet * > & aFacets,
            const Cell< Node * > & aNodes,
                  Matrix< real > & aNodeNormals )
        {
            Cell< Matrix< real > > tFacetNormals( aFacets.size(), {{}});
            Matrix< real > tFacetCenters;
            tFacetCenters.set_size( 3, aFacets.size() );

            // node coordinates
            Vector< real > tX( 6 );
            Vector< real > tY( 6 );
            Vector< real > tZ( 6 );

            // normals
            Vector< real  > tN(3);

            // centers
            Vector< real  > tM(3);

            // coordinaters of node
            Vector< real > tP( 3 );

            // allocate memory
            aNodeNormals.set_size( 3, aNodes.size() );

            // process the facets
            for ( Facet * tFacet : aFacets )
            {
                // grab the coordinates
                for ( uint k=0; k<6; ++k )
                {
                    tX( k ) = tFacet->node( k )->x();
                    tY( k ) = tFacet->node( k )->y();
                    tZ( k ) = tFacet->node( k )->z();
                }

                // get the normal vector
                Matrix< real > & tMat = tFacetNormals( tFacet->index() );

                tMat.set_size( 3, 6 );

                tN(0)=(3.*tZ(0)+tZ(2)-4.*tZ(5))*(tY(1)-tY(2)+4.*(tY(5)-tY(3)))-(3.*tY(0)+tY(2)-4.*tY(5))*(tZ(1)-tZ(2)+4.*(tZ(5)-tZ(3)));
                tN(1)=(3.*tX(0)+tX(2)-4.*tX(5))*(tZ(1)-tZ(2)+4.*(tZ(5)-tZ(3)))-(3.*tZ(0)+tZ(2)-4.*tZ(5))*(tX(1)-tX(2)+4.*(tX(5)-tX(3)));
                tN(2)=(3.*tY(0)+tY(2)-4.*tY(5))*(tX(1)-tX(2)+4.*(tX(5)-tX(3)))-(3.*tX(0)+tX(2)-4.*tX(5))*(tY(1)-tY(2)+4.*(tY(5)-tY(3)));
                tN/=norm(tN);
                tMat.set_col( 0, tN );

                tN(0)=(3.*tY(1)+tY(2)-4.*tY(4))*(tZ(0)-tZ(2)+4.*(tZ(4)-tZ(3)))-(3.*tZ(1)+tZ(2)-4.*tZ(4))*(tY(0)-tY(2)+4.*(tY(4)-tY(3)));
                tN(1)=(3.*tZ(1)+tZ(2)-4.*tZ(4))*(tX(0)-tX(2)+4.*(tX(4)-tX(3)))-(3.*tX(1)+tX(2)-4.*tX(4))*(tZ(0)-tZ(2)+4.*(tZ(4)-tZ(3)));
                tN(2)=(3.*tX(1)+tX(2)-4.*tX(4))*(tY(0)-tY(2)+4.*(tY(4)-tY(3)))-(3.*tY(1)+tY(2)-4.*tY(4))*(tX(0)-tX(2)+4.*(tX(4)-tX(3)));
                tN/=norm(tN);
                tMat.set_col( 1, tN );

                tN(0)=(tY(0)+3.*tY(2)-4.*tY(5))*(tZ(1)+3.*tZ(2)-4.*tZ(4))-(tY(1)+3.*tY(2)-4.*tY(4))*(tZ(0)+3.*tZ(2)-4.*tZ(5));
                tN(1)=(tX(1)+3.*tX(2)-4.*tX(4))*(tZ(0)+3.*tZ(2)-4.*tZ(5))-(tX(0)+3.*tX(2)-4.*tX(5))*(tZ(1)+3.*tZ(2)-4.*tZ(4));
                tN(2)=(tX(0)+3.*tX(2)-4.*tX(5))*(tY(1)+3.*tY(2)-4.*tY(4))-(tX(1)+3.*tX(2)-4.*tX(4))*(tY(0)+3.*tY(2)-4.*tY(5));

                tN/=norm(tN);
                tMat.set_col( 2, tN );

                tN(0)=(tY(0)+tY(2)+2.*(tY(3)-tY(4)-tY(5)))*(tZ(1)+tZ(2)+2.*(tZ(3)-tZ(4)-tZ(5)))-(tY(1)+tY(2)+2.*(tY(3)-tY(4)-tY(5)))*(tZ(0)+tZ(2)+2.*(tZ(3)-tZ(4)-tZ(5)));
                tN(1)=(tX(1)+tX(2)+2.*(tX(3)-tX(4)-tX(5)))*(tZ(0)+tZ(2)+2.*(tZ(3)-tZ(4)-tZ(5)))-(tX(0)+tX(2)+2.*(tX(3)-tX(4)-tX(5)))*(tZ(1)+tZ(2)+2.*(tZ(3)-tZ(4)-tZ(5)));
                tN(2)=(tX(0)+tX(2)+2.*(tX(3)-tX(4)-tX(5)))*(tY(1)+tY(2)+2.*(tY(3)-tY(4)-tY(5)))-(tX(1)+tX(2)+2.*(tX(3)-tX(4)-tX(5)))*(tY(0)+tY(2)+2.*(tY(3)-tY(4)-tY(5)));
                tN/=norm(tN);
                tMat.set_col( 3, tN );

                tN(0)=(tY(1)-tY(2))*(tZ(0)+tZ(2)+2.*(tZ(4)-tZ(3)-tZ(5)))-(tZ(1)-tZ(2))*(tY(0)+tY(2)+2.*(tY(4)-tY(3)-tY(5)));
                tN(1)=(tZ(1)-tZ(2))*(tX(0)+tX(2)+2.*(tX(4)-tX(3)-tX(5)))-(tX(1)-tX(2))*(tZ(0)+tZ(2)+2.*(tZ(4)-tZ(3)-tZ(5)));
                tN(2)=(tX(1)-tX(2))*(tY(0)+tY(2)+2.*(tY(4)-tY(3)-tY(5)))-(tY(1)-tY(2))*(tX(0)+tX(2)+2.*(tX(4)-tX(3)-tX(5)));
                tN/=norm(tN);
                tMat.set_col( 4, tN );

                tN(0)=(tZ(0)-tZ(2))*(tY(1)+tY(2)+2.*(tY(5)-tY(3)-tY(4)))-(tY(0)-tY(2))*(tZ(1)+tZ(2)+2.*(tZ(5)-tZ(3)-tZ(4)));
                tN(1)=(tX(0)-tX(2))*(tZ(1)+tZ(2)+2.*(tZ(5)-tZ(3)-tZ(4)))-(tZ(0)-tZ(2))*(tX(1)+tX(2)+2.*(tX(5)-tX(3)-tX(4)));
                tN(2)=(tY(0)-tY(2))*(tX(1)+tX(2)+2.*(tX(5)-tX(3)-tX(4)))-(tX(0)-tX(2))*(tY(1)+tY(2)+2.*(tY(5)-tY(3)-tY(4)));
                tN/=norm(tN);
                tMat.set_col( 5, tN );

                // center
                tM(0) = (4.*(tX(3)+tX(4)+tX(5))-tX(1)-tX(2)-tX(0))/9.;
                tM(0) = (4.*(tY(3)+tY(4)+tY(5))-tY(1)-tY(2)-tY(0))/9.;
                tM(0) = (4.*(tZ(3)+tZ(4)+tZ(5))-tZ(1)-tZ(2)-tZ(0))/9.;

                tFacetCenters.set_col( tFacet->index(), tM );
            }


            // loop over all nodes
            for ( mesh::Node * tNode : aNodes )
            {
                tN.fill( 0.0 );

                tP( 0 ) = tNode->x();
                tP( 1 ) = tNode->y();
                tP( 2 ) = tNode->z();

                // loop over all facets
                for ( uint f=0; f<tNode->number_of_facets(); ++f )
                {
                    Facet * tFacet = tNode->facet( f );
                    if ( tNode->facet( f )->is_flagged() )
                    {
                        // get the matrix for the facet
                        Matrix< real > & tMat = tFacetNormals( tFacet->index() );

                        // find which node it is
                        for ( uint i=0; i<tFacet->number_of_nodes(); ++i )
                        {
                            if ( tFacet->node( i )->original()->id() == tNode->original()->id() )
                            {
                                // compute the distance of this facet to the node
                                real tR = norm( tP - tFacetCenters.col( tFacet->index() ) );

                                // add the normal value scaled with the distance to the center
                                tN += tMat.col( i ) / ( tR * tR );
                                break ;
                            }
                        }
                    }
                }
                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    Node * tDup = tNode->duplicate( d );
                    for ( uint f=0; f<tDup->number_of_facets(); ++f )
                    {
                        Facet * tFacet = tDup->facet( f );
                        if ( tDup->facet( f )->is_flagged() )
                        {
                            // get the matrix for the facet
                            Matrix< real > & tMat = tFacetNormals( tFacet->index() );

                            // find which node it is
                            for ( uint i=0; i<tFacet->number_of_nodes(); ++i )
                            {
                                if ( tFacet->node( i )->original()->id() == tNode->original()->id() )
                                {
                                    // compute the distance of this facet to the node
                                    real tR = norm( tP - tFacetCenters.col( tFacet->index() ) );

                                    // add the normal value scaled with the distance to the center
                                    tN += tMat.col( i ) / ( tR * tR );
                                    break ;
                                }
                            }
                        }
                    }
                }
                // normalize
                tN /= norm( tN );
                tP /= norm( tP );

                aNodeNormals.set_col( tNode->index(), tN );
            }
        }

        void
        ThinShellFactory::compute_distances( const uint aOrder, const Vector< real > & aThicknesses, Vector< real > & aDistances )
        {
            uint tNumLayers = aOrder * aThicknesses.length() + 1 ;

            aDistances.set_size( tNumLayers, 0.0 );

            if ( aOrder == 1 )
            {
                for ( uint i=0; i<aThicknesses.length(); ++i )
                {
                    aDistances( i+1 ) = aDistances( i ) + aThicknesses( i );
                }
            }
            else if ( aOrder == 2 )
            {
                uint tCount = 1 ;
                real tB = 0.0 ;

                for ( real tT : aThicknesses )
                {
                    real tA = tB ; // shift

                    // center
                    aDistances( tCount++ ) = tA + 0.5 * tT ;

                    // next
                    tB = tA + tT ;
                    aDistances( tCount++ ) = tB ;
                }
            }


            aDistances -= 0.5 * aDistances( tNumLayers-1 );
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_nodes_on_layers(
                  Cell< Node * >  & aNodes,
            const Matrix< real >  & aNodeNormals,
            const Vector< real >  & aDistances,
                  Cell< Layer * > & aLayers )
        {
            index_t tNumNodes = aNodes.size();
            index_t tCount = 0 ;

            // count nodes that have duplicates
            Cell< uint > tNumDuplicates( tNumNodes, 0 );

            const real tTolerance
                = 0.5 * std::abs( aDistances( 0 ) - aDistances( aDistances.length()-1 ) )
                  + BELFEM_MESH_EPSILON ;

            for ( Node * tNode : aNodes )
            {
                Node * tOrg = tNode->original();

                if ( tOrg->number_of_duplicates() > 1 )
                {
                    tCount = 0 ;

                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        if ( tOrg->duplicate( d )->is_flagged() && tOrg->id() != tNode->id() )
                        {
                            ++tCount ;
                        }
                    }
                    tNumDuplicates( tNode->index() ) = tCount ;
                }
            }

            tCount = 0 ;
            for ( index_t k=0; k<tNumNodes; ++k )
            {
                if ( tNumDuplicates( k ) > 0 )
                {
                    ++tCount ;
                }
            }

            // extrude nodes
            id_t & tID = mMaxNodeID ;
            tID = this->max_node_id() ;

            Vector< real > tX( 3, 0.0 );

            // loop over all layers
            for ( uint l=0; l<aLayers.size(); ++l )
            {
                Cell< Node * > & tNodes = aLayers(l)->Nodes ;
                tNodes.set_size( tNumNodes, nullptr );

                // get the distance
                real tD = aDistances( l );

                // create the new nodes
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    // get the original node coordinates
                    for ( uint i=0; i<mNumDimensions; ++i )
                    {
                        tX( i ) = aNodes( k )->x( i );
                    }

                    // add the node normal
                    tX += tD * aNodeNormals.col( k );

                    // create the node

                    Node * tNode = new Node( ++tID, tX( 0 ), tX( 1 ), tX( 2 ) );
                    tNode->set_index( k );
                    tNodes( k ) = tNode  ;
                }

                // link duplicates
                if ( tCount > 0 )
                {
                    for ( index_t k=0; k<tNumNodes; ++k )
                    {
                        if ( tNumDuplicates( k ) > 0 )
                        {
                            Node * tRef = aNodes( k )->original();
                            Node * tOrg = tNodes( k );
                            tOrg->allocate_duplicate_container( tNumDuplicates( k ) );

                            for ( uint d=0; d<tRef->number_of_duplicates(); ++d )
                            {
                                if ( tRef->duplicate( d )->is_flagged() && tRef->id() != aNodes( k )->id() )
                                {
                                    Node * tDup = tNodes( tRef->duplicate( d )->index() );
                                    tOrg->add_duplicate( tDup );
                                    tDup->set_original( tOrg );
                                }
                            }
                        }
                    }
                }
            }

            if ( mMesh->has_periodicity() )
            {
                DynamicBitset tBitset( tNumNodes );
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    if ( aNodes( k )->is_periodic() )
                    {
                        tBitset.set( k );
                    }
                }
                Cell< index_t > tIndices ;
                tBitset.where( tIndices );

                // at this points, the new nodes are not flagged
                // so we flag them during processing to avoid redundancy
                Periodicity * tPeriodicity = mMesh->periodicity();

                for ( uint l=0; l<aLayers.size(); ++l )
                {
                    Cell< Node * > & tNodes = aLayers(l)->Nodes ;
                    for ( index_t k : tIndices )
                    {

                        Node * tC  = tNodes( k );
                        if ( tC->is_flagged() ) continue ;

                        Node * tA  = aNodes( k );
                        Node * tB  = tA->periodic();

                        Node * tD  = tNodes( tB->index() );

                        tC->flag();
                        tD->flag();
                        tC->set_periodic( tD );
                        tD->set_periodic( tC );

                        tPeriodicity->add_node_pair_to_backup( tC, tD, tTolerance );

                        // flags for master and slave sets
                        if ( tA->is_flagged( 1 ) ) tC->flag( 1 );
                        if ( tA->is_flagged( 2 ) ) tC->flag( 2 );
                        if ( tB->is_flagged( 1 ) ) tD->flag( 1 );
                        if ( tB->is_flagged( 2 ) ) tD->flag( 2 );
                    }
                }
            }

            // restore node indices
            /*tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tNodeIndices( tCount++ ) );
                tNode->unflag();
            }*/
        }


//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_elements_on_blocks_line2(
            id_t & aBlockID,
            id_t & aElementID,
            Cell< Facet * > & aFacets,
            Cell< Layer * > & aLayers,
            Cell< Block * > & aBlocks  )
        {
            index_t tNumElements = aFacets.size() ;
            index_t tNumBlocks   = aLayers.size() - 1 ;

            aBlocks.set_size( tNumBlocks, nullptr );

            ElementFactory tFactory ;

            // loop over all facets
            for ( index_t b=0; b<tNumBlocks; ++b )
            {
                Layer * tBottom = aLayers( b );
                Layer * tTop = aLayers( b+1 );

                Block * tBlock = new Block( ++aBlockID, tNumElements);
                tBlock->label() = "edge";

                for ( Facet * tFacet : aFacets )
                {
                    Element * tElement = tFactory.create_element( ElementType::QUAD4TS, ++aElementID );

                    tElement->insert_node( tBottom->Nodes( tFacet->node( 0 )->original()->index() ), 0 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 1 )->original()->index() ), 1 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 1 )->original()->index() ), 2 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 0 )->original()->index() ), 3 );

                    tElement->set_block_id( aBlockID );

                    tBlock->insert_element( tElement );
                }
                aBlocks( b ) = tBlock ;
            }
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_elements_on_blocks_line3(
            id_t & aBlockID,
            id_t & aElementID,
            Cell< Facet * > & aFacets,
            Cell< Layer * > & aLayers,
            Cell< Block * > & aBlocks  )
        {
            index_t tNumElements = aFacets.size() ;
            index_t tNumBlocks   = aLayers.size()/2 - 1 ;

            aBlocks.set_size( tNumBlocks, nullptr );

            uint tCount = 0 ;

            ElementFactory tFactory ;

            // loop over all facets
            for ( index_t b=0; b<tNumBlocks; ++b )
            {
                Layer * tBottom = aLayers( tCount );
                Layer * tMid = aLayers( tCount+1 );
                Layer * tTop = aLayers( tCount+2 );
                tCount += 2 ;

                Block * tBlock = new Block( ++aBlockID, tNumElements);

                for ( Facet * tFacet : aFacets )
                {
                    Element * tElement = tFactory.create_element( ElementType::QUAD9TS, ++aElementID );

                    tElement->insert_node( tBottom->Nodes( tFacet->node( 0 )->original()->index() ), 0 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 1 )->original()->index() ), 1 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 1 )->original()->index() ), 2 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 0 )->original()->index() ), 3 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 2 )->original()->index() ), 4 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 1 )->original()->index() ), 5 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 2 )->original()->index() ), 6 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 0 )->original()->index() ), 7 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 2 )->original()->index() ), 8 );

                    tElement->set_block_id( aBlockID );

                    tBlock->insert_element( tElement );
                }
                aBlocks( b ) = tBlock ;
            }
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_elements_on_blocks_tri3(
            id_t & aBlockID,
            id_t & aElementID,
            Cell< Facet * > & aFacets,
            Cell< Layer * > & aLayers,
            Cell< Block * > & aBlocks  )
        {
            index_t tNumElements = aFacets.size() ;
            index_t tNumBlocks   = aLayers.size() - 1 ;

            aBlocks.set_size( tNumBlocks, nullptr );

            ElementFactory tFactory ;

            for ( index_t b=0; b<tNumBlocks; ++b )
            {
                Layer * tBottom = aLayers( b );
                Layer * tTop = aLayers( b+1 );

                Block * tBlock = new Block( ++aBlockID, tNumElements);

                for ( Facet * tFacet : aFacets )
                {
                    Element * tElement = tFactory.create_element( ElementType::PENTA6TS, ++aElementID );

                    tElement->insert_node( tBottom->Nodes( tFacet->node( 0 )->original()->index() ), 0 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 1 )->original()->index() ), 1 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 2 )->original()->index() ), 2 );

                    tElement->insert_node( tTop->Nodes( tFacet->node( 0 )->original()->index() ), 3 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 1 )->original()->index() ), 4 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 2 )->original()->index() ), 5 );

                    tElement->set_block_id( aBlockID );

                    tBlock->insert_element( tElement );
                }
                aBlocks( b ) = tBlock ;
            }
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_elements_on_blocks_tri6(
              id_t & aBlockID,
              id_t & aElementID,
              Cell< Facet * > & aFacets,
              Cell< Layer * > & aLayers,
              Cell< Block * > & aBlocks  )
        {
            index_t tNumElements = aFacets.size() ;
            index_t tNumBlocks   = aLayers.size()/2 - 1 ;

            aBlocks.set_size( tNumBlocks, nullptr );

            uint tCount = 0 ;

            ElementFactory tFactory ;

            // loop over all facets
            for ( index_t b=0; b<tNumBlocks; ++b )
            {
                Layer * tBottom = aLayers( tCount );
                Layer * tMid = aLayers( tCount+1 );
                Layer * tTop = aLayers( tCount+2 );
                tCount += 2 ;
                Block * tBlock = new Block( ++aBlockID, tNumElements);

                for ( Facet * tFacet : aFacets )
                {
                    Element * tElement = tFactory.create_element( ElementType::PENTA18TS, ++aElementID );

                    tElement->insert_node( tBottom->Nodes( tFacet->node( 0 )->original()->index() ), 0 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 1 )->original()->index() ), 1 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 2 )->original()->index() ), 2 );

                    tElement->insert_node( tTop->Nodes( tFacet->node( 0 )->original()->index() ), 3 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 1 )->original()->index() ), 4 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 2 )->original()->index() ), 5 );

                    tElement->insert_node( tBottom->Nodes( tFacet->node( 3 )->original()->index() ), 6 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 4 )->original()->index() ), 7 );
                    tElement->insert_node( tBottom->Nodes( tFacet->node( 5 )->original()->index() ), 8 );

                    tElement->insert_node( tMid->Nodes( tFacet->node( 0 )->original()->index() ), 9 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 1 )->original()->index() ), 10 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 2 )->original()->index() ), 11 );

                    tElement->insert_node( tTop->Nodes( tFacet->node( 3 )->original()->index() ), 12 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 4 )->original()->index() ), 13 );
                    tElement->insert_node( tTop->Nodes( tFacet->node( 5 )->original()->index() ), 14 );

                    tElement->insert_node( tMid->Nodes( tFacet->node( 3 )->original()->index() ), 15 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 4 )->original()->index() ), 16 );
                    tElement->insert_node( tMid->Nodes( tFacet->node( 5 )->original()->index() ), 17 );

                    tElement->set_block_id( aBlockID );

                    tBlock->insert_element( tElement );
                }
                aBlocks( b ) = tBlock ;
            }
        }

//------------------------------------------------------------------------------
        void
        ThinShellFactory::create_temporary_edges(
                  Cell< Node * >   & aNodes,
                  Cell< Facet * >  & aFacets,
                  Cell< Edge* >    & aEdges )
        {
            key_t tNumNodes = aNodes.size() ;

            uint tNumEdgesPerFacet = aFacets( 0 )->element()->number_of_edges();

            Vector< key_t > tEdgeKeys(
            tNumEdgesPerFacet  * aFacets.size(), 0 );

            index_t tCount = 0 ;
            Cell< Node * > tNodes( tNumEdgesPerFacet, nullptr );

            // create the keys
            for ( Facet * tFacet : aFacets )
            {

                for ( uint e=0; e<tNumEdgesPerFacet; ++e )
                {
                    tFacet->element()->get_nodes_of_edge( e, tNodes );

                    key_t tA = tNodes( 0 )->original()->index();
                    key_t tB = tNodes( 1 )->original()->index();

                    key_t tKey = tA  > tB ? tA * tNumNodes + tB : tB * tNumNodes + tA;

                    tEdgeKeys( tCount++ ) = tKey ;
                }
            }

            unique( tEdgeKeys );

            aEdges.set_size( tEdgeKeys.length(), nullptr );

            uint tNumNodesPerEdge = interpolation_order_numeric(  aFacets( 0 )->element()->type() ) + 1 ;

            tCount = 0 ;

            for ( key_t tKey : tEdgeKeys )
            {
                Node * tA = aNodes( tKey % tNumNodes );
                Node * tB = aNodes( tKey / tNumNodes );

                Edge * tEdge = new Edge();
                tEdge->allocate_node_container( tNumNodesPerEdge );
                if ( tA->id() < tB->id() )
                {
                    tEdge->insert_node( tA, 0 );
                    tEdge->insert_node( tB, 1 );
                }
                else
                {
                    tEdge->insert_node( tB, 0 );
                    tEdge->insert_node( tA, 1 );
                }
                tEdge->set_index( tCount );
                aEdges( tCount++ ) = tEdge ;
            }

            // create a map for the keys
            tCount = 0 ;
            mEdgeMap.clear();
            for ( key_t tKey : tEdgeKeys )
            {
                mEdgeMap[ tKey ] = tCount++ ;
            }

            for ( Facet * tFacet : aFacets )
            {
                tFacet->element()->allocate_edge_container();

                for ( uint e=0; e<tNumEdgesPerFacet; ++e )
                {
                    tFacet->element()->get_nodes_of_edge( e, tNodes );

                    key_t tA = tNodes( 0 )->original()->index();
                    key_t tB = tNodes( 1 )->original()->index();

                    key_t tKey = tA > tB ? tA * tNumNodes + tB : tB * tNumNodes + tA;

                    // grab the edge
                    Edge * tEdge = aEdges( mEdgeMap( tKey ) );

                    // insert the edge into the facet
                    tFacet->element()->insert_edge(  tEdge , e );
                }
            }

            // handle midside nodes if they exist
            if ( tNumNodesPerEdge == 3 )
            {

                for ( Facet * tFacet : aFacets )
                {
                    for ( uint e=0; e<tNumEdgesPerFacet; ++e )
                    {
                        tFacet->element()->get_nodes_of_edge( e, tNodes );

                        key_t tA = tNodes( 0 )->original()->index();
                        key_t tB = tNodes( 1 )->original()->index();

                        key_t tKey = tA > tB ? tA * tNumNodes + tB : tB * tNumNodes + tA;

                        // grab the edge
                        Edge * tEdge = aEdges( mEdgeMap( tKey ) );

                        if ( ! tEdge->is_flagged() )
                        {
                            tEdge->insert_node( tNodes( 2 ) , 2 );
                            tEdge->flag();
                        }
                        tFacet->element()->insert_edge( tEdge, e );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_edges_on_layers(
            const Cell< Edge * >  & aEdges,
                  id_t            & aEdgeID,
                  Cell< Layer * > & aLayers )
        {
            uint n = aLayers.size() ;

            for ( uint l=0; l<n; ++l )
            {
                Layer * tLayer = aLayers( l ) ;

                Cell< Node * > & tNodes = tLayer->Nodes ;

                for ( uint d=0; d<2; ++d )
                {
                    Cell < Edge * > & tEdges = d==0 ? tLayer->Edges : tLayer->EdgeDuplicates;
                    index_t tCount = 0 ;
                    tEdges.set_size( aEdges.size(), nullptr );

                    for ( Edge * tOrg : aEdges )
                    {
                        Edge * tDup = new Edge();

                        // populate edges with node copies
                        tDup->allocate_node_container( tOrg->number_of_nodes() );
                        for ( uint k=0; k<tOrg->number_of_nodes(); ++k )
                        {
                            tDup->insert_node( tNodes( tOrg->node( k )->original()->index() ), k );
                        }

                        // set id
                        tDup->set_id( ++aEdgeID );

                        // add edge to container
                        tEdges( tCount++ ) = tDup ;
                    }

                    if ( ! tLayer->hasDuplicates ) break ;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_faces_on_layers(
            const Cell< Facet * >  & aFacets,
                  id_t            & aFaceID,
                  Cell< Layer * > & aLayers )
        {
            uint n = aLayers.size() ;
            for ( uint l=0; l<n; ++l )
            {
                Layer * tLayer = aLayers( l ) ;

                Cell< Node * > & tNodes = tLayer->Nodes ;

                for ( uint d=0; d<2; ++d )
                {
                    Cell < Face * > & tFaces = d == 0 ? tLayer->Faces : tLayer->FaceDuplicates ;

                    tFaces.set_size( aFacets.size(), nullptr );

                    index_t tCount = 0 ;
                    for ( Facet * tOrg : aFacets )
                    {
                        Face * tDup = new Face( tOrg->master(),
                            tOrg->index_on_master(),
                            tOrg->slave(), tOrg->index_on_slave() );

                        tDup->set_id( ++aFaceID );

                        for ( uint k=0; k<tOrg->number_of_nodes(); ++k )
                        {
                            tDup->insert_node( tNodes( tOrg->node( k )->original()->index() ), k );
                        }

                        tFaces( tCount++ ) = tDup ;
                    }

                    if ( ! tLayer->hasDuplicates ) break ;
                }
            }
        }

        void
        ThinShellFactory::link_elements_with_edges(
            Cell< Facet * > & aFacets,
            Cell< Layer * > & aLayers,
            Cell< Block * > & aBlocks )
        {
            // get interpolation order
            uint tOrder = interpolation_order_numeric( aFacets( 0 )->element()->type() );

            if ( tOrder == 1 )
            {
                uint l = 0 ;
                for ( Block * tBlock : aBlocks )
                {
                    Cell< Edge * > & tBottom = aLayers( l )->hasDuplicates ? aLayers( l )->EdgeDuplicates : aLayers( l )->Edges ;
                    Cell< Edge * > & tTop = aLayers( ++l )->Edges ;

                    index_t tElemCount = 0 ;
                    for ( Facet * tFacet : aFacets  )
                    {
                        Element * tElement = tBlock->element( tElemCount++ );

                        tElement->allocate_edge_container();

                        uint tEdgeCount = 0 ;
                        for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                        {
                            tElement->insert_edge( tBottom( tFacet->element()->edge( e )->index()), tEdgeCount++ );
                        }
                        for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                        {
                            tElement->insert_edge( tTop( tFacet->element()->edge( e )->index()), tEdgeCount++ );
                        }
                    }
                }
            }
            else if ( tOrder == 2 )
            {
                uint l = 0 ;
                for ( Block * tBlock : aBlocks )
                {

                    Cell< Edge * > & tBottom = aLayers( l )->hasDuplicates ? aLayers( l )->EdgeDuplicates : aLayers( l )->Edges ;
                    Cell< Edge * > & tMid = aLayers( ++l )->Edges ;
                    Cell< Edge * > & tTop = aLayers( ++l )->Edges ;

                    index_t tElemCount = 0 ;
                    for ( Facet * tFacet : aFacets  )
                    {
                        Element * tElement = tBlock->element( tElemCount++ );

                        tElement->allocate_edge_container();

                        uint tEdgeCount = 0 ;
                        for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                        {
                            tElement->insert_edge( tBottom( tFacet->element()->edge( e )->index()), tEdgeCount++ );
                        }
                        for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                        {
                            tElement->insert_edge( tMid( tFacet->element()->edge( e )->index()), tEdgeCount++ );
                        }
                        for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                        {
                            tElement->insert_edge( tTop( tFacet->element()->edge( e )->index()), tEdgeCount++ );
                        }
                    }
                }
            }
            else
            {
                BELFEM_ERROR( false, "interpolation order not supported");
            }

        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::link_elements_with_faces(
            Cell< Layer * > & aLayers,
            Cell< Block * > & aBlocks )
        {
            uint l = 0 ;
            for ( Block * tBlock : aBlocks )
            {
                Cell< Element * > & tElements = tBlock->elements() ;
                Cell< Face * > & tBottom = aLayers( l )->hasDuplicates ? aLayers( l )->FaceDuplicates : aLayers( l )->Faces ;
                Cell< Face * > & tMid = aLayers( ++l )->Faces ;
                Cell< Face * > & tTop = aLayers( ++l )->Faces ;

                index_t tCount = 0 ;

                for ( Element * tElement : tElements )
                {
                    tElement->allocate_face_container();
                    tElement->insert_face( tBottom( tCount ), 0 );
                    tElement->insert_face( tTop( tCount ), 1 );
                    tElement->insert_face( tMid( tCount ), 2 );
                    ++tCount ;
                }
            }
        }

        void
        ThinShellFactory::create_buffers(
                const Cell< string > & aMaterials,
                Cell< Layer * >      & aLayers,
                Cell< Block * >      & aBlocks )
        {
            index_t n = aBlocks.size() ;


            for ( Layer * tLayer : aLayers )
            {
                for ( Edge * tEdge : tLayer->Edges )
                {
                    tEdge->unflag();
                }
                for ( Face * tFace : tLayer->Faces )
                {
                    tFace->unflag();
                }
            }

            Cell< Node * > tNodes ;

            for ( index_t b=0; b<n; ++b )
            {

                const string tLowerName = string_to_lower( aMaterials( b ) );

                if ( mMaterialMap->find( tLowerName )->second->have( MaterialProperty::rho ) ) continue ;

                Block * tBlock = aBlocks( b ) ;

                tBlock->set_domain_type( DomainType::Buffer ) ;

                for ( Element * tElement : tBlock->elements() )
                {
                    for ( uint e=0; e<tElement->number_of_edges(); ++e )
                    {
                        Edge * tEdge = tElement->edge( e ) ;

                        if ( tEdge->is_flagged() ) continue ;

                        // Allocate storage for source nodes and their weights
                        tEdge->allocate_source_container(
                            tEdge->number_of_nodes() );

                        // For each node on this master edge, find corresponding slave node
                        for ( uint k = 0; k < tEdge->number_of_nodes(); ++k )
                        {
                            // note: we determine the weights later in
                            // DofData::create_dofwise_t_matrices_master()
                            tEdge->add_source( tEdge->node( k ) );
                        }

                        tEdge->flag();
                    }
                }

            }
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_ghost_facets(
                const uint aOrder,
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks )
        {
            index_t l = 0 ;
            index_t nb = aBlocks.size() ;
            index_t nf  = aFacets.size() ;

            ElementFactory tFactory ;

            if ( nb == 1 ) return ;

            for ( index_t b = 1; b < nb; ++b )
            {
                l += aOrder ;
                Layer * tLayer = aLayers( l );

                if ( ! tLayer->hasDuplicates ) continue ;

                Block * tMasterBlock = aBlocks( b-1 ) ;
                Block * tSlaveBlock = aBlocks( b ) ;


                if ( tMasterBlock->domain_type() == DomainType::Buffer || tSlaveBlock->domain_type() == DomainType::Buffer ) continue ;

                Cell< Facet * > & tFacets = tLayer->GhostFacets ;

                tFacets.set_size( aFacets.size(), nullptr );

                for ( index_t f = 0; f < nf; ++f )
                {
                    // original facet for reference
                    Facet * tOrg = aFacets( f ) ;

                    // new element
                    Element * tElement = tFactory.create_element(
                        tOrg->element()->type(),++mMaxElementID );

                    // new facet
                    Facet * tGhost = new Facet( tElement );

                    // Ghost facet sits at the interface between the lower
                    // (master) and upper (slave) shell blocks: it meets the
                    // master's top face and the slave's bottom face.
                    Element * tMasterElement = tMasterBlock->element( f );
                    Element * tSlaveElement  = tSlaveBlock->element( f );

                    tGhost->set_master(
                        tMasterElement,
                        top_facet_index( tMasterElement->type() ) );
                    tGhost->set_slave(
                        tSlaveElement,
                        bottom_facet_index( tSlaveElement->type() ),
                        1 );

                    tFacets( f ) = tGhost ;
                }
            }
        }


//------------------------------------------------------------------------------

        id_t
        ThinShellFactory::max_node_id()
        {
            id_t aID = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->id() > aID )
                {
                    aID = tNode->id();
                }
            }
            return aID ;
        }

//------------------------------------------------------------------------------

        id_t
        ThinShellFactory::max_element_id()
        {
            id_t aID = 0 ;
            for ( Element * tElement : mMesh->elements() )
            {
                if ( tElement->id() > aID )
                {
                    aID = tElement->id();
                }
            }
            for ( Facet * tFacet : mMesh->facets() )
            {
                if ( tFacet->id() > aID )
                {
                    aID = tFacet->id();
                }
            }
            for ( Edge * tEdge : mMesh->edges() )
            {
                if ( tEdge->id() > aID )
                {
                    aID = tEdge->id();
                }
            }
            for ( Face * tFace : mMesh->faces() )
            {
                if ( tFace->id() > aID )
                {
                    aID = tFace->id();
                }
            }
            return aID ;
        }

//------------------------------------------------------------------------------

        id_t
        ThinShellFactory::max_block_id()
        {
            id_t aID = 0 ;
            for ( Block * tBlock : mMesh->blocks() )
            {
                if ( tBlock->id() > aID )
                {
                    aID = tBlock->id();
                }
            }
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( tSideSet->id() > aID )
                {
                    aID = tSideSet->id();
                }
            }
            return aID ;
        }

//------------------------------------------------------------------------------

        real
        ThinShellFactory::compute_binomial_vectors(
            Protoshell * aProtoShell,
            Curve * aCurve,
            const Matrix< real > & aNormals,
            Matrix< real > & aBinomials )
        {
            Cell< Node * > & tCurveNodes = aCurve->nodes() ;

            // the three-point tangent stencil needs at least three nodes
            BELFEM_ERROR( tCurveNodes.size() >= 3,
                "side curve %lu must contain at least three nodes "
                "to compute the binomial vectors",
                ( long unsigned int ) aCurve->id() );

            // note: the tape-side sign is anchored to a terminal curve below.
            // There is no reliable closure test here — is_closed() is
            // overloaded by the CutFactory ( "not shared with another shell" )
            // and the node list may or may not repeat the start node,
            // depending on which pipeline built the curve. The terminal
            // search below fails loudly if no anchor exists, which also
            // covers closed-loop tapes.

            // Help vectors
            Vector< real > P( 3 );  // point 0
            Vector< real > Q( 3 );  // point 1
            Vector< real > R( 3 );  // point 2
            Vector< real > C(3);    // coefficients
            Vector< real > S(3);    // coefficients
            Vector< real > B( 3 );  // binomial vector
            Vector< real > T( 3 );  // tangential vector
            Vector< real > N( 3 );  // normal vector

            // first, we must figure out on which side of the curve the tape is.
            // we must do this over the terminals

            // get the first node of the curve
            Node * tStart = tCurveNodes.first() ;

            // get the normal direction of the tape
            N = aNormals.col( tStart->index() ) ;

            // get the direction vector of the curve
            tStart->get_coords( P ) ;
            tCurveNodes(1)->get_coords( Q );

            // compute the tangent vector
            T = Q - P ;
            T /= norm( T );

            bool tFound = false ;
            for ( Curve * tOther : aProtoShell->terminal_curves() )
            {
                tOther->nodes().last()->get_coords( P );
                for ( Node * tNode : tOther->nodes() )
                {
                    Q = P ;
                    tNode->get_coords( P );
                    if ( tNode == tStart )
                    {
                        tFound = true ;
                        break ;
                    }
                }
                if ( tFound ) break ;
            }
            BELFEM_ERROR( tFound,
                "Side curve %lu has no terminal curve anchor: its first "
                "node is not contained in any terminal curve of the "
                "protoshell. The side-connector sign computation needs "
                "this anchor. Closed-loop tapes have no terminals and "
                "are not supported yet.",
                ( long unsigned int ) aCurve->id() );

            Q -= P ;
            Q/=norm( Q ) ;
            B = cross( T, N );

            // Q points from the curve start into the tape; the sign test
            // must stay valid on oblique corners, so we use the projection
            // rather than a distance threshold
            real tProjection = dot( Q, B );

            BELFEM_ERROR( std::abs( tProjection ) > BELFEM_EPSILON,
                "cannot determine on which side of side curve %lu the tape is: "
                "the terminal curve direction is perpendicular to the binormal",
                ( long unsigned int ) aCurve->id() );

            real aSign = tProjection > 0.0  ? -1.0 :  1.0 ;

            // create additional nodes
            index_t n = aCurve->nodes().size() ;

            Matrix< real > M(3,3);

            index_t p ;
            index_t q ;
            index_t r ;

            aBinomials.set_size( 3, n );

            for ( index_t i=0; i<n; ++i )
            {
                if ( i == 0 )
                {
                    p = 0 ;
                    q = 1 ;
                    r = 2 ;
                }
                else if ( i == n-1 )
                {
                    p = i-2 ;
                    q = i-1 ;
                    r = i ;
                }
                else
                {
                    p = i-1 ;
                    q = i ;
                    r = i+1 ;
                }

                tCurveNodes( p )->get_coords( P );
                tCurveNodes( q )->get_coords( Q );
                tCurveNodes( r )->get_coords( R );

                real a = -norm(Q-P);
                real b =  norm(R-Q);

                // these are the coefficients we need to compute the tangent
                M(0,0) = 1./(a*(a-b));
                M(1,0) = b/(a*(b-a));
                M(2,0) = 0.0 ;
                M(0,1) = 1./(a*b);
                M(1,1) = -(a+b)/(a*b);
                M(2,1) = 1.0 ;
                M(0,2) = 1./(b*(b-a));
                M(1,2) = a/(b*(a-b));
                M(2,2) = 0.0 ;

                real s = ( i==0 ) ? a : ( i==n-1 ) ? b : 0.0;
                for ( uint j=0; j<3; ++j )
                {
                    S(0) = P(j);
                    S(1) = Q(j);
                    S(2) = R(j);
                    C = M * S ;

                    T( j ) = 2.0 * C(0)* s + C(1) ;
                }
                T/=norm(T);

                // now we get the normal
                N = aNormals.col( tCurveNodes( i )->index() ) ;

                // binomial vector
                B = cross( T, N );
                B/=norm(B);

                // tidy up
                for ( uint j=0; j<3; ++j )
                {
                    if ( std::abs( B( j )) < BELFEM_EPSILON )  B( j ) = 0.0 ;
                }
                aBinomials.set_col( i, B );
            }

            return aSign ;
        }

//------------------------------------------------------------------------------

        void
        ThinShellFactory::create_edge_map( const key_t aNumNodes, Cell< Edge * > & aEdges, Map< key_t, index_t > & aEdgeMap )
        {

            // create edge map
            for ( Edge * tEdge : aEdges )
            {
                key_t tA = tEdge->node( 0)->original()->index();
                key_t tB = tEdge->node( 1)->original()->index();
                key_t tKey = tA > tB ? tA * aNumNodes + tB : tB * aNumNodes + tA;
                aEdgeMap[ tKey ] = tEdge->index() ;
            }
        }

        void
        ThinShellFactory::create_edge_to_face_map(
                        Cell< Edge * > & aEdges,
                        Cell< Facet * > & aFacets,
                        const Cell< index_t > & aIndices,
                        Map< index_t, std::pair< index_t, int8_t > > & aMap )
        {
            // STEP 0 : reset the flags
            for ( Facet * tFacet : aFacets )
            {
                tFacet->unflag( 1 );
                for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                {
                    tFacet->edge( e )->unflag( 1 );
                }
            }

            // STEP 1: flag only the side edges of the current curve
            // ( flagging all temporary edges would make every facet edge
            //   a hit, since the facet containers hold nothing else )
            for ( index_t i : aIndices )
            {
                aEdges( i )->flag( 1 );
            }

            // STEP 2: record each side edge with its facet and slot.
            // A healthy mesh has at most one side edge per facet; a second
            // hit on the same facet means the data is corrupted, so we
            // check every slot instead of breaking on the first
            index_t tCount = 0 ;
            for ( Facet * tFacet : aFacets )
            {
                for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                {
                    if ( tFacet->edge( e )->is_flagged( 1 ) )
                    {
                        BELFEM_ERROR( ! tFacet->is_flagged( 1 ),
                            "Facet %lu is already flagged",
                            ( long unsigned int ) tFacet->id() );

                        aMap[ tFacet->edge( e )->index() ] = { tCount, e } ;
                        tFacet->flag( 1 );
                    }
                }
                ++tCount ;
            }

            // cleanup
            for ( Facet * tFacet : aFacets )
            {
                for ( uint e=0; e<tFacet->number_of_edges(); ++e )
                {
                    tFacet->edge( e )->unflag( 1 );
                }
                tFacet->unflag( 1 );
            }
        }
        
        void
        ThinShellFactory::compute_side_edge_indices(
                Curve * aCurve,
                const Map< key_t, index_t > & aEdgeMap,
                const key_t                   aNumNodes ,
                const uint                    aEdgeOrder,
                Cell< index_t >             & aEdgeIndices )
        {
            // number of nodes
            index_t n = aCurve->nodes().size();

            // number of edges
            index_t m = ( n - 1 )/aEdgeOrder ;

            BELFEM_ERROR( n >= aEdgeOrder+1,
                "side curve %lu must contain at least %u nodes",
                ( long unsigned int ) aCurve->id(),
                ( unsigned int ) aEdgeOrder + 1 );

            aEdgeIndices.set_size( m, 0 );

            Cell< Node * > & tNodes = aCurve->nodes();

            // keys must be original-normalized: the edge map is built from
            // original()->index(), and a side-curve station can be a cut
            // duplicate
            index_t off = 0 ;
            Node * tPrev = tNodes( off );
            key_t p = tPrev->original()->index();
            for ( index_t e=0; e<m; ++e )
            {
                off += aEdgeOrder ;
                key_t q = p ;
                Node * tThis = tNodes( off );
                p = tThis->original()->index() ;
                key_t tKey = p > q ? p * aNumNodes + q : q * aNumNodes + p ;

                // a station pair that is not a facet edge of the shell means
                // the curve walk and the edge map disagree. Report the station
                // instead of letting Map::operator() abort with a bare
                // "Key not found in map"
                BELFEM_ERROR( aEdgeMap.key_exists( tKey ),
                    "side curve %lu, station %lu of %lu: node pair is not an edge of the thin shell\n"
                    "    node A: id %lu original %lu index %lu at ( %g, %g, %g )\n"
                    "    node B: id %lu original %lu index %lu at ( %g, %g, %g )\n"
                    "    curve nodes %lu, edge order %u, radix %lu, map entries %lu",
                    ( long unsigned int ) aCurve->id(),
                    ( long unsigned int ) e,
                    ( long unsigned int ) m,
                    ( long unsigned int ) tPrev->id(),
                    ( long unsigned int ) tPrev->original()->id(),
                    ( long unsigned int ) tPrev->original()->index(),
                    tPrev->x(), tPrev->y(), tPrev->z(),
                    ( long unsigned int ) tThis->id(),
                    ( long unsigned int ) tThis->original()->id(),
                    ( long unsigned int ) tThis->original()->index(),
                    tThis->x(), tThis->y(), tThis->z(),
                    ( long unsigned int ) n,
                    ( unsigned int ) aEdgeOrder,
                    ( long unsigned int ) aNumNodes,
                    ( long unsigned int ) aEdgeMap.size() );

                aEdgeIndices( e ) = aEdgeMap( tKey );
                tPrev = tThis ;
            }
        }

        void
        ThinShellFactory::compute_side_authority(
                Curve * aCurve,
                Cell< Facet * > & aFacets,
                Cell< Edge * >  & aEdges,
                const Cell< index_t > & aEdgeIndices,
                Cell< Node * > & aEdgeSources,
                Cell< Node * > & aNodeSources )
        {
            // Resolve each rim station against the air volume elements that
            // hang_thinshell_edges_on_nodes_bottom/top read later. The volume
            // element->node links are already cut-composed because the
            // CutFactory runs before this factory; matching by original()
            // identity keeps this free of any orientation convention.

            index_t tNumEdges = aEdgeIndices.size() ;

            Cell< Node * > & tCurveNodes = aCurve->nodes() ;
            index_t tNumStations = tCurveNodes.size() ;

            aEdgeSources.set_size( 2 * tNumEdges, nullptr );
            aNodeSources.set_size( tNumStations, nullptr );

            // two curve positions resolving to the same temporary edge means
            // the original-key unique() in create_temporary_edges collapsed
            // twin edges ( cut sides or periodic partners ); the collapsed
            // location then keeps one branch only
            {
                Map< index_t, index_t > tSeen ;
                index_t tNumCollisions = 0 ;
                for ( index_t e : aEdgeIndices )
                {
                    if ( tSeen.key_exists( e ) )
                    {
                        ++tNumCollisions ;
                    }
                    else
                    {
                        tSeen[ e ] = 0 ;
                    }
                }
                if ( tNumCollisions > 0 )
                {
                    message( InfoLevel::Default,
                        "    Warning: side curve %lu: %lu edge positions collapse onto twins",
                        ( long unsigned int ) aCurve->id(),
                        ( long unsigned int ) tNumCollisions );
                }
            }

            // find the facet that carries each side edge
            Map< index_t, std::pair< index_t, int8_t > > tEdgeToFacetMap ;
            this->create_edge_to_face_map( aEdges, aFacets, aEdgeIndices, tEdgeToFacetMap );

            // one authority per station ( original index ), first writer wins
            Map< index_t, Node * > tStationAuthority ;

            Cell< Node * > tNodesOnVolume ;

            for ( index_t p=0; p<tNumEdges; ++p )
            {
                Edge * tOrg = aEdges( aEdgeIndices( p ) );

                BELFEM_ASSERT( tOrg->number_of_nodes() == 2,
                    "side edge authority is only implemented for first order edges" );

                Facet * tFacet = aFacets(
                    tEdgeToFacetMap( aEdgeIndices( p ) ).first );

                // authority side: master, unless it is a conductor
                // ( mirrors the air-master routing of the anchor pass )
                bool tUseMaster = mMesh->block( tFacet->master()->block_id() )
                    ->domain_type() != DomainType::Conductor ;

                if ( tUseMaster )
                {
                    tFacet->master()->get_nodes_of_facet(
                        tFacet->index_on_master(), tNodesOnVolume );
                }
                else
                {
                    BELFEM_ERROR( tFacet->slave() != nullptr &&
                        mMesh->block( tFacet->slave()->block_id() )
                            ->domain_type() != DomainType::Conductor,
                        "no phi authority for side edge fusing: both sides of facet %lu are conductors",
                        ( long unsigned int ) tFacet->id() );

                    tFacet->slave()->get_nodes_of_facet(
                        tFacet->index_on_slave(), tNodesOnVolume );
                }

                for ( uint s=0; s<2; ++s )
                {
                    index_t tStation = tOrg->node( s )->original()->index() ;

                    Node * tHit  = nullptr ;
                    uint   tHits = 0 ;

                    for ( Node * tNode : tNodesOnVolume )
                    {
                        if ( tNode->original()->index() == tStation )
                        {
                            tHit = tNode ;
                            ++tHits ;
                        }
                    }

                    // an unlinked duplicate on the volume side would land here
                    BELFEM_ERROR( tHits == 1,
                        "side edge fusing: station node %lu resolves %u times on facet %lu",
                        ( long unsigned int ) tOrg->node( s )->id(),
                        ( unsigned int ) tHits,
                        ( long unsigned int ) tFacet->id() );

                    aEdgeSources( 2 * p + s ) = tHit ;

                    if ( ! tStationAuthority.key_exists( tStation ) )
                    {
                        tStationAuthority[ tStation ] = tHit ;
                    }
                }
            }

            // hand the station authorities out in curve order
            for ( index_t k=0; k<tNumStations; ++k )
            {
                index_t tStation = tCurveNodes( k )->original()->index() ;

                BELFEM_ERROR( tStationAuthority.key_exists( tStation ),
                    "side edge fusing: no authority for station node %lu on curve %lu",
                    ( long unsigned int ) tCurveNodes( k )->id(),
                    ( long unsigned int ) aCurve->id() );

                aNodeSources( k ) = tStationAuthority( tStation );
            }
        }

        void
        ThinShellFactory::connect_side_edges(
            Layer * aLayer,
            Cell< Edge * > & aEdges,
            const Cell< index_t > & aEdgeIndices,
            const Cell< Node * > & aEdgeSources )
        {

            if ( aEdges.size() == 0 ) return ;

            uint n = aEdges.first()->number_of_nodes();

            BELFEM_ASSERT( aEdgeSources.size() == n * aEdgeIndices.size(),
                "size mismatch of edge authority container" );

            for ( uint d=0; d<2; ++d )
            {
                Cell< Edge * > & tEdges = d == 0 ? aLayer->Edges : aLayer->EdgeDuplicates ;

                index_t p = 0 ;
                for ( index_t e: aEdgeIndices )
                {
                    Edge * tDup = tEdges( e );

                    if ( tDup->is_hanging() ) tDup->reset_source_container() ;

                    tDup->allocate_source_container( n );
                    for ( uint s=0; s<n; ++s )
                    {
                        // authority source: same station, but the cut-composed
                        // branch the outer-interface anchors resolve to
                        tDup->add_source( aEdgeSources( p * n + s ) );
                    }
                    ++p ;
                }

                if ( ! aLayer->hasDuplicates ) break ;
            }
        }

        void
        ThinShellFactory::connect_side_edges(
            SideLayer * aSideLayer,
            Cell< Edge * > & aEdges,
            const Cell< index_t > & aEdgeIndices )
        {

            if ( aEdges.size() == 0 ) return ;

            uint n = aEdges.first()->number_of_nodes();

            Cell< Node * > tNodes( n, nullptr );
            for ( uint d=0; d<2; ++d )
            {
                Cell< Edge * > & tOuterEdges = d == 0 ? aSideLayer->OuterEdges : aSideLayer->OuterEdgeDuplicates ;

                index_t tCount = 0 ;

                // e is the position in the temporary edge container;
                // the layer edge copies carry no index of their own
                for ( index_t e : aEdgeIndices )
                {
                    Edge * tOrg = aEdges( e );
                    Edge * tDup = tOuterEdges( tCount++ );

                    if ( tDup->is_hanging() ) tDup->reset_source_container() ;

                    tDup->allocate_source_container( n );
                    for ( uint s=0; s<n; ++s )
                    {
                        tDup->add_source( tOrg->node( s ) );
                    }
                }

                if ( aSideLayer->OuterEdgeDuplicates.size() == 0) break ;
            }
        }

        void
        ThinShellFactory::connect_side_nodes(
               Cell< Node * > & aSourceNodes,
               const Cell< Node * > & aAuthorityNodes,
               Cell< Node * > & aTargetNodes )
        {
            index_t n = aSourceNodes.size() ;

            BELFEM_ASSERT( aAuthorityNodes.size() == n,
                "size mismatch of node authority container" );

            for ( index_t k=0; k<n; ++k )
            {
                // the authority is the air-side representative of this station,
                // cut-composed like the outer-interface anchors; duplicate
                // stations share one authority ( cf. compute_side_authority )
                Node * tOrg = aAuthorityNodes( k );

                // the target container is positional over the ORIGINAL surface
                // nodes, and a side-curve station can be a cut duplicate whose
                // own index is gNoIndex ( cf. compute_side_edge_indices )
                Node * tDup = aTargetNodes( aSourceNodes( k )->original()->index() );
                if ( tDup->is_hanging() ) tDup->reset_source_container() ;

                if ( tOrg->is_hanging() )
                {
                    tDup->allocate_source_container( tOrg->number_of_sources() );
                    for ( uint s=0; s<tOrg->number_of_sources(); ++s )
                    {
                        tDup->add_source( tOrg->source( s ), tOrg->weight( s ) );
                    }
                }
                else
                {
                    tDup->allocate_source_container( 1 );
                    tDup->add_source( tOrg );
                }
            }
        }

        void
        ThinShellFactory::update_node_edge_tables( Cell< Node * > & aNodes, Cell< Edge * > & aEdges )
        {
            index_t tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->reset_edge_container();
                BELFEM_ASSERT( ! tNode->is_duplicate(), "Node %lu is a duplicate of %lu. All notes must be originals",
                ( long unsigned int ) tNode->id(), ( long unsigned int ) tNode->original()->id() );
            }

            tCount = 0 ;
            for ( Edge * tEdge : aEdges )
            {
                tEdge->set_index( tCount++ );
                Node * tA = tEdge->node( 0 )->original() ;
                Node * tB = tEdge->node( 1 )->original() ;
                tA->increment_edge_counter();
                tB->increment_edge_counter();
            }

            for ( Node * tNode : aNodes )
            {
                tNode->allocate_edge_container();
            }

            for ( Edge * tEdge : aEdges )
            {
                Node * tA = tEdge->node( 0 )->original() ;
                Node * tB = tEdge->node( 1 )->original() ;
                tA->add_edge( tEdge );
                tB->add_edge( tEdge );
            }
        }



        void
        ThinShellFactory::update_node_facet_tables( Cell< Node * > & aNodes, Cell< Facet * > & aFacets )
        {
            Vector< index_t > tIndices( aNodes.size() );

            index_t  tCount = 0 ;

            // backup original indices
            for ( Node * tNode : aNodes )
            {
                tIndices( tCount++ ) = tNode->index() ;
                Node * tOrg = tNode->original() ;
                tOrg->reset_facet_container();
                tOrg->unflag( 2 );
                for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                {
                    Node * tDup = tOrg->duplicate( d ) ;
                    tDup->reset_facet_container();
                    tDup->unflag( 2 );
                }
            }

            for ( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    Node * tOrg = tFacet->node( k )->original();
                    tOrg->increment_facet_counter();
                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        tOrg->duplicate( d )->increment_facet_counter();
                    }
                }
            }

            for ( Node * tNode : aNodes )
            {
                Node * tOrg = tNode->original() ;

                if ( ! tOrg->is_flagged( 2 ) )
                {
                    tOrg->allocate_facet_container();
                    tOrg->flag( 2 );
                }
                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    Node * tDup = tNode->duplicate( d ) ;
                    if ( ! tDup->is_flagged( 2 ) )
                    {
                        tDup->allocate_facet_container();
                        tDup->flag( 2 );
                    }
                }
            }

            for ( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    Node * tOrg = tFacet->node( k )->original();
                    tOrg->add_facet( tFacet );

                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        tOrg->duplicate( d )->add_facet( tFacet );
                    }
                }
            }

            // restore original indices and flag states
            tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tIndices( tCount++ ) );
                Node * tOrg = tNode->original() ;
                tOrg->unflag( 2 );
                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    tNode->duplicate( d )->unflag( 2 );
                }
            }
        }


        void
        ThinShellFactory::flag_periodic_nodes()
        {
            if ( ! mMesh->has_periodicity() ) return ;

            this->unflag_periodic_nodes();
            Cell< Node * > & tMasterNodes = mMesh->periodicity()->master_nodes() ;
            for ( Node * tNode : tMasterNodes )
            {
                tNode->flag( 4 );
            }
            Cell< Node * > & tSlaveNodes = mMesh->periodicity()->slave_nodes() ;
            for ( Node * tNode : tSlaveNodes )
            {
                tNode->flag( 5 );
            }
        }

        void
        ThinShellFactory::unflag_periodic_nodes()
        {
            mMesh->unflag_all_nodes( 4 );
            mMesh->unflag_all_nodes( 5 );
        }


        void
        ThinShellFactory::flag_layer_nodes( Cell< Node * > & aNodes, Cell< Layer * > & aLayers )
        {
            index_t n = aNodes.size() ;

            for ( Layer * tLayer : aLayers )
            {
                Cell< Node * > & tNodes = tLayer->Nodes ;

                for ( index_t k=0; k<n; ++k )
                {
                    if ( aNodes( k )->is_flagged( 4 ) )
                    {
                        tNodes( k )->flag( 4 );
                    }
                    if ( aNodes( k )->is_flagged( 5 ) )
                    {
                        tNodes( k )->flag( 5 );
                    }
                }
            }
        }

        SideSet *
        ThinShellFactory::create_periodic_sideset(  Cell< Block * > & aBlocks, const bool aMaster )
        {
            BELFEM_ERROR( mMesh->max_element_order() == 1, "only linear elements supported for periodic sidesets" );
            uint tFlag = aMaster ? 4 : 5 ;


            id_t tGroup = ++mMaxGroupID ;

            Cell< Facet * > tFacets ;

            ElementFactory tFactory ;

            Cell< Node * > tNodes ;
            for ( Block * tBlock : aBlocks )
            {
                BELFEM_ASSERT( geometry_type( tBlock->element_type() ) == GeometryType::PENTA, "Expect Penta Element" );

                for ( Element * tElement : tBlock->elements() )
                {
                    for ( uint f=0; f<3; ++f )
                    {
                        tElement->get_nodes_of_facet( f, tNodes );

                        bool tIsPeriodic = true ;

                        for ( Node * tNode : tNodes )
                        {
                            if ( ! tNode->is_flagged( tFlag ) )
                            {
                                tIsPeriodic = false ;
                                break ;
                            }
                        }

                        if ( ! tIsPeriodic ) continue ;

                        ElementType tType = element_type_of_facet( tElement->type(), f );

                        Facet * tFacet = new Facet( tFactory.create_element( tType, ++mMaxElementID ) );

                        tFacet->set_master( tElement, f );
                        tFacet->set_sideset_id( tGroup );
                        tFacets.push( tFacet );
                    }
                }
            }

            if ( tFacets.size() == 0 )
            {
                --mMaxGroupID ;
                return nullptr ;
            }

            tFacets.shrink_to_fit();
            SideSet * aSideSet = new SideSet( tGroup, tFacets.size() );

            for ( Facet * tFacet : tFacets )
            {
                aSideSet->insert_facet( tFacet );
            }

            aSideSet->set_domain_type( DomainType::Periodic );
            aSideSet->hide( true );
            return aSideSet;
        }

    }
}
