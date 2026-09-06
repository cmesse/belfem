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

#ifndef CL_FEM_THINSHELLFACTORY_HPP
#define CL_FEM_THINSHELLFACTORY_HPP

#include "cl_Mesh.hpp"
#include "cl_ThinShell.hpp"
#include "cl_Input_Section.hpp"
#include "cl_Protoshell.hpp"
#include "cl_Material.hpp"

namespace belfem
{
    namespace mesh
    {
        class ThinShellFactory
        {
            const proc_t            mCommRank ;
            const uint              mNumDimensions ;
            Mesh *                  mMesh ;
            Map< string, Material * > * mMaterialMap ;
            id_t                    mMaxGroupID ;
            id_t                    mMaxElementID ;
            id_t                    mMaxNodeID ;

            Cell< Node * >    & mMasterNodes ;
            Cell< Node * >    & mSlaveNodes ;
            Cell< index_t >     mNodeIndices ;

            struct Layer
            {
                bool hasDuplicates = false ;
                Cell< Node * >  Nodes ;
                Cell< Edge * >  Edges ;
                Cell< Edge * >  EdgeDuplicates ;
                Cell< Face * >  Faces ;
                Cell< Face * >  FaceDuplicates ;
                Cell< Facet * > GhostFacets ;
            };

            struct SideLayer
            {
                const real Sign ;
                Cell< Node * >    InnerNodes ;
                Cell< Node * >    InnerNodeReference ;
                Cell< Node * >    OuterNodes ;
                Cell< Edge * >    InnerEdges ;
                Cell< Edge * >    InnerEdgeReference ;
                Cell< Edge * >    OuterEdges ;
                Cell< Edge * >    InnerEdgeDuplicates ;
                Cell< Edge * >    InnerEdgeDuplicateReference ;
                Cell< Edge * >    OuterEdgeDuplicates ;

                SideLayer(
                    const real aSign,
                    Curve * aCurve,
                    Layer * aLayer,
                    Cell< Edge * > & aReferenceEdges,
                    const Cell< index_t > & aEdgeIndices,
                    const Matrix< real >  & aBinomials,
                    const real aCoatingThickness,
                          id_t & aMaxNodeID,
                          id_t & aMaxEdgeID,
                    const bool aFuseSideEdges )
                    : Sign ( aSign )
                {
                    // collect the inner edges and their UNIQUE nodes: shared
                    // edge endpoints must appear exactly once, otherwise the
                    // slot containers get null holes and the index backup
                    // stores an already-overwritten temp index on the second
                    // occurrence
                    uint tNumNodesPerEdge = aLayer->Edges.first()->number_of_nodes() ;
                    InnerNodeReference.reserve( aEdgeIndices.size() * tNumNodesPerEdge );
                    InnerEdgeReference.reserve( aEdgeIndices.size() );
                    for ( index_t e : aEdgeIndices )
                    {
                        aLayer->Edges( e )->unflag_nodes( 7 ); // flag 7 is arbitrary, just want to make sure
                                                                     // we don't collide with other flags
                    }
                    for ( index_t e : aEdgeIndices )
                    {
                        Edge * tEdge = aLayer->Edges( e );
                        InnerEdgeReference.push( tEdge );
                        for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                        {
                            Node * tNode = tEdge->node( k );
                            if ( ! tNode->is_flagged( 7 ) )
                            {
                                tNode->flag( 7 );
                                InnerNodeReference.push( tNode );
                            }
                        }
                    }
                    for ( index_t e : aEdgeIndices )
                    {
                        aLayer->Edges( e )->unflag_nodes( 7 );
                    }

                    // collect inner edge duplicates if they exist
                    if ( aLayer->hasDuplicates )
                    {
                        InnerEdgeDuplicateReference.reserve( aEdgeIndices.size() );
                        for ( index_t e : aEdgeIndices )
                        {
                            InnerEdgeDuplicateReference.push( aLayer->EdgeDuplicates( e ) );
                        }
                    }

                    // backup node indices and unflag
                    index_t tCount = 0 ;
                    Cell< index_t > tNodeIndices( InnerNodeReference.size(), 0 );
                    for ( Node * tNode : InnerNodeReference )
                    {
                        tNodeIndices( tCount ) = tNode->index();
                        tNode->set_index( tCount++ );
                    }

                    // create the node duplicates
                    InnerNodes.set_size( InnerNodeReference.size(), nullptr );
                    OuterNodes.set_size( InnerNodeReference.size(), nullptr );

                    index_t tEdgeOrder = InnerEdgeReference.first()->number_of_nodes() - 1 ;

                    Vector< real > tX( 3 );
                    Vector< real > tB( 3 );

                    tCount = 0 ;
                    index_t off = 0 ;

                    // temporary work container for node duplicates
                    Cell< Node * > tSwap ;

                    for ( Edge * tEdge : InnerEdgeReference )
                    {
                        Edge * tRef = aReferenceEdges( aEdgeIndices( tCount++ ) );

                        bool tIsReversed = tRef->node( 1 )->original() == aCurve->nodes()( off )->original();

                        BELFEM_ASSERT( tIsReversed || tRef->node( 0 )->original() == aCurve->nodes()( off )->original(),
                            "Can't determine edge orientation along curve %lu",
                            ( long unsigned int ) aCurve->id() );

                        for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                        {
                            Node * tNode = tEdge->node( k );

                            if ( OuterNodes( tNode->index() ) == nullptr )
                            {
                                switch ( k )
                                {
                                    case 0 :
                                    {
                                        tB = tIsReversed ? aBinomials.col( off + tEdgeOrder ) : aBinomials.col( off );
                                        break ;
                                    }
                                    case 1 :
                                    {
                                        tB = tIsReversed ? aBinomials.col( off ) : aBinomials.col( off + tEdgeOrder );
                                        break ;
                                    }
                                    case 2 :
                                    {
                                        tB = aBinomials.col(  off + 1 );
                                        break ;
                                    }
                                    default:
                                    {
                                        BELFEM_ERROR( false, "not implemented" );
                                    }
                                }

                                tNode->get_coords( tX );

                                Node * tDupI = new Node( ++aMaxNodeID, tX( 0 ), tX( 1 ), tX( 2 ) );
                                tX += aSign * tB * aCoatingThickness ;

                                Node * tOrg0 = tNode->original();

                                uint n = tOrg0->number_of_duplicates();
                                if ( n > 0 )
                                {
                                    tSwap.set_size( n, nullptr );
                                    for ( uint d=0; d<n; ++d )
                                    {
                                        tSwap( d ) = tOrg0->duplicate( d );
                                    }
                                    tOrg0->reset_duplicate_container();
                                    tOrg0->allocate_duplicate_container( n+1 );
                                    for ( Node * tDup : tSwap )
                                    {
                                        tOrg0->add_duplicate( tDup );
                                    }
                                }
                                else
                                {
                                    tOrg0->allocate_duplicate_container( 1 );
                                }
                                tOrg0->add_duplicate( tDupI );
                                tDupI->set_original( tOrg0 );

                                Node * tDupO = new Node( ++aMaxNodeID, tX( 0 ), tX( 1 ), tX( 2 ) );


                                Node * tOrg = tRef->node( k );

                                if ( tOrg->is_hanging() )
                                {
                                    tDupI->allocate_source_container( tOrg->number_of_sources() );
                                    tDupO->allocate_source_container( tOrg->number_of_sources() );
                                    for ( uint s=0; s<tOrg->number_of_sources(); ++s )
                                    {
                                        tDupI->add_source( tOrg->source( s ), tOrg->weight( s ) );
                                        tDupO->add_source( tOrg->source( s ), tOrg->weight( s ) );
                                    }
                                }
                                else
                                {
                                    tDupO->allocate_source_container( 1 );
                                    tDupO->add_source( tOrg );
                                    tDupI->allocate_source_container( 1 );
                                    tDupI->add_source( tOrg );
                                }
                                InnerNodes( tNode->index() ) = tDupI ;
                                OuterNodes( tNode->index() ) = tDupO ;
                            }
                        }

                        off += tEdgeOrder ;
                    }

                    // create the decoupled edge copies. The references are the
                    // SOURCE and must never receive pushes ( self-pushing a
                    // range-iterated Cell reallocates and dangles the loop )
                    for ( uint d=0; d<2; ++d )
                    {
                        Cell< Edge * > & tReferenceEdges = d==0 ? InnerEdgeReference : InnerEdgeDuplicateReference ;
                        Cell< Edge * > & tInnerEdges = d==0 ? InnerEdges : InnerEdgeDuplicates ;
                        Cell< Edge * > & tOuterEdges = d==0 ? OuterEdges : OuterEdgeDuplicates ;

                        tInnerEdges.reserve( tReferenceEdges.size() );
                        tOuterEdges.reserve( tReferenceEdges.size() );

                        tCount = 0 ;
                        for ( Edge * tEdge : tReferenceEdges )
                        {
                            Edge * tDupI = new Edge();
                            tDupI->set_id( ++aMaxEdgeID );
                            tDupI->allocate_node_container( tEdge->number_of_nodes() );
                            for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                            {
                                // duplicate-sheet edges ( d == 1 ) reference the
                                // duplicate-sheet nodes, which were not in the
                                // gather: resolve the slot through their
                                // per-layer original ( the primary-sheet node ).
                                // Both twin edges share the wall node copies
                                Node * tSlotNode = d == 0 ?
                                    tEdge->node( k ) : tEdge->node( k )->original() ;
                                tDupI->insert_node( InnerNodes( tSlotNode->index() ), k );
                            }
                            if ( tEdge->is_hanging() )
                            {
                                tDupI->allocate_source_container( tEdge->number_of_sources() );
                                for ( uint s=0; s<tEdge->number_of_sources(); ++s )
                                {
                                    tDupI->add_source( tEdge->source( s ), tEdge->weight( s ) );
                                }
                            }
                            else
                            {
                                tDupI->allocate_source_container( 1 );
                                tDupI->add_source( tEdge, 1.0 );
                            }

                            tInnerEdges.push( tDupI );

                            Edge * tDupO = new Edge();
                            tDupO->set_id( ++aMaxEdgeID );
                            tDupO->allocate_node_container( tEdge->number_of_nodes() );
                            for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                            {
                                // slot resolution as for the inner copy above
                                Node * tSlotNode = d == 0 ?
                                    tEdge->node( k ) : tEdge->node( k )->original() ;
                                tDupO->insert_node( OuterNodes( tSlotNode->index() ), k );
                            }
                            tOuterEdges.push( tDupO );

                            // fusing the outer side edges
                            if( aFuseSideEdges )
                            {
                                // layer edge copies carry no index of their own
                                Edge * tOrg = aReferenceEdges( aEdgeIndices( tCount++ ) );
                                tDupO->allocate_source_container( tOrg->number_of_nodes() );
                                for ( uint k=0; k<tOrg->number_of_nodes(); ++k )
                                {
                                    // weights are computed later in
                                    // DofData::create_dofwise_t_matrices_master()
                                    tDupO->add_source( tOrg->node( k ) );
                                }
                            }
                        }

                        if ( ! aLayer->hasDuplicates ) break ;
                    }

                    // restore node indices; outer nodes mirror the
                    // shell-local index of their inner partner
                    tCount = 0 ;
                    for ( Node * tNode : InnerNodeReference )
                    {
                        tNode->set_index( tNodeIndices( tCount ) );
                        InnerNodes( tCount )->set_index( tNodeIndices( tCount ) );
                        OuterNodes( tCount )->set_index( tNodeIndices( tCount ) );
                        ++tCount ;
                    }

                }

                ~SideLayer() = default ;
            };

            //! edge-coating wall width: 0 = derive from the outer stabilizer
            //! layer thickness ( side coating deposits in the same plating
            //! step ); overridden per tape by "edge coating width" in the input
            real mConnectorWidth = 0.0 ;
            //! set from the constructor argument; see there
            bool mCreateGhostFacets = false ;
            //! set per tape in create() from the input key "edge coating : on"
            bool mCreateSideConnectors = false ;
            bool mConnectorsForAllLayers = true ;

            Map< key_t , index_t > mEdgeMap ;

            // experimental switch that fuses the side edges to the air trace.
            // off: free rims ( validated legacy ). The fuse is branch-coherent
            // since the single-authority fix, but suppresses the through-
            // thickness branch transition where the upper and lower cuts
            // differ — see todo/side_edge_fusing_cut_aware_plan.md, O1.
            // Physics position ( C. Messe / Prof. Sirous, 2026-08-13 ):
            // fusing is the mathematically cleaner continuity statement, yet
            // measured runs converge slower for no better result — the fuse
            // appears to overconstrain the rim, like enforcing a B·n = 0
            // that the formulation already fulfills at the boundary. Both
            // flags therefore stay false
            bool mFuseEdges = false ;
            bool mFuseEdgesWhenHavingSideConnectors = false ;
//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            /**
             * @param aMesh               the mesh the shell layers are built into
             * @param aMasterNodes        the original ( non-duplicate ) nodes of the tape
             *        surfaces, as collected by CutFactory::thin_shell_master_nodes()
             * @param aSlaveNodes         per master node, its cut duplicate, or the master
             *        node itself where no duplicate exists ( same ordering )
             * @param aMaterialMap        lower-case label -> material map the layer
             *        materials of a protoshell are resolved against
             * @param aCreateGhostFacets  duplicate the edge / face dofs at every
             *        interface between two different resistive layers and
             *        create the ghost facets that couple them ( the Nitsche
             *        ghost, mt_maxwell_h.cpp ). false = the layers share their
             *        interface entities. Decided by the deck
             *        ( fn_FEM_ghost_switch.hpp ), off by default
             */
            ThinShellFactory(
                Mesh * aMesh,
                Cell< Node * > & aMasterNodes,
                Cell< Node * > & aSlaveNodes,
                Map< string, Material * > * aMaterialMap = nullptr,
                const bool aCreateGhostFacets = false );

            ~ThinShellFactory();

            ThinShell *
            create( Protoshell * aProtoShell );

            void
            flag_periodic_nodes();

            void
            unflag_periodic_nodes();

//------------------------------------------------------------------------------
            private :
//------------------------------------------------------------------------------

            void
            create_side_connectors(
                const uint aOrder,
                Protoshell * aProtoShell,
                const Cell< id_t > & aBlockIDs,
                Cell< Layer * > & aLayers,
                Matrix< real > & aNodeNormals,
                Cell< Edge* > & aEdges,
                Cell< Facet * > & aFacets,
                ThinShell * aThinShell );

//------------------------------------------------------------------------------

            void
            reset_node_indices();

//------------------------------------------------------------------------------

            void
            write_normals_to_mesh( Cell< Node * > & aNodes, const Matrix< real > & aNodeNormals );

//------------------------------------------------------------------------------

            void
            collect_sidesets( const Vector< id_t > & aSideSetIDs, Cell< SideSet * > & aSideSets );

//------------------------------------------------------------------------------

            ElementType
            collect_facets( Cell< SideSet * > & aSideSets, Cell< Facet * > & aFacets );

//------------------------------------------------------------------------------

            void
            collect_nodes( Cell< Facet * > & aFacets, Cell< Node * > & aNodes );

//------------------------------------------------------------------------------
// Begin new Side connector functions
//------------------------------------------------------------------------------

            void
            create_edge_map(
                const key_t aNumNodes,
                Cell< Edge * > & aEdges,
                Map< key_t, index_t > & aEdgeMap );

            void
            create_edge_to_face_map(
                Cell< Edge * > & aEdges,
                Cell< Facet * > & aFacets,
                const Cell< index_t > & aIndices,
                Map< index_t, std::pair< index_t, int8_t > > & aMap );

            void
            compute_side_edge_indices(
                Curve * aCurve,
                const Map< key_t, index_t > & aEdgeMap,
                const key_t                 aNumNodes ,
                const uint                    aEdgeOrder,
                Cell< index_t >             & aEdgeIndices );

            /**
             * resolve each side-curve station ONCE against the air volume
             * elements that the outer-interface anchors read later, so that
             * every interior level and the node fusing tie to the same
             * cut-composed branch of the air trace g
             *
             * @param aCurve        side curve, stations in traversal order
             * @param aFacets       mid-surface facets
             * @param aEdges        temporary side edges
             * @param aEdgeIndices  curve position -> temporary edge index
             * @param aEdgeSources  out: source pair per curve position,
             *                      flat ( 2*p+s ), ordered like the
             *                      temporary edge's own nodes
             * @param aNodeSources  out: one authority node per curve station
             */
            void
            compute_side_authority(
                Curve * aCurve,
                Cell< Facet * > & aFacets,
                Cell< Edge * >  & aEdges,
                const Cell< index_t > & aEdgeIndices,
                Cell< Node * > & aEdgeSources,
                Cell< Node * > & aNodeSources );

            void
            connect_side_edges(
                Layer * aLayer,
                Cell< Edge * > & aEdges,
                const Cell< index_t > & aEdgeIndices,
                const Cell< Node * > & aEdgeSources );

            void
            connect_side_edges(
                SideLayer * aSideLayer,
                Cell< Edge * > & aEdges,
                const Cell< index_t > & aEdgeIndices );

            void
            connect_side_nodes(
               Cell< Node * > & aSourceNodes,
               const Cell< Node * > & aAuthorityNodes,
               Cell< Node * > & aTargetNodes );

//------------------------------------------------------------------------------
// End new Side connector functions
//------------------------------------------------------------------------------

            void
            process_nodes_line2(
                const Cell< Facet * > & aFacets,
                const Cell< Node * > & aNodes,
                      Matrix< real > & aNodeNormals );

//------------------------------------------------------------------------------

            void
            process_nodes_line3(
                const Cell< Facet * > & aFacets,
                const Cell< Node * > & aNodes,
                      Matrix< real > & aNodeNormals );

//------------------------------------------------------------------------------

            void
            process_nodes_tri3(
                const Cell< Facet * > & aFacets,
                const Cell< Node * > & aNodes,
                      Matrix< real > & aNodeNormals );

//------------------------------------------------------------------------------

            void
            process_nodes_tri6(
                const Cell< Facet * > & aFacets,
                const Cell< Node * > & aNodes,
                      Matrix< real > & aNodeNormals );

//------------------------------------------------------------------------------

            void
            compute_distances(
                const uint aOrder,
                const Vector< real > & aThicknesses,
                Vector< real > & aDistances );

//------------------------------------------------------------------------------

            void
            create_nodes_on_layers(
                      Cell< Node * >  & aNodes,
                const Matrix< real >  & aNodeNormals,
                const Vector< real >  & aDistances,
                      Cell< Layer * > & aLayers );

//------------------------------------------------------------------------------

            void
            create_edges_on_layers(
                const Cell< Edge * >  & aEdges,
                      id_t            & aEdgeID,
                      Cell< Layer * > & aLayers );

//------------------------------------------------------------------------------

            void
            create_faces_on_layers(
                const Cell< Facet * >  & aFacets,
                      id_t            & aFaceID,
                      Cell< Layer * > & aLayers );

//------------------------------------------------------------------------------

            void
            create_temporary_edges(
               Cell< Node * >   & aNodes,
               Cell< Facet * >  & aFacets,
               Cell< Edge* >    & aEdges );

//------------------------------------------------------------------------------

            void
            create_elements_on_blocks_line2(
                id_t & aBlockID,
                id_t & aElementID,
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks ) ;

//------------------------------------------------------------------------------

            void
            create_elements_on_blocks_line3(
                id_t & aBlockID,
                id_t & aElementID,
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks  ) ;

//------------------------------------------------------------------------------

            void
            create_elements_on_blocks_tri3(
                id_t & aBlockID,
                id_t & aElementID,
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks ) ;

//------------------------------------------------------------------------------

            void
            create_elements_on_blocks_tri6(
                id_t & aBlockID,
                id_t & aElementID,
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks  ) ;

//------------------------------------------------------------------------------

            void
            link_elements_with_edges(
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks );

//------------------------------------------------------------------------------

            void
            link_elements_with_faces(
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks );

//------------------------------------------------------------------------------

            void
            create_buffers(
                const Cell< string > & aMaterials,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks );

//------------------------------------------------------------------------------

            // ghost facets for normal directions
            void
            create_ghost_facets(
                const uint aOrder,
                Cell< Facet * > & aFacets,
                Cell< Layer * > & aLayers,
                Cell< Block * > & aBlocks );

//------------------------------------------------------------------------------

            id_t
            max_node_id();

//------------------------------------------------------------------------------

            id_t
            max_element_id();

//------------------------------------------------------------------------------

            id_t
            max_block_id();

//------------------------------------------------------------------------------

            real
            compute_binomial_vectors(
                Protoshell * aProtoShell,
                Curve * aCurve,
                const Matrix< real > & aNormals,
                      Matrix< real > & aBinomials );

            void
            update_node_edge_tables( Cell< Node * > & aNodes, Cell< Edge * > & aEdges );

            void
            update_node_facet_tables( Cell< Node * > & aNodes, Cell< Facet * > & aFacets );

            void
            flag_layer_nodes( Cell< Node * > & aNodes, Cell< Layer * > & aLayers );

            SideSet *
            create_periodic_sideset(
                Cell< Block * > & aBlocks,
                const bool aMaster );
        };
    }
}
#endif //CL_FEM_THINSHELLFACTORY_HPP
