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

#include "cl_CurveFactory.hpp"

#include "fn_unique.hpp"
#include "cl_Element_Factory.hpp"
#include "commtools.hpp"
#include "fn_min.hpp"
#include "cl_Queue.hpp"

#include "op_Node_Index.hpp"
#include "op_Segment_Index.hpp"
#include "cl_EdgeFactory.hpp"
#include "intpoints.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        CurveFactory::CurveFactory( Mesh * aMesh ) :
            mCommRank( comm_rank() ),
            mMesh( aMesh ),
            mMaxElementID( aMesh->max_element_id() ),
            mMaxCurveID( aMesh->max_block_and_sideset_id() )
        {

        }

//------------------------------------------------------------------------------

        bool
        CurveFactory::intersection_exists( const id_t aSideSetA, const id_t aSideSetB )
        {
            BELFEM_ERROR( comm_rank() == 0, "At this time, curves are not intended to be used in parallel" );

            SideSet * tSideSetA = mMesh->sideset( aSideSetA );
            SideSet * tSideSetB = mMesh->sideset( aSideSetB );

            tSideSetA->unflag_all_nodes() ;
            tSideSetB->unflag_all_nodes() ;
            tSideSetA->flag_corner_nodes() ;

            index_t tCount = 0 ;
            for ( Node * tNode : tSideSetB->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    // we want to have at least two nodes
                    if ( ++tCount > 1 )
                    {
                        return true ;
                    }
                }
            }

            return false ;
        }

//------------------------------------------------------------------------------

        Cell< Curve * >
        CurveFactory::thin_shell_side_curves( const Vector< id_t > & aThinShellSideSets, const Vector< id_t > & aDomainBoundaries )
        {
            BELFEM_ERROR( comm_rank() == 0, "At this time, curves are not intended to be used in parallel" );

            // all nodes that sit on thin shells and are used to create the edges
            Cell< Node * >  tNodes ;
            this->collect_nodes( aThinShellSideSets, tNodes );

            Cell< Edge * > tAllEdges ;
            Map< key_t, Edge * > tMap ;

            this->create_edges_on_sidesets( aThinShellSideSets, tNodes, tAllEdges, tMap );

            this->connect_facets_to_edges( aThinShellSideSets, tNodes.size(), tAllEdges, tMap );

            Cell< Edge * > tEdges ;
            this->select_edges( aThinShellSideSets, aDomainBoundaries, tAllEdges, tEdges );
            this->select_node_subset( tEdges,  tNodes );
            this->connect_edges_to_nodes( tEdges, tNodes );
            this->identify_subchains( tNodes );

            Cell< ProtoCurve * > tProtoCurves ;
            this->create_protocurves( tNodes, tEdges, tProtoCurves );

            Cell< Curve * > aCurves( tProtoCurves.size(), nullptr );

            index_t tCount = 0 ;

            for ( ProtoCurve * tProtoCurve : tProtoCurves )
            {
                ElementType tType = interpolation_order_numeric(
                    tProtoCurve->mEdges( 0 )->facet( 0 )->element()->type() ) == 2 ?
                        ElementType::LINE3 : ElementType::LINE2 ;

                Cell< Node * > & tEndNodes = tProtoCurve->mNodes ;

                //Create node pairs on the side curve
                Cell< std::pair< Node *, Node * > > tPairs(tProtoCurve->mEdges.size(),std::make_pair( nullptr, nullptr )) ;
                for (uint i = 0 ; i < tProtoCurve->mEdges.size() ; ++i )
                {
                    Edge * tEdge = tProtoCurve->mEdges(i) ;
                    tPairs(i).first = tEdge->node(0) ;
                    tPairs(i).second = tEdge->node(1) ;
                }

                Cell< Node * >   tMidNodes ;
                if ( tType == ElementType::LINE3 )
                {
                    this->collect_midnodes( tPairs, tMidNodes );
                }

                // create the new curve
                Curve * tCurve = new Curve( ++mMaxCurveID, tType );

                // create the segments
                this->create_segments( tCurve, tPairs, tMidNodes ) ;

                Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
                this->sort_end_nodes( tCurve, tAdjacency, tEndNodes );

                // with the nodes sorted, we reorient the segments
                this->orient_segments( tEndNodes, tCurve->segments() ) ;

                // Step 8: we can now order the segments in the order they are connected
                this->sort_segments( tEndNodes( 0 ), tAdjacency, tCurve->segments() ) ;

                // Step 9: collect and sort the nodes
                this->collect_nodes( tCurve );

                this->compute_coordinates( tCurve );

                aCurves( tCount++ ) = tCurve ;
                delete tProtoCurve ;
            }

            for ( Node * tNode : tNodes )
            {
                tNode->reset_vertex_container() ;
            }
            for ( Edge * tEdge : tAllEdges )
            {
                delete tEdge ;
            }

            // restore node indices on mesh, since we have messed with them
            mMesh->update_node_indices() ;

            return aCurves ;
        }

        void
        CurveFactory::collect_nodes( const Vector< id_t > & aSideSets, Cell< Node * > & aNodes )
        {
            mMesh->unflag_all_nodes() ;
            for ( id_t tID : aSideSets )
            {
                mMesh->sideset( tID )->flag_corner_nodes() ;
            }

            index_t tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
                else
                {
                    tNode->set_index( gNoIndex );
                }
            }
            aNodes.set_size( tCount, nullptr ) ;
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    aNodes( tCount++ ) = tNode ;
                }
            }
        }

        void
        CurveFactory::create_edges_on_sidesets(
            const Vector< id_t > & aSideSets ,
                  Cell< Node * > & aNodes,
                  Cell< Edge * >       & aEdges,
                  Map< key_t, Edge * > & aMap  )
        {
            index_t tCount = 0 ;
            for ( id_t tID : aSideSets )
            {
                tCount += mMesh->sideset( tID )->number_of_facets() * number_of_edges( mMesh->sideset( tID )->element_type() );
            }

            Vector< key_t > tKeys( tCount ) ;

            tCount = 0 ;
            Cell< Node * > tNodes ;

            key_t tNumNodes = aNodes.size();

            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    for ( uint e=0; e<tFacet->element()->number_of_edges(); ++e )
                    {
                        tFacet->element()->get_nodes_of_edge( e, tNodes ) ;

                        // create the key
                        key_t tA = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 0 )->index() : tNodes( 1 )->index() ;
                        key_t tB = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 1 )->index() : tNodes( 0 )->index() ;
                        tKeys( tCount++ ) = tA * tNumNodes + tB ;
                    }
                }
            }

            unique( tKeys );

            tCount = 0 ;
            aEdges.set_size( tKeys.length(), nullptr ) ;
            for ( key_t tKey : tKeys )
            {
                // compute first and second node ids
                key_t tIndexA = tKey % tNumNodes ;
                key_t tIndexB = ( tKey - tIndexA ) / tNumNodes ;

                Node * tNodeA = aNodes( tIndexA );
                Node * tNodeB = aNodes( tIndexB );

                Edge * tEdge = new Edge ;
                tEdge->set_index( tCount );
                tEdge->allocate_node_container( 2 );
                tEdge->insert_node( tNodeA, 0 );
                tEdge->insert_node( tNodeB, 1 );
                aEdges( tCount++ ) = tEdge ;
                aMap[ tKey ] = tEdge ;
            }
        }

        void
        CurveFactory::connect_facets_to_edges(
            const Vector< id_t > & aSideSets,
            const key_t aNumNodes,
                Cell< Edge * > & aEdges,
                Map< key_t, Edge * > & aMap  )
        {
            Cell< Node * > tNodes ;

            // count facets per edge
            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;

                for ( Facet * tFacet : tFacets )
                {
                    for ( uint e=0; e<tFacet->element()->number_of_edges(); ++e )
                    {
                        // get the nodes of the edge
                        tFacet->element()->get_nodes_of_edge( e, tNodes ) ;

                        // create the key
                        key_t tA = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 0 )->index() : tNodes( 1 )->index() ;
                        key_t tB = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 1 )->index() : tNodes( 0 )->index() ;
                        key_t tKey = tA * aNumNodes + tB ;

                        // get the edge
                        aMap( tKey )->increment_facet_counter() ;
                    }
                }
            }

            for ( Edge * tEdge : aEdges )
            {
                tEdge->allocate_facet_container() ;
            }

            // connect facets to edges
            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;

                for ( Facet * tFacet : tFacets )
                {
                    for ( uint e=0; e<tFacet->element()->number_of_edges(); ++e )
                    {
                        // get the nodes of the edge
                        tFacet->element()->get_nodes_of_edge( e, tNodes ) ;

                        // create the key
                        key_t tA = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 0 )->index() : tNodes( 1 )->index() ;
                        key_t tB = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 1 )->index() : tNodes( 0 )->index() ;
                        key_t tKey = tA * aNumNodes + tB ;

                        // add facet to edge
                        aMap( tKey )->add_facet( tFacet ) ;
                    }
                }
            }
        }

        void
        CurveFactory::select_edges(
            const Vector< id_t > & aThinShellSideSets,
            const Vector< id_t > & aDomainBoundaries,
                      Cell< Edge * > & aAllEdges,
                      Cell< Edge * > & aSelectedEdges )
        {
            for ( id_t tID : aThinShellSideSets )
            {
                mMesh->sideset( tID )->flag_corner_nodes();
            }
            for ( id_t tID : aDomainBoundaries )
            {
                mMesh->sideset( tID )->unflag_all_nodes();
            }

            index_t tCount = 0 ;

            for ( Edge * tEdge : aAllEdges )
            {
                if ( tEdge->number_of_facets() == 1 )
                {
                    if ( tEdge->node( 0 )->is_flagged() or tEdge->node( 1 )->is_flagged() )
                    {
                        ++tCount ;
                    }
                }
            }

            aSelectedEdges.set_size( tCount, nullptr ) ;
            tCount = 0 ;

            for ( Edge * tEdge : aAllEdges )
            {
                if ( tEdge->number_of_facets() == 1 )
                {

                    if ( tEdge->node( 0 )->is_flagged() or tEdge->node( 1 )->is_flagged() )
                    {
                        aSelectedEdges( tCount++ ) = tEdge ;
                        tEdge->flag() ;
                    }
                    else
                    {
                        tEdge->unflag() ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::select_node_subset( Cell< Edge * > & aEdges, Cell< Node * > & aNodes )
        {
            // first, let's make sure that we only select the nodes we are still interested in
            for ( Node * tNode : aNodes )
            {
                tNode->unflag();
            }
            for ( Edge * tEdge : aEdges )
            {
                tEdge->node( 0 )->flag() ;
                tEdge->node( 1 )->flag() ;
            }

            // count flagged nodes
            index_t tCount = 0 ;

            for ( Node * tNode : aNodes )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
                else
                {
                    tNode->set_index( gNoIndex );
                }
            }

            // temporary node container
            Cell< Node * > tNodes;
            tNodes.vector_data() = std::move( aNodes.vector_data() );

            aNodes.set_size( tCount, 0 );
            tCount = 0 ;

            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    aNodes( tCount++ ) = tNode ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::connect_edges_to_nodes(
            Cell< Edge * > & aEdges,
            Cell< Node * > & aNodes )
        {
            // note: usually, we would use the edge containers on the nodes.
            // we can't do this here because edges might already exist on the mesh.
            // we use the vertex container instead which is not used by the mesh

            for ( Edge * tEdge : aEdges )
            {
                tEdge->node( 0 )->increment_vertex_counter() ;
                tEdge->node( 1 )->increment_vertex_counter() ;
            }

            for ( Node * tNode : aNodes )
            {
                tNode->init_vertex_container();
            }

            for ( Edge * tEdge : aEdges )
            {
                tEdge->node( 0 )->insert_vertex( tEdge );
                tEdge->node( 1 )->insert_vertex( tEdge );
            }
        }

        index_t
        CurveFactory::identify_subchains( Cell< Node * > & aNodes )
        {
            Queue< Node * > tQueue ;
            Queue< Node * > tPending ;

            for ( Node * tNode : aNodes )
            {
                tNode->set_index( gNoIndex );
                tPending.push( tNode );
            }


            index_t aGroup = 0 ;
            index_t tIteration = 0 ;

            while ( ! tPending.empty() )
            {
                Node * tStart = tPending.pop();

                if ( tStart->index() != gNoIndex ) continue;

                tQueue.push( tStart ) ;
                tStart->set_index( ++aGroup ) ;

                while ( ! tQueue.empty() )
                {
                    Node * tNode = tQueue.pop() ;

                    for ( uint e=0; e<tNode->number_of_vertices(); ++e )
                    {
                        Edge * tEdge = reinterpret_cast< Edge * >( tNode->vertex( e ) ) ;

                        // put nodes on queue if they have not been visited
                        if ( tEdge->node( 0 )->index() == gNoIndex )
                        {
                            tEdge->set_index( aGroup );
                            tEdge->node( 0 )->set_index( aGroup );
                            tQueue.push( tEdge->node( 0 ) );
                        }
                        if ( tEdge->node( 1 )->index() == gNoIndex )
                        {
                            tEdge->set_index( aGroup );
                            tEdge->node( 1 )->set_index( aGroup );
                            tQueue.push( tEdge->node( 1 ) );
                        }
                    }
                }
                BELFEM_ERROR( tIteration++ < aNodes.size(), "Failed to identify subgraphs" );
            }

            return aGroup ;
        }


        void
        CurveFactory::create_protocurves(
            Cell< Node * >     & aNodes,
            Cell< Edge * >     & aEdges,
            Cell< ProtoCurve *  > & aProtoCurves )
        {
            index_t tNumChains = this->identify_subchains( aNodes );

            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tNode->index() - 1 );
                tNode->unflag();
            }
            for ( Edge * tEdge : aEdges )
            {
                tEdge->set_index( tEdge->index() - 1 );
                tEdge->unflag();
            }

            // count nodes per chain
            Vector< index_t > tNodeCount( tNumChains, 0 );

            for ( Node * tNode : aNodes )
            {
                ++tNodeCount( tNode->index() );
            }

            // count edges per chain
            Vector< index_t > tEdgeCount( tNumChains, 0 );
            for ( Edge * tEdge : aEdges )
            {
                ++tEdgeCount( tEdge->index() );
            }

            aProtoCurves.set_size( tNumChains, nullptr );
            for ( index_t k=0; k<tNumChains; ++k )
            {
                ProtoCurve * tProtoCurve = new ProtoCurve;
                tProtoCurve->mNodes.set_size( tNodeCount( k ), nullptr );
                tProtoCurve->mEdges.set_size( tEdgeCount( k ), nullptr );
                tNodeCount( k ) = 0 ;
                tEdgeCount( k ) = 0 ;
                aProtoCurves( k ) = tProtoCurve;
            }

            for ( Node * tNode : aNodes )
            {
                aProtoCurves( tNode->index() )->mNodes( tNodeCount( tNode->index() ) ) = tNode ;
                tNode->set_index( tNodeCount( tNode->index() )++ );

            }

            for ( Edge * tEdge : aEdges )
            {
                aProtoCurves( tEdge->index() )->mEdges( tEdgeCount( tEdge->index() ) ) = tEdge ;
                tEdge->set_index( tEdgeCount( tEdge->index() )++ );
            }
        }


//------------------------------------------------------------------------------

        Curve *
        CurveFactory::intersect( const id_t aSideSetA, const id_t aSideSetB, id_t aID )
        {
            BELFEM_ERROR( comm_rank() == 0, "Curves are not intended to be used in parallel" );

            // Step 1 : create the curve object
            Curve * aCurve = new Curve(
                aID,
                mMesh->sideset( aSideSetA ),
                mMesh->sideset( aSideSetB ) );

            ++mMaxCurveID ;

            // Step 2: collects shared corner nodes between side sets
            Cell< Node * > tEndNodes ;
            this->collect_end_nodes_from_intersection( aCurve->sideset_a(), aCurve->sideset_b(), tEndNodes ) ;

            // Step 3: forms edges from corner nodes
            Cell< std::pair< Node *, Node * > > tPairs ;
            this->create_pairs( tEndNodes, tPairs ) ;

            // only nodes on sideset A have been flagged so far
            // let's reset them to maintain a clean dataset
            aCurve->sideset_a()->unflag_all_nodes() ;

            // Step 4 : Gathers mid-side nodes for quadratic elements
            Cell< Node * > tMidNodes ;
            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                this->collect_midnodes( tPairs, tMidNodes ) ;
            }

            // Step 5: create the segments
            this->create_segments( aCurve, tPairs, tMidNodes ) ;

            // Step 6: put the end nodes into the appropriate order
            Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
            this->sort_end_nodes( aCurve, tAdjacency, tEndNodes );

            // Step 7: with the nodes sorted, we reorient the segments
            this->orient_segments( tEndNodes, aCurve->segments() ) ;

            // Step 8: we can now order the segments in the order they are connected
            this->sort_segments( tEndNodes( 0 ), tAdjacency, aCurve->segments() ) ;

            // Step 9: collect and sort the nodes
            this->collect_nodes( aCurve );

            // Step 10 : compute the lengths of the curve
            this->compute_coordinates( aCurve );

            // Step 11 : restore node indices on mesh, since we have messed with them
            mMesh->update_node_indices();

            return aCurve ;
        }

        Curve *
        CurveFactory::intersect( const id_t aThinShellSideSet, const id_t aBoundarySideSet, DynamicBitset & aBoundaryNodeBitset )
        {
            // Step 1 : create the curve object
            Curve * aCurve = new Curve(
                ++mMaxCurveID,
                mMesh->sideset( aBoundarySideSet ),
                mMesh->sideset( aThinShellSideSet ) );

            // Step 1: collects shared corner nodes between side sets
            DynamicBitset tEndNodeBitset( mMesh->nodes().size() );
            DynamicBitset tMidNodeBitset( mMesh->nodes().size() );

            Cell< Facet * > & tFacet = mMesh->sideset( aThinShellSideSet )->facets() ;
            for ( Facet * tFacet : tFacet )
            {
                for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                {
                    tEndNodeBitset.set( tFacet->node( k )->index() ) ;
                }
            }

            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                for ( Facet * tFacet : tFacet )
                {
                    for ( uint k=tFacet->number_of_corner_nodes(); k<tFacet->number_of_nodes(); ++k )
                    {
                        tMidNodeBitset.set( tFacet->node( k )->index() ) ;
                    }
                }
            }

            // get intersection with boundary nodes
            tEndNodeBitset &= aBoundaryNodeBitset ;

            index_t tNumEndNodes = tEndNodeBitset.count() ;

            if ( tNumEndNodes == 0 )
            {
                return nullptr ;
            }


            index_t tNumMidNodes = 0 ;
            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                tMidNodeBitset &= aBoundaryNodeBitset ;
                tNumMidNodes = tMidNodeBitset.count() ;
            }


            // reset all node indices
            Cell< Node * > & tNodes = mMesh->nodes() ;
            for ( Node * tNode : tNodes )
            {
                tNode->set_index( gNoIndex );
            }

            // set node indices for end nodes
            Cell< index_t > tIndices ;
            tEndNodeBitset.where( tIndices );

            index_t tCount = 0 ;
            Cell< Node * > tEndNodes( tNumEndNodes, nullptr );
            for ( index_t tIndex : tIndices )
            {
                tEndNodes( tCount ) = tNodes( tIndex ) ;
                tNodes( tIndex )->set_index( tCount++ );
            }

            // set node indices for mid nodes
            Cell< Node * > tMidNodes ;
            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                tIndices.clear() ;
                tMidNodeBitset.where( tIndices );
                tCount = 0 ;
                tMidNodes.set_size( tNumMidNodes, nullptr );
                for ( index_t tIndex : tIndices )
                {
                    tMidNodes( tCount++ ) = tNodes( tIndex ) ;
                }
            }

            // Step 3: forms edges from corner nodes
            Cell< std::pair< Node *, Node * > > tPairs ;
            this->create_pairs( tEndNodes, tPairs ) ;

            // Step 5: create the segments
            this->create_segments( aCurve, tPairs, tMidNodes ) ;

            // Step 6: put the end nodes into the appropriate order
            Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
            this->sort_end_nodes( aCurve, tAdjacency, tEndNodes );

            // Step 7: with the nodes sorted, we reorient the segments
            this->orient_segments( tEndNodes, aCurve->segments() ) ;

            // Step 8: we can now order the segments in the order they are connected
            this->sort_segments( tEndNodes( 0 ), tAdjacency, aCurve->segments() ) ;

            // Step 9: collect and sort the nodes
            this->collect_nodes( aCurve );

            // Step 10 : compute the lengths of the curve
            this->compute_coordinates( aCurve );

            // Step 11 : restore node indices on mesh, since we have messed with them
            mMesh->update_node_indices();

            return aCurve ;
        }

        Curve *
        CurveFactory::from_2d_sidesets( const Vector< id_t > & aSideSets, id_t aID )
        {
            BELFEM_ERROR( comm_rank() == 0, "Curves are not intended to be used in parallel" );

            // Step 1: Collect the facets
            Cell< Facet * > tFacets ;
            index_t tCount = 0 ;
            for ( id_t tID : aSideSets )
            {
                tCount+=mMesh->sideset( tID )->number_of_facets() ;
            }
            tFacets.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tSideSet = mMesh->sideset( tID )->facets() ;
                for ( Facet * tFacet : tSideSet )
                {
                    tFacets( tCount++ ) = tFacet ;
                }
            }

            // Step 2: create the curve object
            Curve * aCurve = new Curve( aID, mMesh->sideset( aSideSets(0) ), nullptr );
            ++mMaxCurveID ;

            // Step 3: collect the corner nodes
            Cell< Node * > tEndNodes( 2*tFacets.size(), nullptr );
            tCount = 0 ;
            //tEndNodes( tCount++ ) = tFacets( 0 )->node( 0 );
            for ( Facet * tFacet : tFacets )
            {
                tEndNodes( tCount++ ) = tFacet->node( 0 );
                tEndNodes( tCount++ ) = tFacet->node( 1 );
            }
            //remove the doubles
            unique( tEndNodes );
            for ( Node * tNode : tEndNodes )
            {
                for ( uint k=0; k<tNode->number_of_nodes(); ++k )
                {
                    tNode->node( k )->set_index( gNoIndex );
                }
            }
            tCount = 0 ;
            for ( Node * tNode : tEndNodes )
            {
                tNode->set_index( tCount++ );
                tNode->unflag() ;
            }

            // Step 4 : Gathers mid-side nodes for quadratic elements
            Cell< Node * > tMidNodes ;
            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                tCount = 0 ;
                tMidNodes.set_size( tFacets.size(), nullptr );
                for ( Facet * tFacet : tFacets )
                {
                    tMidNodes( tCount++ ) = tFacet->node( 2 );
                }
            }

            // Step 4: Create the Pairs
            Cell< std::pair< Node *, Node * > > tPairs ;
            this->create_pairs( tEndNodes, tPairs ) ;

            // Step 5: create the segments
            this->create_segments( aCurve, tPairs, tMidNodes ) ;

            // Step 6: put the end nodes into the appropriate order
            Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
            this->sort_end_nodes( aCurve, tAdjacency, tEndNodes );

            // Step 7: with the nodes sorted, we reorient the segments
            this->orient_segments( tEndNodes, aCurve->segments() ) ;

            // Step 8: we can now order the segments in the order they are connected
            this->sort_segments( tEndNodes( 0 ), tAdjacency, aCurve->segments() ) ;

            // Step 9: collect and sort the nodes
            this->collect_nodes( aCurve );

            // Step 10 : compute the lengths of the curve
            this->compute_coordinates( aCurve );

            // Step 11 : restore node indices on mesh, since we have messed with them
            mMesh->update_node_indices();

            return aCurve ;
        }
//------------------------------------------------------------------------------

        void
        CurveFactory::collect_end_nodes_from_intersection(
                SideSet * aSidesetA,
                SideSet * aSideSetB,
                Cell< Node * > & aEndNodes )
        {
            // voiding the node indices will provoke a crash
            // if we did something wrong.
            // the value of gNoIndex is 2^32-1
            for ( Node * tNode : mMesh->nodes() )
            {
                tNode->unflag();
                tNode->set_index( gNoIndex );
            }

            aSidesetA->flag_corner_nodes() ;

            index_t tCount = 0 ;

            for ( Node * tNode : aSideSetB->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }

            aEndNodes.set_size( tCount, nullptr ) ;
            tCount = 0 ;
            for ( Node * tNode : aSideSetB->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    aEndNodes( tCount++ ) = tNode ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::create_pairs( Cell< Node * > & aEndNodes, Cell< std::pair< Node *, Node * > > & aPairs )
        {
            // count possible edges
            index_t tCount = 0 ;
            for ( Node * tNode : aEndNodes )
            {
                for ( uint k=0; k<tNode->number_of_nodes(); ++k )
                {
                    if ( tNode->node( k )->index() != gNoIndex )
                    {
                        ++tCount ;
                    }
                }
            }
            key_t tNumNodes = aEndNodes.size();

            Vector< key_t > tKeys( tCount ) ;
            tCount = 0 ;
            for ( Node * tNodeA : aEndNodes )
            {
                for ( uint k=0; k<tNodeA->number_of_nodes(); ++k )
                {
                    Node * tNodeB = tNodeA->node( k ) ;
                    if ( tNodeB->index() != gNoIndex )
                    {
                        key_t tA = tNodeA->index() > tNodeB->index() ? tNodeA->index() : tNodeB->index() ;
                        key_t tB = tNodeA->index() > tNodeB->index() ? tNodeB->index() : tNodeA->index() ;
                        tKeys( tCount++ ) = tA * tNumNodes + tB ;
                    }
                }
            }

            unique( tKeys );


            aPairs.reserve( tKeys.length() );

            for ( key_t tKey : tKeys )
            {
                // recover first and second node indices
                key_t tIndexA = tKey % tNumNodes ;
                key_t tIndexB = ( tKey - tIndexA ) / tNumNodes ;

                Node * tNodeA = aEndNodes( tIndexA );
                Node * tNodeB = aEndNodes( tIndexB );

                aPairs.push( std::make_pair( tNodeA, tNodeB ) );
            }

        }

//------------------------------------------------------------------------------

        void
        CurveFactory::collect_midnodes( Cell< std::pair< Node *, Node * > > & aPairs, Cell< Node * > & aMidNodes )
        {
            aMidNodes.set_size( aPairs.size(), nullptr ) ;

            index_t tCount = 0 ;

            id_t tA ;
            id_t tB ;
            bool tFound ;

            Cell< Node * > tNodes ;

            for ( auto tPair : aPairs )
            {
                tFound  = false ;
                tA = tPair.first->id() ;
                tB = tPair.second->id() ;

                for ( uint e=0; e<tPair.first->number_of_elements(); ++e )
                {
                    Element * tElement = tPair.first->element( e ) ;

                    for ( uint d=0; d<tElement->number_of_edges(); ++d )
                    {
                        tElement->get_nodes_of_edge( d, tNodes ) ;

                        if (   ( tNodes( 0 )->id() == tA && tNodes( 1 )->id() == tB )
                            || ( tNodes( 0 )->id() == tB && tNodes( 1 )->id() == tA ) )
                        {
                            // setting the index to none is just a precaution, we want to provoke
                            // a crash in the sorting algorithm in which the midnodes are not relevant
                            tNodes( 2 )->set_index( gNoIndex );

                            // add node to container
                            aMidNodes( tCount++ ) = tNodes( 2 ) ;
                            tFound = true ;
                            break ;
                        }
                    }

                    if ( tFound ) break ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::create_segments( Curve * aCurve, Cell< std::pair< Node *, Node * > > & aPairs, Cell< Node * > & aMidNodes )
        {
            Cell< Segment * > & tSegments = aCurve->segments() ;

            tSegments.set_size( aPairs.size(), nullptr ) ;

            ElementFactory tFactory ;

            index_t tCount = 0 ;

            for (  auto tPair : aPairs )
            {
                Element * tElement = tFactory.create_element( aCurve->element_type(), ++mMaxElementID );
                tElement->insert_node( tPair.first, 0 );
                tElement->insert_node( tPair.second, 1 );

                tElement->set_geometry_tag( aCurve->id() );
                Segment * tSegment = new Segment( tElement );
                tSegment->set_index( tCount );

                tSegments( tCount++ ) = tSegment ;
            }
            // add midside node of we are second order
            tCount = 0 ;
            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                for ( Segment * tSegment : tSegments )
                {
                    tSegment->insert_node( aMidNodes( tCount++ ), 2 );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::sort_end_nodes( Curve * aCurve, Matrix< index_t > & aAdjacency, Cell< Node * > & aEndNodes )
        {
            Cell< Segment * > & tSegments = aCurve->segments() ;

            // counter for nodes
            Vector< uint > tNumSegments( aEndNodes.size(), 0 ) ;

            // the first loop is to determine if the loop is closed and sane
            for ( Segment * tSegment : tSegments )
            {
                tNumSegments( tSegment->node( 0 )->index() )++ ;
                tNumSegments( tSegment->node( 1 )->index() )++ ;
            }


            Node * tStart = nullptr ;
            Node * tEnd = nullptr ;

            index_t tOne = 0 ;
            index_t tTwo = 0 ;
            index_t tCount = 0 ;

            for ( index_t tN : tNumSegments )
            {
                if ( tN == 1 )
                {
                    if ( tStart == nullptr )
                    {
                        tStart = aEndNodes( tCount );
                    }
                    else
                    {
                        if ( tStart->id() > aEndNodes( tCount )->id() )
                        {
                            tEnd = tStart ;
                            tStart = aEndNodes( tCount );
                        }
                        else
                        {
                            tEnd = aEndNodes( tCount );
                        }
                    }
                    ++tOne ;
                }
                else if ( tN == 2 )
                {
                    ++tTwo ;
                }
                ++tCount ;
            }
            BELFEM_ERROR( tOne + tTwo == tCount, "The Curve seems to branch or is not properly connected. This is not allowed: 1: %lu 2: %lu 1+2: %lu",
                ( luint ) tOne, ( luint ) tTwo, ( luint ) tCount );

            // reset the counter
            tNumSegments.fill( 0 );

            // create the adjacency, filled with gNoIndex: next() tests the
            // slots against gNoIndex. The fill must go through set_size, whose
            // fill value is typed T -- the constructor fill is typed real, and
            // a 64-bit gNoIndex does not fit double's 53-bit mantissa
            Matrix< index_t > tAdjacency ;
            tAdjacency.set_size( aEndNodes.size(), 2, gNoIndex );

            for ( Segment * tSegment : tSegments )
            {
                tAdjacency( tSegment->node( 0 )->index(), tNumSegments( tSegment->node( 0 )->index() )++ ) = tSegment->index() ;
                tAdjacency( tSegment->node( 1 )->index(), tNumSegments( tSegment->node( 1 )->index() )++ ) = tSegment->index() ;
            }

            // if the loop is closed, we still have no starting point, let's find the one with the smallest id
            aCurve->set_closed_flag( tOne == 0 );

            if ( tStart == nullptr )
            {
                // find the point with the smallest id as starting point
                tStart = aEndNodes( 0 ) ;
                for ( Node * tNode : aEndNodes )
                {
                    if ( tStart->id() > tNode->id() )
                    {
                        tStart = tNode ;
                    }
                }
                tEnd = tStart ;
            }

            // we know that all nodes are unflagged
            // let's flag the end
            tEnd->flag();

            index_t tNumCornerNodes = aEndNodes.size();

            // the counter for the segments

            tCount = 0 ;
            Node * tNode = tStart ;

            Vector< index_t > tNewIndices( tNumCornerNodes, gNoIndex ) ;

            // set the new index of the first node
            tNewIndices( tNode->index() ) = tCount++ ;

            for ( index_t k=0; k<tNumCornerNodes; ++k )
            {
                // flag the node
                tNode->flag();

                // increment the node index
                tNewIndices( tNode->index() ) = k ;

                // get the next node
                tNode = this->next( tAdjacency, tSegments, tNode ) ;


                if ( tNode->is_flagged() )
                {
                    break ;
                }
            }

            // note that if we are in a loop, the index of the start node must be
            // set to zero again
            tNewIndices( tEnd->index() ) = tNumCornerNodes-1 ;
            tNewIndices( tStart->index() ) = 0 ;

            // next, let's write the indices into the nodes
            tCount = 0 ;
            for ( Node * tNode : aEndNodes )
            {
                tNode->set_index( tNewIndices( tCount++ ) );
            }

            // now, let's sort the nodes
            sort( aEndNodes, opNodeIndex );

            // finally, we must reorganize the adjacency
            aAdjacency.set_size( tNumCornerNodes, 2, gNoIndex );
            tCount = 0 ;
            for ( index_t k=0; k<tNumCornerNodes; ++k )
            {
                aAdjacency( tNewIndices( k ), 0 ) = tAdjacency( k, 0 ) ;
                aAdjacency( tNewIndices( k ), 1 ) = tAdjacency( k, 1 ) ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::orient_segments(  Cell< Node * > & aEndNodes, Cell< Segment * > & aSegments )
        {
            // with the nodes sorted, we can now adjust the segments
            // first, we orient them
            for ( Segment * tSegment : aSegments )
            {
                if ( tSegment->node( 0 )->index() > tSegment->node( 1 )->index() )
                {
                    Node * tSwap = tSegment->node( 0 ) ;

                    // the segment doesn't have an insert function
                    // we access the wrapped line element
                    tSegment->element()->insert_node(  tSegment->node( 1 ), 0 ) ;
                    tSegment->element()->insert_node(  tSwap, 1 ) ;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        CurveFactory::sort_segments( Node * aStart, const Matrix< index_t > & aAdjacency, Cell< Segment * > & aSegments )
        {
            // get the first segment
            Segment * tSegment = aSegments( aAdjacency( aStart->index(), 0 ) );
            tSegment->set_index( 0 );
            tSegment->flag();

            index_t tNumSegments = aSegments.size();
            index_t tCount = 1 ;

            while ( tCount < tNumSegments )
            {
                // get the next segment
                Segment * tNext = this->next( aAdjacency, aSegments, tSegment ) ;

                // if the segment has been visited, we know that we have closed the loop
                if ( tNext->is_flagged() ) break ;

                // write index of next segment
                tNext->set_index( tCount++ );

                // mark segment as visited
                tNext->flag();

                // shift segment
                tSegment = tNext ;
            }

            sort( aSegments, opSegmentIndex );
        }


//------------------------------------------------------------------------------

         /**
         * Determines the next node in the curve traversal based on the current node and its adjacency information.
         *
         * This function is part of a curve traversal algorithm. It handles two cases:
         * - For nodes with one segment (endpoints of an open curve), it returns the only neighboring node.
         * - For nodes with two segments (internal nodes or nodes in a closed curve), it selects the unflagged neighbor
         *   to continue the traversal. If both neighbors are unflagged, it picks the one with the larger ID for consistency.
         *
         * @param aAdjacency Adjacency matrix where each row corresponds to a node and lists indices of connected segments.
         * @param aSegments List of all segments in the curve.
         * @param aNode Pointer to the current node in the traversal.
         * @return Pointer to the next node in the traversal.
         */
        Node *
        CurveFactory::next( const Matrix< index_t > & aAdjacency, Cell< Segment * > & aSegments, const Node * aNode )
        {
            // Get the indices of the segments connected to the current node from the adjacency matrix
            index_t tS1 = aAdjacency( aNode->index(), 0 );  // Index of the first connected segment
            index_t tS2 = aAdjacency( aNode->index(), 1 );  // Index of the second connected segment (gNoIndex if none)

            if ( tS2 == gNoIndex )  // Case 1: Node has only one segment (an endpoint of an open curve)
            {
                Segment * tSegment = aSegments( tS1 );  // Retrieve the single connected segment
                // Return the neighboring node of the segment, which is the only possible next step
                if ( tSegment->node( 0 )->id() == aNode->id() )
                {
                    return tSegment->node( 1 );  // Current node is node 0, so return node 1
                }
                else
                {
                    return tSegment->node( 0 );  // Current node is node 1, so return node 0
                }
            }
            else  // Case 2: Node has two segments (internal node or part of a closed curve)
            {
                Segment * tSegment1 = aSegments( tS1 );  // First connected segment
                Segment * tSegment2 = aSegments( tS2 );  // Second connected segment

                // Identify the neighboring nodes from each segment
                Node * aX = tSegment1->node( 0 )->id() == aNode->id() ? tSegment1->node( 1 ) : tSegment1->node( 0 );
                Node * aY = tSegment2->node( 0 )->id() == aNode->id() ? tSegment2->node( 1 ) : tSegment2->node( 0 );

                // Decide the next node based on flagging and ID:
                // - Prefer the unflagged neighbor to move forward in the traversal.
                // - If both are unflagged (e.g., at the start), choose the one with the larger ID for consistency.
                if ( aX->id() > aY->id() )
                {
                    // aX has a larger ID:
                    // - If aX is flagged (visited), return aY (the unvisited neighbor).
                    // - If aX is unflagged, return aX (arbitrary choice when both unflagged, overridden by traversal logic later).
                    return aX->is_flagged() ? aY : aX;
                }
                else
                {
                    // aY has a larger or equal ID:
                    // - If aY is flagged (visited), return aX (the unvisited neighbor).
                    // - If aY is unflagged, return aY.
                    return aY->is_flagged() ? aX : aY;
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * Finds the next unflagged segment in the curve traversal.
         *
         * This function examines segments connected to the current segment's endpoints (nodes) and selects
         * an unflagged (unvisited) segment to continue the traversal. It supports both open and closed curves,
         * returning the current segment if no valid next segment exists to signal the end of traversal.
         *
         * @param aAdjacency Adjacency matrix mapping each node to indices of connected segments.
         * @param aSegments List of all segments in the curve.
         * @param aSegment The current segment being processed.
         * @return Pointer to the next unflagged segment, or the current segment if traversal should end.
         */
        Segment *
        CurveFactory::next( const Matrix< index_t > & aAdjacency, Cell< Segment * > & aSegments, Segment * aSegment )
        {
            // get the indices of the neigboring segments
            index_t tIndices[ 4 ];

            // First segment at start node (node 0)
            tIndices[ 0 ] = aAdjacency( aSegment->node( 0 )->index(), 0 );

            // Second segment at start node (node 0)
            tIndices[ 1 ] = aAdjacency( aSegment->node( 0 )->index(), 1 );

            // First segment at end node (node 1)
            tIndices[ 2 ] = aAdjacency( aSegment->node( 1 )->index(), 0 );

            // Second segment at end node (node 1)
            tIndices[ 3 ] = aAdjacency( aSegment->node( 1 )->index(), 1 );

            // Search for a candidate segment connected to the start node (node 0)
            Segment * aCandidateA = nullptr ;
            for ( uint s=0; s<2; ++s )
            {
                // an empty slot means the node carries only one segment: nothing
                // to consider. The current segment cannot stand in for it:
                // sort_segments has already overwritten aSegment->index() with its
                // traversal counter, so aSegments( aSegment->index() ) can be an
                // unrelated segment
                if ( tIndices[ s ] == gNoIndex ) continue ;

                // index-as-position is still valid here: the cell is not reordered
                // until the final sort in sort_segments
                Segment * tSegment = aSegments( tIndices[ s ] );

                // Select this segment if it’s unflagged (unvisited) and not the current segment
                if ( ! tSegment->is_flagged() && tSegment->id() != aSegment->id() )
                {
                    aCandidateA = tSegment ;
                }
            }

            // Search for a candidate segment connected to the end node (node 1)
            Segment * aCandidateB = nullptr ;
            for ( uint s=2; s<4; ++s )
            {
                // empty slot: skip, see above
                if ( tIndices[ s ] == gNoIndex ) continue ;

                Segment * tSegment = aSegments( tIndices[ s ] );

                // Select this segment if it’s unflagged (unvisited) and not the current segment
                if ( ! tSegment->is_flagged() && tSegment->id() != aSegment->id() )
                {
                    aCandidateB = tSegment ;
                }
            }

            // Determine which segment to return based on available candidates
            if ( aCandidateA == nullptr )
            {
                if ( aCandidateB == nullptr )
                {
                    // No unflagged neighbors found; return the current segment
                    // itself to end the traversal loop. It is flagged, so the
                    // caller breaks before writing another index
                    return aSegment ;
                }
                else
                {
                    // Only end node has a valid candidate; return it
                    return aCandidateB ;
                }
            }
            else if ( aCandidateB == nullptr )
            {
                // Only start node has a valid candidate; return it
                return aCandidateA ;
            }
            else
            {
                // Both nodes have candidates; return the one with the smaller ID for consistency
                return aCandidateA->id() < aCandidateB->id() ? aCandidateA : aCandidateB ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::collect_nodes( Curve * aCurve )
        {
            Cell< Segment * > & tSegments = aCurve->segments() ;
            Cell< Node * >    & tNodes    = aCurve->nodes() ;

            // let's collect the nodes, but first put them in the correct order
            index_t tCount = 0 ;
            tSegments( 0 )->node( 0 )->set_index( tCount++ ) ;
            switch (  aCurve->element_type() )
            {
                case( ElementType::LINE2 ):
                {
                    for ( Segment * tSegment : tSegments )
                    {
                        tSegment->node( 1 )->set_index( tCount++ ) ;
                    }
                    break ;
                }
                case( ElementType::LINE3 ):
                {
                    for ( Segment * tSegment : tSegments )
                    {
                        // note that in a second order line,
                        //  the nodes areenumerated as follows:
                        // (0) --- (2) --- (1)
                        // therefore, the 0-2-1 order below is NOT a bug!
                        tSegment->node( 2 )->set_index( tCount++ ) ;
                        tSegment->node( 1 )->set_index( tCount++ ) ;
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "The element type is not supported" );
                }
            }

            // finally, we need to fix the last node if we have a closed loop
            if (  aCurve->is_closed() )
            {
                BELFEM_ASSERT( tSegments.first()->node( 0 )->id() == tSegments.last()->node( 1 )->id(),
                    "The loop is supposed to be closed but something went wrong");

                // fix the first node id
                tSegments( 0 )->node( 0 )->set_index( 0 ) ;
            }

            // now we can put the nodes into the container
            tNodes.set_size( tCount, nullptr ) ;
            tCount = 0 ;

            // add the first node
            tNodes( tCount++ ) = tSegments( 0 )->node( 0 ) ;

            switch (  aCurve->element_type() )
            {
                case( ElementType::LINE2 ):
                {
                    for ( Segment * tSegment : tSegments )
                    {
                        tNodes( tCount++ ) = tSegment->node( 1 ) ;
                    }
                    break ;
                }
                case( ElementType::LINE3 ):
                {
                    for ( Segment * tSegment : tSegments )
                    {
                        // note that in a second order line,
                        //  the nodes areenumerated as follows:
                        // (0) --- (2) --- (1)
                        // therefore, the 0-2-1 order below is NOT a bug!
                        tNodes( tCount++ ) = tSegment->node( 2 ) ;
                        tNodes( tCount++ ) = tSegment->node( 1 ) ;
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "You should not never read this!" );
                }
            }

            BELFEM_ASSERT( tCount == tNodes.size(), "The number of nodes is not correct (is %u expect %u)", ( uint ) tCount, ( uint ) tNodes.size());
        }

        void
        CurveFactory::compute_coordinates( Curve * aCurve )
        {
            real s = 0 ;
            uint tCount = 0 ;

            Vector< real >    & tS = aCurve->arclength();
            Cell< Segment * > & tSegments = aCurve->segments() ;

            tS.set_size( aCurve->nodes().size() ) ;

            tS( tCount++ ) = s ;

            Vector< real > d( 3 );

            if ( aCurve->element_type() == ElementType::LINE2 )
            {
                for ( Segment * tSegment : tSegments )
                {
                    for ( uint i=0; i<3; ++i )
                    {
                        d( i ) = tSegment->node( 1 )->x( i ) - tSegment->node( 0 )->x( i ) ;
                    }
                    s += norm( d );
                    tS( tCount++ ) = s ;
                }
            }
            else if ( aCurve->element_type() == ElementType::LINE3 )
            {
                // number of integgration points
                int n = 7 ;
                Vector< double > w( n );
                Vector< double > xi( n );

                intpoints_gauss( &n, w.data(), xi.data() ) ;

                Vector< real > a( 3 );
                Vector< real > b( 3 );

                for ( Segment * tSegment : tSegments )
                {
                    for ( uint i=0; i<3; ++i )
                    {
                        a( i ) =       tSegment->node( 0 )->x( i )
                                 +     tSegment->node( 1 )->x( i )
                                 - 2 * tSegment->node( 2 )->x( i ) ;

                        b( i ) = 0.5 * ( tSegment->node( 1 )->x( i )
                                       - tSegment->node( 0 )->x( i ) ) ;
                    }

                    real l = 0 ;
                    for ( int k=0; k<n; ++k )
                    {
                        d( 0 ) = a( 0 ) * xi( k ) + b( 0 );
                        d( 1 ) = a( 1 ) * xi( k ) + b( 1 ) ;
                        d( 2 ) = a( 2 ) * xi( k ) + b( 2 ) ;

                        l += w( k ) * norm( d );
                    }

                    tS( tCount++ ) = s + 0.5 * l ;
                    s += l ;
                    tS( tCount++ ) = s ;
                }
            }
        }
//------------------------------------------------------------------------------
    }
}