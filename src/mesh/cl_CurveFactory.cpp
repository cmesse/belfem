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

#include <memory>
#include <string>

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
            const index_t tNumChains = this->identify_subchains( tNodes );

            // undoes the scratch state of this function on the mesh nodes
            auto tReleaseScratch = [ & ]()
            {
                for ( Node * tNode : tNodes )
                {
                    tNode->reset_vertex_container() ;
                }
                for ( Edge * tEdge : tAllEdges )
                {
                    delete tEdge ;
                }
                mMesh->update_node_indices() ;
            };

            // a component with at least as many edges as nodes contains a
            // loop. create_protocurves() cannot place the edge that closes
            // it, so a closed side curve is refused here. identify_subchains()
            // numbers the components from 1
            {
                Cell< index_t > tNodeCount( tNumChains, 0 );
                Cell< index_t > tDegreeSum( tNumChains, 0 );

                for ( Node * tNode : tNodes )
                {
                    ++tNodeCount( tNode->index() - 1 );
                    tDegreeSum( tNode->index() - 1 ) += tNode->number_of_vertices() ;
                }

                for ( Node * tNode : tNodes )
                {
                    const index_t tChain = tNode->index() - 1 ;
                    if ( tDegreeSum( tChain ) >= 2 * tNodeCount( tChain ) )
                    {
                        const id_t tNodeID = tNode->id() ;
                        tReleaseScratch() ;
                        BELFEM_ERROR( false,
                            "The thin-shell side curve through node %lu is closed; closed side curves are not supported",
                            ( long unsigned int ) tNodeID );
                    }
                }
            }

            Cell< ProtoCurve * > tProtoCurves ;
            this->create_protocurves( tNodes, tEdges, tProtoCurves );

            // the midpoints of every quadratic side curve, read from the one
            // facet that carries each edge; gathered before any curve exists,
            // so that a refusal has only the scratch state to undo
            Cell< ElementType > tTypes( tProtoCurves.size(), ElementType::EMPTY );
            Cell< Cell< Node * > > tMidNodes( tProtoCurves.size(), Cell< Node * >() );
            {
                Cell< Node * > tWork ;

                for ( index_t k=0; k<tProtoCurves.size(); ++k )
                {
                    Cell< Edge * > & tChainEdges = tProtoCurves( k )->mEdges ;

                    tTypes( k ) = interpolation_order_numeric(
                        tChainEdges( 0 )->facet( 0 )->element()->type() ) == 2 ?
                            ElementType::LINE3 : ElementType::LINE2 ;

                    if ( tTypes( k ) != ElementType::LINE3 ) continue ;

                    tMidNodes( k ).set_size( tChainEdges.size(), nullptr );

                    for ( index_t e=0; e<tChainEdges.size(); ++e )
                    {
                        Edge * tEdge = tChainEdges( e );
                        tMidNodes( k )( e ) = this->midnode_of_facet_edge(
                            tEdge->facet( 0 ), tEdge->node( 0 ), tEdge->node( 1 ), tWork );

                        if ( tMidNodes( k )( e ) == nullptr )
                        {
                            const id_t tA = tEdge->node( 0 )->id() ;
                            const id_t tB = tEdge->node( 1 )->id() ;
                            for ( ProtoCurve * tProtoCurve : tProtoCurves )
                            {
                                delete tProtoCurve ;
                            }
                            tReleaseScratch() ;
                            BELFEM_ERROR( false,
                                "The side curve is quadratic, but the facet edge joining nodes %lu and %lu has no midpoint",
                                ( long unsigned int ) tA,
                                ( long unsigned int ) tB );
                        }
                    }
                }
            }

            Cell< Curve * > aCurves( tProtoCurves.size(), nullptr );

            index_t tCount = 0 ;

            for ( ProtoCurve * tProtoCurve : tProtoCurves )
            {
                Cell< Node * > & tEndNodes = tProtoCurve->mNodes ;

                //Create node pairs on the side curve
                Cell< std::pair< Node *, Node * > > tPairs(tProtoCurve->mEdges.size(),std::make_pair( nullptr, nullptr )) ;
                for (uint i = 0 ; i < tProtoCurve->mEdges.size() ; ++i )
                {
                    Edge * tEdge = tProtoCurve->mEdges(i) ;
                    tPairs(i).first = tEdge->node(0) ;
                    tPairs(i).second = tEdge->node(1) ;
                }

                Curve * tCurve = new Curve( ++mMaxCurveID, tTypes( tCount ) );

                this->create_segments( tCurve, tPairs, tMidNodes( tCount ) ) ;

                Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
                this->sort_end_nodes( tCurve, tAdjacency, tEndNodes );

                this->orient_segments( tCurve->is_closed(), tEndNodes, tCurve->segments() ) ;

                this->sort_segments( tCurve->segments() ) ;

                this->collect_nodes( tCurve );

                this->compute_coordinates( tCurve );

                aCurves( tCount++ ) = tCurve ;
                delete tProtoCurve ;
            }

            // restore node indices on mesh, since we have messed with them
            tReleaseScratch() ;

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

            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;

                for ( Facet * tFacet : tFacets )
                {
                    for ( uint e=0; e<tFacet->element()->number_of_edges(); ++e )
                    {
                        tFacet->element()->get_nodes_of_edge( e, tNodes ) ;

                        key_t tA = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 0 )->index() : tNodes( 1 )->index() ;
                        key_t tB = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 1 )->index() : tNodes( 0 )->index() ;
                        key_t tKey = tA * aNumNodes + tB ;

                        aMap( tKey )->increment_facet_counter() ;
                    }
                }
            }

            for ( Edge * tEdge : aEdges )
            {
                tEdge->allocate_facet_container() ;
            }

            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;

                for ( Facet * tFacet : tFacets )
                {
                    for ( uint e=0; e<tFacet->element()->number_of_edges(); ++e )
                    {
                        tFacet->element()->get_nodes_of_edge( e, tNodes ) ;

                        key_t tA = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 0 )->index() : tNodes( 1 )->index() ;
                        key_t tB = tNodes( 0 )->index() > tNodes( 1 )->index() ? tNodes( 1 )->index() : tNodes( 0 )->index() ;
                        key_t tKey = tA * aNumNodes + tB ;

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

            Vector< index_t > tNodeCount( tNumChains, 0 );

            for ( Node * tNode : aNodes )
            {
                ++tNodeCount( tNode->index() );
            }

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

            // Step 1 : create the curve object; owned here until the return
            std::unique_ptr< Curve > tOwner( new Curve(
                aID,
                mMesh->sideset( aSideSetA ),
                mMesh->sideset( aSideSetB ) ) );
            Curve * tCurve = tOwner.get() ;

            ++mMaxCurveID ;

            // Step 2: collects shared corner nodes between side sets
            Cell< Node * > tEndNodes ;
            this->collect_end_nodes_from_intersection( tCurve->sideset_a(), tCurve->sideset_b(), tEndNodes ) ;

            // only nodes on sideset A have been flagged so far
            // let's reset them to maintain a clean dataset
            tCurve->sideset_a()->unflag_all_nodes() ;

            // Step 3: one pair per edge that a facet of A and a facet of B
            // share, and for a quadratic curve its midpoint
            Cell< std::pair< Node *, Node * > > tPairs ;
            Cell< Node * > tMidNodes ;
            this->create_pairs_from_facets( tCurve->sideset_a()->facets(), tEndNodes,
                tCurve->element_type() == ElementType::LINE3, tPairs, tMidNodes,
                & tCurve->sideset_b()->facets() ) ;

            if ( tPairs.size() == 0 )
            {
                mMesh->update_node_indices() ;
                BELFEM_ERROR( false,
                    "Curve %lu: sidesets %lu and %lu share no facet edge, so the curve would be empty",
                    ( long unsigned int ) aID,
                    ( long unsigned int ) aSideSetA,
                    ( long unsigned int ) aSideSetB );
            }

            // Step 5: create the segments
            this->create_segments( tCurve, tPairs, tMidNodes ) ;

            // Step 6: put the end nodes into the appropriate order
            Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
            this->sort_end_nodes( tCurve, tAdjacency, tEndNodes );

            // Step 7: with the nodes sorted, we reorient the segments
            this->orient_segments( tCurve->is_closed(), tEndNodes, tCurve->segments() ) ;

            // Step 8: we can now order the segments in the order they are connected
            this->sort_segments( tCurve->segments() ) ;

            // Step 9: collect and sort the nodes
            this->collect_nodes( tCurve );

            // Step 10 : compute the lengths of the curve
            this->compute_coordinates( tCurve );

            // Step 11 : restore node indices on mesh, since we have messed with them
            mMesh->update_node_indices();

            return tOwner.release() ;
        }

        Curve *
        CurveFactory::from_2d_sidesets( const Vector< id_t > & aSideSets, id_t aID )
        {
            BELFEM_ERROR( comm_rank() == 0, "Curves are not intended to be used in parallel" );

            BELFEM_ERROR( aSideSets.length() > 0, "Curve %lu: no sideset given", ( long unsigned int ) aID );

            // Step 1: Collect the facets
            Cell< Facet * > tFacets ;
            index_t tCount = 0 ;
            for ( id_t tID : aSideSets )
            {
                // an empty sideset would make the type check below depend
                // on the order of the sidesets
                BELFEM_ERROR( mMesh->sideset( tID )->number_of_facets() > 0,
                    "Curve %lu: sideset %lu has no facets",
                    ( long unsigned int ) aID,
                    ( long unsigned int ) tID );

                tCount+=mMesh->sideset( tID )->number_of_facets() ;
            }

            // SideSet::element_type() reads the first facet only, so every
            // facet is checked against it
            const ElementType tType = mMesh->sideset( aSideSets( 0 ) )->element_type() ;

            tFacets.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( id_t tID : aSideSets )
            {
                Cell< Facet * > & tSideSet = mMesh->sideset( tID )->facets() ;
                for ( Facet * tFacet : tSideSet )
                {
                    const ElementType tFacetType = tFacet->element()->type() ;

                    BELFEM_ERROR( tFacetType == ElementType::LINE2 || tFacetType == ElementType::LINE3,
                        "Curve %lu: facet %lu of sideset %lu is of type %s; a 2-D curve needs LINE2 or LINE3 facets",
                        ( long unsigned int ) aID,
                        ( long unsigned int ) tFacet->id(),
                        ( long unsigned int ) tID,
                        to_string( tFacetType ).c_str() );

                    BELFEM_ERROR( tFacetType == tType,
                        "Curve %lu: facet %lu of sideset %lu is of type %s, but sideset %lu is of type %s; "
                        "all sidesets of one curve must have the same element type",
                        ( long unsigned int ) aID,
                        ( long unsigned int ) tFacet->id(),
                        ( long unsigned int ) tID,
                        to_string( tFacetType ).c_str(),
                        ( long unsigned int ) aSideSets( 0 ),
                        to_string( tType ).c_str() );

                    tFacets( tCount++ ) = tFacet ;
                }
            }

            // Step 2: create the curve object; owned here until the return
            std::unique_ptr< Curve > tCurve( new Curve( aID, mMesh->sideset( aSideSets(0) ), nullptr ) );
            Curve * aCurve = tCurve.get() ;
            ++mMaxCurveID ;

            // Step 3: collect the corner nodes
            Cell< Node * > tEndNodes( 2*tFacets.size(), nullptr );
            tCount = 0 ;
            for ( Facet * tFacet : tFacets )
            {
                tEndNodes( tCount++ ) = tFacet->node( 0 );
                tEndNodes( tCount++ ) = tFacet->node( 1 );
            }
            unique( tEndNodes );
            tCount = 0 ;
            for ( Node * tNode : tEndNodes )
            {
                tNode->set_index( tCount++ );
                tNode->unflag() ;
            }

            // Step 4: one pair per facet, and for a quadratic curve its midpoint
            Cell< std::pair< Node *, Node * > > tPairs ;
            Cell< Node * > tMidNodes ;
            this->create_pairs_from_facets( tFacets, tEndNodes,
                aCurve->element_type() == ElementType::LINE3, tPairs, tMidNodes ) ;

            // Step 5: create the segments
            this->create_segments( aCurve, tPairs, tMidNodes ) ;

            // Step 6: put the end nodes into the appropriate order
            Matrix< index_t > tAdjacency( tEndNodes.size(), 2, BELFEM_UINT_MAX );
            this->sort_end_nodes( aCurve, tAdjacency, tEndNodes );

            // Step 7: with the nodes sorted, we reorient the segments
            this->orient_segments( aCurve->is_closed(), tEndNodes, aCurve->segments() ) ;

            // Step 8: we can now order the segments in the order they are connected
            this->sort_segments( aCurve->segments() ) ;

            // Step 9: collect and sort the nodes
            this->collect_nodes( aCurve );

            // Step 10 : compute the lengths of the curve
            this->compute_coordinates( aCurve );

            // Step 11 : restore node indices on mesh, since we have messed with them
            mMesh->update_node_indices();

            return tCurve.release() ;
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
        CurveFactory::create_pairs_from_facets(
                const Cell< Facet * > & aFacets,
                Cell< Node * >        & aEndNodes,
                const bool              aQuadratic,
                Cell< std::pair< Node *, Node * > > & aPairs,
                Cell< Node * > & aMidNodes,
                const Cell< Facet * > * aSharedWith )
        {
            const key_t tNumNodes = aEndNodes.size() ;

            Cell< Node * > tNodes ;

            // refuses two edges that join the same end nodes through
            // different midpoints
            auto tCheckMidNode = [ & ]( const Node * aStored, const Node * aMidNode )
            {
                if ( aQuadratic && aStored != aMidNode )
                {
                    mMesh->update_node_indices() ;
                    BELFEM_ERROR( false,
                        "Two facet edges join nodes %lu and %lu but carry different midpoints",
                        ( long unsigned int ) tNodes( 0 )->id(),
                        ( long unsigned int ) tNodes( 1 )->id() );
                }
            };

            // refuses an edge of a quadratic curve that has no midpoint
            auto tCheckHasMidNode = [ & ]( const Node * aMidNode )
            {
                if ( aQuadratic && aMidNode == nullptr )
                {
                    mMesh->update_node_indices() ;
                    BELFEM_ERROR( false,
                        "The curve is quadratic, but the facet edge joining nodes %lu and %lu has no midpoint",
                        ( long unsigned int ) tNodes( 0 )->id(),
                        ( long unsigned int ) tNodes( 1 )->id() );
                }
            };

            // the edges of aSharedWith between two end nodes, with their midpoints
            Map< key_t, Node * > tSharedMap ;

            if ( aSharedWith != nullptr )
            {
                for ( Facet * tFacet : *aSharedWith )
                {
                    Element * tElement = tFacet->element() ;

                    for ( uint e=0; e<tElement->number_of_edges(); ++e )
                    {
                        tElement->get_nodes_of_edge( e, tNodes );

                        const index_t tA = tNodes( 0 )->index() ;
                        const index_t tB = tNodes( 1 )->index() ;

                        if ( tA == gNoIndex || tB == gNoIndex )
                        {
                            continue ;
                        }

                        const key_t tKey = tA > tB ? tA * tNumNodes + tB : tB * tNumNodes + tA ;

                        Node * tMidNode = tNodes.size() > 2 ? tNodes( 2 ) : nullptr ;

                        if ( tSharedMap.key_exists( tKey ) )
                        {
                            tCheckMidNode( tSharedMap( tKey ), tMidNode );
                        }
                        else
                        {
                            tCheckHasMidNode( tMidNode );
                            tSharedMap[ tKey ] = tMidNode ;
                        }
                    }
                }
            }

            index_t tCount = 0 ;
            for ( Facet * tFacet : aFacets )
            {
                tCount += tFacet->element()->number_of_edges() ;
            }

            Cell< key_t > tKeys ;
            tKeys.reserve( tCount );

            Map< key_t, Node * > tMidNodeMap ;

            for ( Facet * tFacet : aFacets )
            {
                Element * tElement = tFacet->element() ;

                for ( uint e=0; e<tElement->number_of_edges(); ++e )
                {
                    tElement->get_nodes_of_edge( e, tNodes );

                    const index_t tA = tNodes( 0 )->index() ;
                    const index_t tB = tNodes( 1 )->index() ;

                    if ( tA == gNoIndex || tB == gNoIndex )
                    {
                        continue ;
                    }

                    const key_t tKey = tA > tB ? tA * tNumNodes + tB : tB * tNumNodes + tA ;

                    Node * tMidNode = tNodes.size() > 2 ? tNodes( 2 ) : nullptr ;

                    if ( aSharedWith != nullptr )
                    {
                        if ( ! tSharedMap.key_exists( tKey ) )
                        {
                            continue ;
                        }
                        tCheckMidNode( tSharedMap( tKey ), tMidNode );
                    }

                    if ( tMidNodeMap.key_exists( tKey ) )
                    {
                        tCheckMidNode( tMidNodeMap( tKey ), tMidNode );
                    }
                    else
                    {
                        tCheckHasMidNode( tMidNode );
                        tMidNodeMap[ tKey ] = tMidNode ;
                        tKeys.push( tKey );
                    }
                }
            }

            unique( tKeys );

            aPairs.clear() ;
            aPairs.reserve( tKeys.size() );

            if ( aQuadratic )
            {
                aMidNodes.set_size( tKeys.size(), nullptr );
            }

            tCount = 0 ;
            for ( key_t tKey : tKeys )
            {
                key_t tIndexA = tKey % tNumNodes ;
                key_t tIndexB = ( tKey - tIndexA ) / tNumNodes ;

                aPairs.push( std::make_pair( aEndNodes( tIndexA ), aEndNodes( tIndexB ) ) );

                if ( aQuadratic )
                {
                    aMidNodes( tCount ) = tMidNodeMap( tKey ) ;
                }
                ++tCount ;
            }
        }

//------------------------------------------------------------------------------

        Node *
        CurveFactory::midnode_of_facet_edge(
                Facet * aFacet,
                const Node * aA,
                const Node * aB,
                Cell< Node * > & aWork )
        {
            Element * tElement = aFacet->element() ;

            for ( uint e=0; e<tElement->number_of_edges(); ++e )
            {
                tElement->get_nodes_of_edge( e, aWork );

                if ( ( aWork( 0 ) == aA && aWork( 1 ) == aB ) || ( aWork( 0 ) == aB && aWork( 1 ) == aA ) )
                {
                    return aWork.size() == 3 ? aWork( 2 ) : nullptr ;
                }
            }

            return nullptr ;
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
            tCount = 0 ;
            if ( aCurve->element_type() == ElementType::LINE3 )
            {
                // a segment has no node container of its own; the midpoint
                // goes into the wrapped line element
                for ( Segment * tSegment : tSegments )
                {
                    tSegment->element()->insert_node( aMidNodes( tCount++ ), 2 );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::sort_end_nodes( Curve * aCurve, Matrix< index_t > & aAdjacency, Cell< Node * > & aEndNodes )
        {
            Cell< Segment * > & tSegments = aCurve->segments() ;

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

            tNumSegments.fill( 0 );

            // create the adjacency, filled with gNoIndex: next() tests the
            // slots against gNoIndex. Both the constructor and set_size fill
            // with a value typed T, so a 64-bit gNoIndex survives the fill
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

            tCount = 0 ;
            Node * tNode = tStart ;

            Vector< index_t > tNewIndices( tNumCornerNodes, gNoIndex ) ;

            tNewIndices( tNode->index() ) = tCount++ ;

            for ( index_t k=0; k<tNumCornerNodes; ++k )
            {
                tNode->flag();

                tNewIndices( tNode->index() ) = k ;

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

            // a node the walk never reached still carries gNoIndex, which the
            // permutation below would use as a row index into the adjacency.
            // The check belongs after the two forced writes: on an open chain
            // the end node is assigned only by the force, the start by the reset
            index_t tUnassigned = 0 ;
            for ( index_t tIndex : tNewIndices )
            {
                if ( tIndex == gNoIndex )
                {
                    ++tUnassigned ;
                }
            }

            if ( tUnassigned != 0 )
            {
                std::string tSideSets ;
                if ( aCurve->sideset_a() != nullptr )
                {
                    tSideSets += ", sideset A: " + std::to_string( aCurve->sideset_a()->id() ) ;
                }
                if ( aCurve->sideset_b() != nullptr )
                {
                    tSideSets += ", sideset B: " + std::to_string( aCurve->sideset_b()->id() ) ;
                }

                BELFEM_ERROR( false,
                    "Curve %lu is not one connected chain: indices were assigned for %lu of %lu nodes%s",
                    ( luint ) aCurve->id(),
                    ( luint ) ( tNumCornerNodes - tUnassigned ),
                    ( luint ) tNumCornerNodes,
                    tSideSets.c_str() );
            }

            tCount = 0 ;
            for ( Node * tNode : aEndNodes )
            {
                tNode->set_index( tNewIndices( tCount++ ) );
            }

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
        CurveFactory::orient_segments( const bool aClosed, Cell< Node * > & aEndNodes, Cell< Segment * > & aSegments )
        {
            const index_t tLast = aEndNodes.size() - 1 ;

            for ( Segment * tSegment : aSegments )
            {
                const index_t tI = tSegment->node( 0 )->index() ;
                const index_t tJ = tSegment->node( 1 )->index() ;

                // an open chain of two nodes also has a segment { 0, tLast },
                // so the closing segment is told apart by the closed flag
                const bool tClosing = aClosed && ( ( tI == 0 && tJ == tLast ) || ( tI == tLast && tJ == 0 ) );

                if ( tClosing ? tI == 0 : tI > tJ )
                {
                    Node * tSwap = tSegment->node( 0 ) ;

                    // only corners 0 and 1 swap: a LINE3 midpoint stays node 2.
                    // The segment has no insert function of its own
                    tSegment->element()->insert_node(  tSegment->node( 1 ), 0 ) ;
                    tSegment->element()->insert_node(  tSwap, 1 ) ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        CurveFactory::sort_segments( Cell< Segment * > & aSegments )
        {
            for ( Segment * tSegment : aSegments )
            {
                tSegment->set_index( tSegment->node( 0 )->index() );
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

                tSegments( 0 )->node( 0 )->set_index( 0 ) ;
            }

            tNodes.set_size( tCount, nullptr ) ;
            tCount = 0 ;

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