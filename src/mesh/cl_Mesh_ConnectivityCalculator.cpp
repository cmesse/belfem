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
#include "cl_Logger.hpp"

#include "cl_Mesh_ConnectivityCalculator.hpp"

#include "cl_Element_Factory.hpp"
#include "cl_Mesh_Partitioner.hpp"

#include "fn_max.hpp"
#include "op_Graph_Vertex_ID.hpp"
#include "fn_Graph_symrcm.hpp"
#include "cl_HDF5.hpp"
#include "cl_OrderedMap.hpp"
#include "fn_to_master_orientation.hpp"
#include "cl_Timer.hpp"
#include "op_Graph_Vertex_ID.hpp"

namespace belfem
{
    namespace mesh
    {
        ConnectivityCalculator::ConnectivityCalculator( Mesh * aMesh ) :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mMesh( aMesh ),
            mNodes( aMesh->nodes() ),
            mEdges( aMesh->edges() ),
            mFaces( aMesh->faces() ),
            mElements( aMesh->elements() ),
            mFacets( aMesh->facets()),
            mControlPoints( aMesh->control_points() )
        {

        }

        ConnectivityCalculator::~ConnectivityCalculator()
        {

        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_nodes_to_elements()
        {
            for( Node * tNode: mNodes )
            {
                tNode->reset_element_container();
            }

            // loop over all elements
            for( Element * tElement: mElements )
            {
                for( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    tElement->node( k )->increment_element_counter();
                }
            }

            // allocate container for nodes
            for( Node * tNode: mNodes )
            {
                tNode->allocate_element_container();
            }

            for( Element * tElement: mElements )
            {
                for( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    tElement->node( k )->add_element( tElement );
                }
            }

            mMesh->set_connectivity( Connectivity::NodeToElement );
        }

        void
        ConnectivityCalculator::connect_control_points_to_elements()
        {
            for ( ControlPoint * tControlPoint : mMesh->control_points() )
            {
                tControlPoint->reset_element_container();
            }

            for( Element * tElement: mElements )
            {
                for( uint k=0; k<tElement->number_of_control_points(); ++k )
                {
                    tElement->control_point( k )->increment_element_counter();
                }
            }

            for ( ControlPoint * tControlPoint : mMesh->control_points() )
            {
                tControlPoint->allocate_element_container();
            }

            for( Element * tElement: mElements )
            {
                for( uint k=0; k<tElement->number_of_control_points(); ++k )
                {
                    tElement->control_point( k )->add_element( tElement );
                }
            }

            mMesh->set_connectivity( Connectivity::ControlPointToElement );
        }

        void
        ConnectivityCalculator::connect_control_points_to_control_points()
        {
            mMesh->unflag_all_control_points();
            for ( ControlPoint * tControlPoint : mMesh->control_points() )
            {
                tControlPoint->reset_control_point_container();
                tControlPoint->flag();

                for ( uint e=0; e<tControlPoint->number_of_elements(); ++e )
                {
                    Element * tElement = tControlPoint->element( e ) ;
                    for ( uint k=0; k<tControlPoint->element( e )->number_of_control_points(); ++k )
                    {
                        if ( ! tElement->control_point( k )->is_flagged() )
                        {
                            tElement->control_point( k )->flag();
                            tControlPoint->increment_control_point_counter();
                        }
                    }
                }

                tControlPoint->allocate_control_point_container();
                tControlPoint->unflag();
                for ( uint e=0; e<tControlPoint->number_of_elements(); ++e )
                {
                    Element * tElement = tControlPoint->element( e ) ;
                    for ( uint k=0; k<tControlPoint->element( e )->number_of_control_points(); ++k )
                    {
                        if ( tElement->control_point( k )->is_flagged() )
                        {
                            tElement->control_point( k )->unflag();
                            tControlPoint->add_control_point( tElement->control_point( k ) );
                        }
                    }
                }
            }

            mMesh->set_connectivity( Connectivity::ControlPointToControlPoint );

        }

//------------------------------------------------------------------------------

    void
    ConnectivityCalculator::connect_facets_to_elements()
    {
         // loop over all facets
        for ( Facet * tFacet : mFacets )
        {
            // facets arriving from the distributor or a .bfm carry their
            // inherited master/slave links; trust them instead of re-deriving
            // by node identity, which fails where cut relinking made facet
            // and master node lists diverge ( thin shell tapes at cut
            // termini ). No set_master here: its node relink would re-inject
            // the cut duplicates the restoration deliberately removed
            if ( tFacet->has_master() )
            {
                tFacet->element()->allocate_element_container(
                    tFacet->has_slave() ? 2 : 1 );
                tFacet->element()->insert_element( tFacet->master() );
                if ( tFacet->has_slave() )
                {
                    tFacet->element()->insert_element( tFacet->slave() );
                }
                continue ;
            }

            // get number of nodes from this facet
            uint tNumNodes = tFacet->element()->number_of_corner_nodes();

            // count max size of candidates
            uint tCount = 0;
            for ( uint k = 0; k < tNumNodes; ++k )
            {
                Node * tNode = tFacet->node( k );
                tCount += tNode->number_of_elements();
            }

            BELFEM_ASSERT( tCount>0,
                "Could not find candidate elements that could be master of facet %u. Try running 'gmsh -check %s' to check if the mesh is valid.",
                          ( unsigned int ) tFacet->id(), mMesh->path().c_str() );

            // create list of Element candidates
            Cell<Element *> tCandidates( tCount, nullptr );

            // reset counter
            tCount = 0;

            // populate candidates
            for ( uint k = 0; k < tNumNodes; ++k )
            {
                Node * tNode = tFacet->node( k );

                for ( uint e = 0; e < tNode->number_of_elements(); ++e )
                {
                    tCandidates( tCount++ ) = tNode->element( e );
                }
            }

            // make result unique
            unique( tCandidates );


            // loop over all candidates
            Cell<Node *> tNodes;

            Element * tMaster = nullptr;
            Element * tSlave = nullptr;

            uint tMasterFaceIndex = BELFEM_UINT_MAX;
            uint tSlaveFaceIndex = BELFEM_UINT_MAX;

            for ( Element * tElement : tCandidates )
            {
                uint tIndex = compute_facet_index( tFacet, tElement, tNodes );

                if ( tIndex < BELFEM_UINT_MAX )
                {
                    if ( tMaster == nullptr )
                    {
                        tMaster = tElement;
                        tMasterFaceIndex = tIndex ;
                    }
                    else if ( tElement->id() < tMaster->id() )
                    {
                        tSlave = tMaster ;
                        tSlaveFaceIndex = tMasterFaceIndex ;
                        tMaster = tElement ;
                        tMasterFaceIndex = tIndex ;
                        break ;
                    }
                    else
                    {
                        tSlave = tElement;
                        tSlaveFaceIndex = tIndex ;
                        break ;
                    }
                }
            } // end loop over candidates

            BELFEM_ASSERT( tMaster != nullptr,
                "Could not find master for facet %u",
                          ( unsigned int ) tFacet->id() );

            // check if slave is found
            if ( tSlave == nullptr )
            {
                tFacet->set_master( tMaster, tMasterFaceIndex );
                tFacet->element()->allocate_element_container( 1 );
                tFacet->element()->insert_element( tMaster );
            }
            else
            {
                tFacet->set_master( tMaster, tMasterFaceIndex );
                tFacet->set_slave( tSlave, tSlaveFaceIndex );
                tFacet->element()->allocate_element_container( 2 );
                tFacet->element()->insert_element( tMaster );
                tFacet->element()->insert_element( tSlave );
            }
        }

        mMesh->set_connectivity( Connectivity::FacetToElement );
    }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_nodes_to_facets()
        {
            for( Node * tNode: mNodes )
            {
                tNode->reset_facet_container();
            }

            // loop over all elements
            for( Facet * tFacet: mFacets )
            {
                if( tFacet->is_flagged() )
                {
                    for ( uint k = 0; k < tFacet->element()->number_of_nodes(); ++k )
                    {
                        tFacet->element()->node( k )->increment_facet_counter();
                    }
                }
            }

            // allocate container for nodes
            for( Node * tNode: mNodes )
            {
                tNode->allocate_facet_container();
            }

            for( Facet * tFacet: mFacets )
            {
                if( tFacet->is_flagged() )
                {
                    for ( uint k = 0; k < tFacet->element()->number_of_nodes(); ++k )
                    {
                        tFacet->element()->node( k )->add_facet( tFacet );
                    }
                }
            }

            mMesh->set_connectivity( Connectivity::NodeToFacet );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_nodes_to_nodes()
        {
            // container for indices
            Cell< index_t > tIndices ;

            tIndices.reserve( 256 );

            // loop over all nodes
            for( Node * tNode : mNodes )
            {
                // clear the container
                tIndices.clear() ;

                // flag connected nodes
                for ( uint e=0; e<tNode->number_of_elements(); ++e )
                {
                    Element * tElement = tNode->element( e ) ;

                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        if ( tElement->node( k )->id() != tNode->id() )
                        {
                            tIndices.push( tElement->node( k )->index() );
                        }
                    }
                }

                unique( tIndices );

                index_t tNumNodes = tIndices.size()  ;

                if ( tNumNodes > 0 )
                {
                    tNode->allocate_node_container( tNumNodes );
                    index_t i = 0 ;

                    for ( index_t k : tIndices )
                    {
                        tNode->insert_node(  mNodes( k ), i++ ) ;
                    }
                }
            }

            mMesh->set_connectivity( Connectivity::NodeToNode );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_facets_to_facets()
        {
            Cell< index_t > tIndices ;

            // recompute from scratch: drop containers from an earlier pass
            // (a second kernel may recompute after the mesh has changed)
            for ( Facet * tFacet : mFacets )
            {
                tFacet->reset_facet_container() ;
            }

            mMesh->update_facet_indices() ;

            if ( mMesh->number_of_dimensions() == 2 )
            {
                DynamicBitset tBitset( mMesh->number_of_facets() );

                for ( Facet * tFacet : mFacets )
                {
                    tBitset.reset();
                    for ( uint k=0; k< tFacet->number_of_corner_nodes(); ++k )
                    {
                        Node * tOrg = tFacet->node( k )->original() ;
                        for ( uint f=0; f<tOrg->number_of_facets(); ++f )
                        {
                            tBitset.set( tOrg->facet( f )->index() );
                        }
                        for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                        {
                            Node * tDup = tOrg->duplicate( d ) ;
                            for ( uint f=0; f<tDup ->number_of_facets(); ++f )
                            {
                                tBitset.set( tDup->facet( f )->index() );
                            }
                        }
                    }
                    tBitset.reset( tFacet->index() );
                    tBitset.where( tIndices );
                    tFacet->allocate_facet_container( tIndices.size() );
                    for ( index_t f : tIndices )
                    {
                        tFacet->add_facet( mFacets( f ) );
                    }
                }
            }
            else
            {
                index_t tNumKeys = 0 ;
                for ( Facet * tFacet : mFacets )
                {
                    tNumKeys += tFacet->element()->number_of_facets();
                }
                Cell< key128_t > tKeys( tNumKeys, 0 );
                tNumKeys = 0 ;

                Cell< Node * > tNodes ;
                for ( Facet * tFacet : mFacets )
                {

                    Element * tElement = tFacet->element() ;

                    // facets now are actually edges
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tNodes );
                        tKeys( tNumKeys++ ) = this->facet_key2( tNodes );
                    }
                }

                unique( tKeys );
                tNumKeys = tKeys.size() ;
                Graph tEdges( tNumKeys, nullptr );
                for ( index_t k=0; k<tNumKeys; ++k )
                {
                    tEdges( k ) = new graph::Vertex() ;
                    tEdges( k )->set_index( k );
                }


                for ( Facet * tFacet : mFacets )
                {

                    Element * tElement = tFacet->element() ;

                    // facets now are actually edges
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tNodes );
                        tEdges( find_index_in_unique_cell( tKeys, this->facet_key2( tNodes ) ) )->increment_vertex_counter() ;
                    }
                }
                for ( graph::Vertex * tEdge : tEdges )
                {
                    tEdge->init_vertex_container();
                }

                for ( Facet * tFacet : mFacets )
                {

                    Element * tElement = tFacet->element() ;

                    // facets now are actually edges
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tNodes );
                        tEdges( find_index_in_unique_cell( tKeys, this->facet_key2( tNodes ) ) )->insert_vertex( tFacet );
                    }
                }
                Vector< index_t > tCounters( mMesh->number_of_facets(), 0 );


                for ( graph::Vertex * tEdge : tEdges )
                {
                    if ( tEdge->number_of_vertices() < 2 ) continue ;

                    // loop over all facets
                    for ( uint f=0; f<tEdge->number_of_vertices(); ++f )
                    {
                        tCounters( tEdge->vertex( f )->index() ) += tEdge->number_of_vertices() - 1 ;
                    }
                }
                index_t tCount = 0 ;
                for ( Facet * tFacet : mFacets )
                {
                    tFacet->init_vertex_container( tCounters( tCount++ ));
                }

                for ( graph::Vertex * tEdge : tEdges )
                {
                    uint n = tEdge->number_of_vertices() ;
                    if ( n < 2 ) continue ;

                    for ( uint i = 0; i < n; ++i )
                    {
                        graph::Vertex * A = tEdge->vertex( i ) ;
                        for ( uint j = i+1 ; j<n; ++j )
                        {
                            graph::Vertex * B = tEdge->vertex( j ) ;
                            A->insert_vertex( B ) ;
                            B->insert_vertex( A ) ;
                        }
                    }
                }

                Cell< Facet * > tFacets ;
                for ( Facet * tFacet : mFacets )
                {
                    tFacets.set_size( tFacet->number_of_vertices(), nullptr ) ;
                    for ( uint f=0; f<tFacet->number_of_vertices(); ++f )
                    {
                        tFacets( f ) = reinterpret_cast< Facet * >( tFacet->vertex( f ) ) ;
                    }
                    tFacet->reset_vertex_container() ;
                    tFacet->reset_facet_container() ;
                    sort( tFacets.begin(), tFacets.end(), []( Facet * a, Facet * b ) { return a->index() < b->index() ; } ) ;
                    tFacet->allocate_facet_container( tFacets.size() ) ;
                    for ( Facet * tOtherFacet : tFacets )
                    {
                        tFacet->add_facet( tOtherFacet );
                    }
                }
                for ( graph::Vertex * tEdge : tEdges )
                {
                    delete tEdge ;
                }
                tEdges.clear() ;
            }

            mMesh->set_connectivity( Connectivity::FacetToFacet );
        }


//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_nodes_to_edges()
        {
            for( Node * tNode : mNodes )
            {
                tNode->reset_edge_container();
            }

            for( Edge * tEdge : mEdges )
            {
                for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                {
                    tEdge->node( k )->increment_edge_counter();
                }
            }

            for( Node * tNode : mNodes )
            {
                tNode->allocate_edge_container();
            }

            for( Edge * tEdge : mEdges )
            {
                for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                {
                    tEdge->node( k )->add_edge( tEdge );
                }
            }

            mMesh->set_connectivity( Connectivity::NodeToEdge );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_edges_to_elements( Cell< Element * > & aElements )
        {
            for( Edge * tEdge: mEdges )
            {
                tEdge->reset_element_container();
            }

            // loop over all elements
            for( Element * tElement: aElements )
            {
                if( tElement->has_edges() )
                {
                    for ( uint k = 0; k < tElement->number_of_edges(); ++k )
                    {
                        tElement->edge( k )->increment_element_counter();
                    }
                }
            }

            // allocate container for nodes
            for( Edge * tEdge: mEdges )
            {
                tEdge->allocate_element_container();
            }

            for( Element * tElement: aElements )
            {
                if( tElement->has_edges() )
                {
                    for ( uint k = 0; k < tElement->number_of_edges(); ++k )
                    {
                        tElement->edge( k )->add_element( tElement );
                    }
                }
            }

            mMesh->set_connectivity( Connectivity::EdgeToElement );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_edges_to_ghost_facets()
        {
            if ( mMesh->thin_shells().size() == 0 || ! mMesh->edges_exist() ) return ;

            DynamicBitset tBitset( mMesh->number_of_edges() );

            // first we collect all edges that are connected to ghost facets
            for( ThinShell * tThinShell : mMesh->thin_shells() )
            {
                Cell< Facet * > & tFacets = tThinShell->ghost_facets();

                for( Facet * tFacet : tFacets )
                {
                    Element * tMaster = tFacet->master();
                    Element * tSlave = tFacet->slave();

                    BELFEM_ASSERT( tMaster != nullptr, "Ghost facets need a master element" );
                    BELFEM_ASSERT( tSlave != nullptr, "Ghost facets need a slave element" );

                    for ( uint e=0; e<tMaster->number_of_edges(); ++e )
                    {
                        tMaster->edge( e )->increment_facet_counter();
                        tBitset.set( tMaster->edge( e )->index() );
                    }
                    for ( uint e=0; e<tSlave->number_of_edges(); ++e )
                    {
                        tSlave->edge( e )->increment_facet_counter();
                        tBitset.set( tSlave->edge( e )->index() );
                    }
                }
            }

            Cell< index_t > tIndices ;
            tBitset.where( tIndices );

            Cell< Edge * > & tEdges = mMesh->edges();

            for( index_t k : tIndices )
            {
                tEdges( k )->allocate_facet_container();
            }

            for( ThinShell * tThinShell : mMesh->thin_shells() )
            {
                Cell< Facet * > & tFacets = tThinShell->ghost_facets();

                for( Facet * tFacet : tFacets )
                {
                    Element * tMaster = tFacet->master();
                    Element * tSlave = tFacet->slave();

                    BELFEM_ASSERT( tMaster != nullptr, "Ghost facets need a master element" );
                    BELFEM_ASSERT( tSlave != nullptr, "Ghost facets need a slave element" );

                    for ( uint e=0; e<tMaster->number_of_edges(); ++e )
                    {
                        tMaster->edge( e )->add_facet( tFacet );
                    }
                    for ( uint e=0; e<tSlave->number_of_edges(); ++e )
                    {
                        tSlave->edge( e )->add_facet( tFacet );
                    }
                }
            }

            // once we are done, we make sure that the facets are unique
            Cell< Facet * > tFacets ;

            for( index_t k : tIndices )
            {
                Edge * tEdge = tEdges( k ) ;

                tFacets.set_size( tEdge->number_of_facets(), nullptr ) ;
                for ( uint f=0; f<tEdge->number_of_facets(); ++f )
                {
                    tFacets( f ) = tEdge->facet( f ) ;
                }
                tEdge->reset_facet_container() ;
                unique( tFacets );
                tEdge->allocate_facet_container( tFacets.size() ) ;
                for ( Facet * tFacet : tFacets )
                {
                    tEdge->add_facet( tFacet );
                }
            }

            mMesh->set_connectivity( Connectivity::EdgeToFacet );
        }

        void ConnectivityCalculator::connect_faces_to_ghost_facets()
        {
            if ( mMesh->thin_shells().size() == 0 || ! mMesh->faces_exist() ) return ;

            for( ThinShell * tThinShell : mMesh->thin_shells() )
            {
                Cell< Facet * > & tFacets = tThinShell->ghost_facets();

                for( Facet * tFacet : tFacets )
                {
                    Face * tMaster = tFacet->master()->face( tFacet->index_on_master() );
                    Face * tSlave = tFacet->slave()->face( tFacet->index_on_slave() );

                    tMaster->allocate_facet_container( 1 );
                    tMaster->add_facet( tFacet );
                    tSlave->allocate_facet_container( 1 );
                    tSlave->add_facet( tFacet );
                }
            }
            mMesh->set_connectivity( Connectivity::FaceToFacet );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_edges_to_edges()
        {
            Cell< Edge * > tEdges ;
            mMesh->unflag_all_edges();

            for( Edge * tEdge : mEdges )
            {
                tEdge->reset_edge_container();
                tEdge->flag();
                uint tCount = 0 ;
                for ( uint e=0; e<tEdge->number_of_elements(); ++e )
                {
                    Element * tElement = tEdge->element( e );
                    for ( uint k=0; k<tElement->number_of_edges(); ++k )
                    {
                        Edge * tOther = tElement->edge( k );
                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            ++tCount ;
                        }
                    }
                }
                // ghost facets: traverse master/slave elements to discover
                // cross-layer edge neighbors (ghost facet placeholder has no edges)
                for ( uint f=0; f<tEdge->number_of_facets(); ++f )
                {
                    Facet * tFacet = tEdge->facet( f );
                    for ( uint k=0; k<tFacet->master()->number_of_edges(); ++k )
                    {
                        Edge * tOther = tFacet->master()->edge( k );
                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            ++tCount ;
                        }
                    }
                    for ( uint k=0; k<tFacet->slave()->number_of_edges(); ++k )
                    {
                        Edge * tOther = tFacet->slave()->edge( k );
                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            ++tCount ;
                        }
                    }
                }

                if ( tCount == 0 ) continue ;
                tEdges.set_size( tCount, nullptr );
                tCount = 0 ;
                tEdge->unflag();

                for ( uint e=0; e<tEdge->number_of_elements(); ++e )
                {
                    Element * tElement = tEdge->element( e );
                    for ( uint k=0; k<tElement->number_of_edges(); ++k )
                    {
                        Edge * tOther = tElement->edge( k );
                        if ( tOther->is_flagged() )
                        {
                            tOther->unflag();
                            tEdges( tCount++ ) = tOther ;
                        }
                    }
                }

                for ( uint f=0; f<tEdge->number_of_facets(); ++f )
                {
                    Facet * tFacet = tEdge->facet( f );
                    for ( uint k=0; k<tFacet->master()->number_of_edges(); ++k )
                    {
                        Edge * tOther = tFacet->master()->edge( k );
                        if ( tOther->is_flagged() )
                        {
                            tOther->unflag();
                            tEdges( tCount++ ) = tOther ;
                        }
                    }
                    for ( uint k=0; k<tFacet->slave()->number_of_edges(); ++k )
                    {
                        Edge * tOther = tFacet->slave()->edge( k );
                        if ( tOther->is_flagged() )
                        {
                            tOther->unflag();
                            tEdges( tCount++ ) = tOther ;
                        }
                    }
                }

                sort( tEdges, opVertexID );
                tEdge->allocate_edge_container( tEdges.size() );
                for ( Edge * tOtherEdge : tEdges )
                {
                    tEdge->add_edge( tOtherEdge );
                }
            }

            mMesh->set_connectivity( Connectivity::EdgeToEdge );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_faces_to_edges_and_edges_to_faces()
        {
            if( mMesh->number_of_dimensions() == 2 )
            {
                return ;
            }

            mMesh->unflag_all_nodes() ;
            mMesh->unflag_all_edges();
            mMesh->unflag_all_faces();
            mMesh->unflag_all_elements() ;

            Cell< Edge * > tEdgesOnMaster ;
            Cell< Edge * > tEdgesOnSlave ;

            for( Edge * tEdge : mEdges )
            {
                tEdge->reset_face_container();
            }
            for( Face * tFace : mFaces )
            {
                tFace->reset_edge_container();
            }

            for( Face * tFace : mFaces )
            {
                tEdgesOnMaster.clear();

                if ( tFace->master() != nullptr )
                {
                    if ( tFace->master()->has_edges() )
                    {
                        tFace->master()->get_edges_of_facet( tFace->index_on_master(),
                                         tEdgesOnMaster );
                    }
                }
                else if ( tFace->slave() != nullptr && tEdgesOnMaster.size() >= 0 )
                {
                    tFace->slave()->get_edges_of_facet( tFace->index_on_slave(), tEdgesOnSlave );
                    to_master_orientation( tFace, tEdgesOnSlave, tEdgesOnMaster );
                }


                tFace->reset_edge_container();
                if ( tEdgesOnMaster.size() > 0 )
                {
                    tFace->allocate_edge_container( tEdgesOnMaster.size() ) ;

                    uint tCount = 0 ;
                    for( Edge * tEdge : tEdgesOnMaster )
                    {
                        tFace->insert_edge( tEdge, tCount++ );
                        tEdge->increment_face_counter();
                    }
                }
            }

            for( Edge * tEdge : mEdges )
            {
                tEdge->allocate_face_container();
            }

            for( Face * tFace : mFaces )
            {
                for ( uint e=0; e<tFace->number_of_edges(); ++e )
                {
                    tFace->edge( e )->add_face( tFace );
                }
            }

            mMesh->set_connectivity( Connectivity::FaceToEdge );
            mMesh->set_connectivity( Connectivity::EdgeToFace );
        }

        void
        ConnectivityCalculator::connect_faces_to_faces()
        {
            Cell< Face * > tFaces ;
            for ( Face * tFace : mFaces )
            {
                tFace->reset_face_container();
                index_t tCount = 0 ;
                tFace->flag();

                for ( uint e=0; e<tFace->number_of_edges(); ++e )
                {
                    for ( uint f=0; f<tFace->edge( e )->number_of_faces(); ++f )
                    {
                        Face * tOther = tFace->edge( e )->face( f );
                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            ++tCount ;
                        }
                    }
                }

                // for ghost facets
                if ( tFace->number_of_facets() == 1 )
                {
                    Element * tMaster = tFace->master();
                    for ( uint f=0; f<tMaster->number_of_faces(); ++f )
                    {
                        Face * tOther = tMaster->face( f );
                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            ++tCount ;
                        }
                    }
                    Element * tSlave = tFace->slave();
                    for ( uint f=0; f<tSlave->number_of_faces(); ++f )
                    {
                        Face * tOther = tSlave->face( f );
                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            ++tCount ;
                        }
                    }
                }


                if ( tCount == 0 ) continue ;


                tFaces.set_size( tCount, nullptr );
                tCount = 0 ;
                tFace->unflag();

                for ( uint e=0; e<tFace->number_of_edges(); ++e )
                {
                    Edge * tEdge = tFace->edge( e );
                    for ( uint f=0; f<tEdge->number_of_faces(); ++f )
                    {
                        Face * tOther = tFace->edge( e )->face( f );
                        if ( tOther->is_flagged() )
                        {
                            tOther->unflag();
                            tFaces( tCount++ ) = tOther ;
                        }
                    }
                }

                if ( tFace->number_of_facets() == 1 )
                {
                    Element * tMaster = tFace->master();
                    for ( uint f=0; f<tMaster->number_of_faces(); ++f )
                    {
                        Face * tOther = tMaster->face( f );
                        if ( tOther->is_flagged() )
                        {
                            tOther->unflag();
                            tFaces( tCount++ ) = tOther ;
                        }
                    }
                    Element * tSlave = tFace->slave();
                    for ( uint f=0; f<tSlave->number_of_faces(); ++f )
                    {
                        Face * tOther = tSlave->face( f );
                        if ( tOther->is_flagged() )
                        {
                            tOther->unflag();
                            tFaces( tCount++ ) = tOther ;
                        }
                    }
                }

                sort( tFaces, opVertexID );
                tFace->allocate_face_container( tFaces.size() ) ;
                for ( Face * tNeighbor : tFaces )
                {
                    tFace->add_face( tNeighbor );
                }
            }

            mMesh->set_connectivity( Connectivity::FaceToFace );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::connect_elements_to_elements()
        {
            Timer tTimer ;

            // computing element-to-element connectivity is expensive. We skip this if it is already set
            if ( mMesh->test_connectivity( Connectivity::ElementToElement )
                 && ( mMesh->test_connectivity( Connectivity::TsElementToTsElement )
                     || mMesh->thin_shells().size() == 0 )
            ) return;



            DynamicBitset tBitset( mMesh->number_of_elements() );
            mMesh->update_node_indices();
            mMesh->update_element_indices();

            Cell< Element * > & tElements = mMesh->elements() ;
            if ( ! mMesh->test_connectivity( Connectivity::ElementToElement ) )
            {
                if ( mCommRank == 0 )
                {
                    message( InfoLevel::Default, "    computing element connectivities ..." );
                }

                for ( Element * tElement : tElements )
                {
                    tElement->flag();
                }
                this->connect_elements_to_elements_sub( tElements );

                mMesh->set_connectivity( Connectivity::ElementToElement );
            }

            if ( mMesh->thin_shells().size() > 0 && ! mMesh->test_connectivity( Connectivity::TsElementToTsElement ) )
            {
                if ( mCommRank == 0 )
                {
                    message( InfoLevel::Default, "    computing element connectivities ( thin shells ) ..." );
                }

                // select elements that sit on thin shells
                mMesh->unflag_all_elements() ;
                for ( ThinShell * tShell : mMesh->thin_shells() )
                {
                    for ( Block * tBlock : tShell->blocks() )
                    {
                        for ( Element * tElement : tBlock->elements() )
                        {
                            tElement->flag();
                            tBitset.set( tElement->index() );
                        }
                    }
                }

                Cell< index_t > tIndices ;
                tBitset.where( tIndices );
                Cell< Element * > tTsElements( tIndices.size() , nullptr ) ;

                index_t tCount = 0 ;
                for ( index_t e : tIndices )
                {
                    Element * tElement = tElements( e ) ;
                    tElement->set_index( tCount );
                    tTsElements( tCount++ ) = tElement ;
                }
                this->connect_elements_to_elements_sub( tTsElements );
                mMesh->update_element_indices();
                mMesh->set_connectivity( Connectivity::TsElementToTsElement );
            }

            if ( mCommRank == 0 )
            {
                message( InfoLevel::Default, "    ... done in %u ms.", tTimer.stop() );
            }



        }

        void
        ConnectivityCalculator::connect_elements_to_elements_sub(
            Cell< Element * > & aElements )
        {
            index_t tNumKeys = 0 ;
            index_t tNumElements = 0  ;
            uint tDim = mMesh->number_of_dimensions()  ;

            mKey =  tDim == 2 ? &
                    ConnectivityCalculator::facet_key2
                : & ConnectivityCalculator::facet_key3;

            // count keys
            for ( Element * tElement : aElements )
            {
                if ( tElement->dimension() == tDim )
                {
                    tElement->set_index( tNumElements++ ) ;
                    tNumKeys += tElement->number_of_facets() ;
                }
                else
                {
                    tElement->set_index( gNoIndex ) ;
                }
            }
            Cell< key128_t > tKeys( tNumKeys, 0 ) ;

            tNumKeys = 0 ;
            Cell< Node * > tNodes ;
            for ( Element * tElement : aElements )
            {
                if ( tElement->dimension() == tDim )
                {
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tNodes ) ;
                        tKeys( tNumKeys++ ) = (this->*mKey)( tNodes ) ;

                    }
                }
            }
            unique( tKeys ) ;
            tNumKeys = tKeys.size() ;
            Matrix< index_t > tConnectivity( 2, tNumKeys, gNoIndex ) ;
            for ( Element * tElement : aElements )
            {
                if ( tElement->dimension() == tDim )
                {
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tNodes ) ;
                        index_t tIndex
                            = find_index_in_unique_cell( tKeys,  (this->*mKey)( tNodes ) );
                        if ( tConnectivity( 0, tIndex ) == gNoIndex )
                        {
                            tConnectivity( 0, tIndex ) = tElement->index() ;
                        }
                        else
                        {
                            tConnectivity( 1, tIndex ) = tElement->index() ;
                        }
                    }
                }
            }

            Cell< index_t > tCounters( tNumElements, 0 ) ;

            for ( index_t k = 0 ; k<tNumKeys; ++k )
            {
                index_t A = tConnectivity( 0, k ) ;
                index_t B = tConnectivity( 1, k ) ;
                if ( A != gNoIndex && B != gNoIndex )
                {
                    if ( aElements( A )->is_flagged() && aElements( B )->is_flagged() )
                    {
                        ++tCounters( A ) ;
                        ++tCounters( B ) ;
                    }
                }
            }
            tNumElements = 0 ;
            for ( Element * tElement : aElements )
            {
                index_t tCount = tCounters( tNumElements++ ) ;
                tElement->reset_element_container();
                if ( tCount > 0 )
                {
                    tElement->allocate_element_container( tCount );
                }
            }
            tCounters.clear() ;
            for ( index_t k = 0 ; k<tNumKeys; ++k )
            {
                index_t A = tConnectivity( 0, k ) ;
                index_t B = tConnectivity( 1, k ) ;


                if ( A != gNoIndex && B != gNoIndex )
                {
                    if ( aElements( A )->is_flagged() && aElements( B )->is_flagged() )
                    {
                        aElements( A )->insert_element( aElements( B ) ) ;
                        aElements( B )->insert_element( aElements( A ) ) ;
                    }
                }
            }
            Cell< Element * > tElements ;
            mMesh->update_element_indices() ;
            for ( Element * tElement : aElements )
            {
                tElements.set_size( tElement->number_of_elements(), nullptr ) ;
                for ( uint e=0; e<tElement->number_of_elements(); ++e )
                {
                    tElements( e ) = tElement->element( e ) ;
                }
                sort( tElements.begin(), tElements.end(),
                    []( Element * a, Element * b ) { return a->index() < b->index(); } ) ;
                tElement->allocate_element_container( tElements.size() ) ;
                for ( Element * tOtherElement : tElements )
                {
                    tElement->insert_element( tOtherElement ) ;
                }
            }
        }

        void
        ConnectivityCalculator::connect_thin_shells_to_thin_shells()
        {
            if ( mMesh->thin_shells().size() == 0 ) return ;
            if ( mMesh->test_connectivity( Connectivity::ShellToShell ) ) return ;

            mMesh->update_element_map() ;

            // creating a temporary mesh
            Mesh * tMesh = mMesh->extract_thin_shell_mesh();

            // with the mesh created, we can now set neighbors and elements
            Cell< Element * > & tElements = tMesh->elements() ;

            for ( Element * tDup : tElements )
            {
                // get element on original mesh
                Element * tOrg = mMesh->element( tDup->id() ) ;

                // link elements
                tOrg->allocate_element_container( tDup->number_of_elements() );
                for ( uint e=0; e<tDup->number_of_elements(); ++e )
                {
                    tOrg->insert_element( mMesh->element( tDup->element( e )->id() ) );
                }
            }

            if ( comm_rank() == 0 )
            {

                // neighbors will be needed later to create the thermal mesh
                tMesh->populate_element_neighbors() ;

                for ( Element * tDup : tElements )
                {
                    // get element on original mesh
                    Element * tOrg = mMesh->element( tDup->id() ) ;

                    tOrg->allocate_neighbor_container();

                    for ( uint f=0; f<tOrg->number_of_facets(); ++f )
                    {
                        if ( tDup->neighbor( f ) != nullptr )
                        {
                            tOrg->insert_neighbor( mMesh->element( tDup->neighbor( f )->id() ), f );
                        }
                    }
                }
            }

            delete tMesh ;
            mMesh->set_connectivity( Connectivity::ShellToShell );
        }

//------------------------------------------------------------------------------

        void
        ConnectivityCalculator::save_connectivity_data( const string & aPath )
        {
            if ( mCommRank != 0 ) return ;

            // count memory
            index_t tCount = 0 ;
            sort( mElements, opVertexID );

            for ( Element * tElement : mElements )
            {
                tCount += 2 + tElement->number_of_elements() ;
            }

            Vector< id_t > tData( tCount );

            tCount = 0 ;
            Cell< id_t > tIDs ;

            for ( Element * tElement : mElements )
            {
                tData( tCount++ ) = tElement->id() ;
                tData( tCount++ ) = tElement->number_of_elements();

                tIDs.set_size( tElement->number_of_elements(), 0 );

                for ( uint e=0; e<tElement->number_of_elements(); ++e )
                {
                    tIDs( e ) = tElement->element( e )->id();
                }
                sort( tIDs );
                for ( id_t tID : tIDs )
                {
                    tData( tCount++ ) = tID ;
                }
            }

            HDF5 tFile( aPath, FileMode::NEW );
            tFile.save_data( "data",  tData );
            tFile.close();
        }

    }
}
