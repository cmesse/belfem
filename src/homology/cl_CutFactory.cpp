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

#include "cl_Logger.hpp"
#include "cl_Timer.hpp"
#include "cl_CutFactory.hpp"

#include "commtools.hpp"
#include "fn_unique.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IwgFactory.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_dot.hpp"
#include "cl_BeltedTree.hpp"
#include "cl_CutProcessor.hpp"

#include "cl_DynamicBitset.hpp"
#include "fn_combine.hpp"
#include "cl_InterfaceProcessor.hpp"
#include "cl_CurveFactory.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_Graph_find_connected_partitions.hpp"
#include "op_Graph_Vertex_Owner.hpp"
#include "cl_Queue.hpp"
#include "fn_Mesh_symrcm_nodes.hpp"

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        CutFactory::CutFactory(
            Mesh * aMesh,
            Topology * aTopology,
            Cell< Protoshell * >  & aProtoshells,
            const CutAlgorithm aAlgorithm,
            const bool aUseEnrichment ) :
                mCommRank( comm_rank()),
                mMesh( aMesh ),
                mTopology( aTopology ),
                mThinShellSidesets( aTopology->groups( DomainType::ThinShell ) ),
                mAbstractNodes( aMesh->abstract_nodes() ),
                mOrphanedNodes( aMesh->orphaned_nodes() ),
                mProtoshells( aProtoshells ),
                mAlgorithm( aAlgorithm ),
                mUseEnrichment( aUseEnrichment )
        {

        }

//-----------------------------------------------------------------------------

        CutFactory::~CutFactory()
        {
            for (  SideSet * tSideSet: mTemporaryThinShellSidesets  )
            {
                Cell< Facet * > & tFacets = tSideSet->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    delete tFacet ;
                }
                delete tSideSet ;
            }
            for ( auto tPair : mIndicesOfOriginalTerminalNodes )
            {
                delete tPair.second ;
            }

            if( mCohomology != nullptr )
            {
                delete mCohomology ;
            }
            if( mRelativeHomology != nullptr )
            {
                delete mRelativeHomology ;
            }
            if( mSimplicialComplex != nullptr )
            {
                delete mSimplicialComplex ;
            }

        }

//-----------------------------------------------------------------------------

        void
        CutFactory::set_terminals( const Cell<Cell< id_t >> & aTerminals, const Cell<id_t> & aThinShellIndices )
        {
            mTerminals = aTerminals ;
            mThinShellIndices = aThinShellIndices ;
        }

//-----------------------------------------------------------------------------

        void CutFactory::set_periodicity( PeriodicityFactory *aPeriodicFactory )
        {
            mPeriodicFactory = aPeriodicFactory ;
        }


//-----------------------------------------------------------------------------

        void
        CutFactory::run()
        {
            if( mCommRank == mMesh->master() )
            {
                message( InfoLevel::Default, "    Creating cuts ...");
                mMesh->unfinalize() ;
                mMesh->finalize() ;
                mMesh->create_edges() ;
                mMesh->finalize_edges();

                if( mMesh->number_of_dimensions() == 3  )
                {
                    mMesh->create_faces() ;
                    mMesh->finalize_faces() ;
                }

                //Create periodicity with full mesh
                if ( mPeriodicFactory != nullptr )
                {
                    mMesh->set_periodicity( mPeriodicFactory->create_periodicity() ) ;
                    mMesh->periodicity()->backup_node_pairs();
                }

                // assign edges to segments
                for ( Curve * tCurve : mMesh->curves() )
                {
                    tCurve->assign_edges();
                }
            }

            if ( mUsePoissonInsteadOfRcm )
            {
                this->compute_poisson_problem();
            }
            else
            {
                this->compute_rcm_problem();
            }

            if( mCommRank == mMesh->master() )
            {
                this->compute_cohomologies() ;

                Timer tTimer;

                gLog.message( InfoLevel::Detailed, "    creating thin cuts ... " );


                this->compute_thin_cuts_and_duplicate_interface_nodes() ;
                this->restore_thin_shell_sidesets();

                gLog.message( InfoLevel::Detailed, "    ... time %.3f s\n", tTimer.stop() * 1e-3 );

                // before we delete the meshes, we write the element neighborhoods
                // which will be needed for the thin shell factory
                this->compute_element_adjacencies();


                mMesh->unfinalize() ;
                mMesh->reset_edges() ;
                mMesh->reset_faces() ;
                mMesh->finalize() ;

                //mMesh->save("cut.exo");


            }
            comm_barrier() ;

        }

//-----------------------------------------------------------------------------

        void
        CutFactory::compute_rcm_problem()
        {
            Timer tTimer ;

            if( mCommRank == 0 )
            {
                gLog.message( InfoLevel::Detailed, "    solving Reverse Cuthill McKee ... " );

                Cell< Node * > & tNodes = mMesh->nodes();
                symrcm( tNodes );

                Vector< real > & tPhi = mMesh->field_exists("phi") ? mMesh->field_data("phi") : mMesh->create_field("phi");
                tPhi.fill( BELFEM_QUIET_NAN );

                for ( Node * tNode : tNodes )
                {
                    tNode->reset_vertex_container();
                    tPhi( tNode->index() ) = static_cast< real >( tNode->index() );
                }

                gLog.message( InfoLevel::Detailed, "    ... solving time : %.3f s\n", tTimer.stop() * 1e-3 );
            }
            comm_barrier() ;
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::compute_cohomologies()
        {
            mMesh->update_node_indices();
            mMesh->update_edge_indices();
            mMesh->update_face_indices();

            Homology * tSuggestedHomology = nullptr;

            if( mMesh->number_of_dimensions() == 2 )
            {
                if( mSuggestHomologies )
                {
                    tSuggestedHomology = new Homology(
                            mMesh,
                            mTerminals,
                            mThinShellIndices) ;

                    // Fix the sense of each Amperian loop, which fixes the
                    // sign of the imposed current ( see reorient_generators )
                    tSuggestedHomology->reorient_generators() ;
                }
                mMesh->unflag_everything() ;
                for( id_t tID : mTopology->phi_block_ids() )
                {
                    mMesh->block( tID )->flag_elements() ;

                    for( Element * tElement : mMesh->block( tID )->elements() )
                    {
                        tElement->flag_corner_nodes() ;
                        tElement->flag_edges() ;
                    }
                }
            }
            else
            {
                if( mSuggestHomologies )
                {
                    tSuggestedHomology = new Homology(
                            mMesh,
                            mTerminals,
                            mThinShellIndices) ;

                    // Fix the sense of each Amperian loop, which fixes the
                    // sign of the imposed current ( see reorient_generators )
                    tSuggestedHomology->reorient_generators() ;
                }

                mMesh->unflag_everything() ;
                for( id_t tID : mTopology->phi_block_ids() )
                {
                    mMesh->block( tID )->flag_elements() ;

                    for( mesh::Element * tElement : mMesh->block( tID )->elements() )
                    {
                        tElement->flag_corner_nodes() ;
                        tElement->flag_edges() ;
                        tElement->flag_faces() ;
                    }
                }
            }

            this->unflag_symmetry_sidesets() ;
            mSimplicialComplex = new SimplicialComplex( mMesh, true );

            switch( mAlgorithm )
            {

                case( CutAlgorithm::Pellikka ) :
                {
                    if ( ! mSuggestHomologies )
                    {
                        mSimplicialComplex->reduce_complexPellikka();
                    }
                    mSimplicialComplex->coreduce_complexPellikka();

                    Timer tTimer ;
                    gLog.message( InfoLevel::Detailed, "    Computing Smith normal form ... " );
                    mCohomology = new Cohomology( mSimplicialComplex, mMesh );
                    gLog.message( InfoLevel::Detailed, "    ... done. Computation time : %.3f s\n", tTimer.stop() * 1e-3 );


                    break ;
                }
                case( CutAlgorithm::CCR ) :
                {

                    if ( ! mSuggestHomologies )
                    {
                        mSimplicialComplex->reduce_complexCCR();
                    }

                    mSimplicialComplex->coreduce_complexCCR();

                    Timer tTimer ;
                    gLog.message( InfoLevel::Detailed, "    Computing Smith normal for ... " );

                    mCohomology = new Cohomology( mSimplicialComplex, mMesh );

                    gLog.message( InfoLevel::Detailed, "    ... done. Creation time : %.3f s\n", tTimer.stop() * 1e-3 );

                    break;
                }
                case( CutAlgorithm::BeltedTree ) :
                {
                    BELFEM_ERROR(mSuggestHomologies, "Belted Tree algorithm requires suggesting homology");

                    BeltedTree * tBTree = new mesh::BeltedTree(mMesh,  mSimplicialComplex, tSuggestedHomology->get_Generators()(1)) ;

                    Timer tTimer ;
                    gLog.message( InfoLevel::Detailed, "    Computing Smith normal for ... " );

                    mCohomology = new Cohomology( mSimplicialComplex, mMesh, tBTree);

                    gLog.message( InfoLevel::Detailed, "    ... done. Creation time : %.3f s\n", tTimer.stop() * 1e-3 );

                    delete tBTree;
                    break ;
                }
                case( CutAlgorithm::PellikkaGeneralized ) :
                {
                    if ( ! mSuggestHomologies )
                    {
                        mSimplicialComplex->reduce_complexPellikkaGeneralized();
                    }
                    mSimplicialComplex->coreduce_complexPellikkaGeneralized();

                    Timer tTimer ;
                    gLog.message( InfoLevel::Detailed, "    Computing Smith normal form ... " );
                    mCohomology = new Cohomology( mSimplicialComplex, mMesh );
                    gLog.message( InfoLevel::Detailed, "    ... done. Creation time : %.3f s\n", tTimer.stop() * 1e-3 );

                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Undefined Cut algorithm.");
                }
            } // end switch


            // for debugging

            //this->write_debug_cohomology( tSuggestedHomology ) ;

            //Isolating the cuts
            if ( mSuggestHomologies )
            {
                //Isolating the cuts
                mCohomology->updatekGeneratorsFromHomology( tSuggestedHomology->get_Generators()( 1 ), 1 );

                mCohomology->clean();
            }

            // for debugging
            //this->write_debug_cohomology( tSuggestedHomology ) ;

            // tidy up
            if( tSuggestedHomology != nullptr )
            {
                delete tSuggestedHomology;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::write_debug_cohomology( Homology * aSuggestedHomology )
        {

            Mesh * tMeshEdge = new Mesh( mMesh->number_of_dimensions() , 0, false );
            Cell< Node * > & tNodes = tMeshEdge->nodes();
            Map< id_t, Node * > tNodeMap ;
            id_t tNumNodes = mMesh->number_of_nodes();
            id_t tNumEdges = mMesh->number_of_edges();

            // allocate node memory
            tNodes.set_size( tNumNodes, nullptr );
            for ( uint k=0; k<tNumNodes; ++k)
            {
                // grab data
                id_t id = mMesh->nodes()(k)->id();
                real x = mMesh->nodes()(k)->x();
                real y = mMesh->nodes()(k)->y();
                real z = mMesh->nodes()(k)->z();

                // create a new node
                Node * tNode =  new Node( id, x, y, z );

                // add node to container
                tNodes( k ) = tNode ;

                // add node to map
                tNodeMap[ id ] = tNode ;

            }

            Block * tBlock = new Block( 1, tNumEdges );
            tBlock->label() = "Edges";
            // create element factory
            ElementFactory tElementFactory ;

            // grab element container
            Cell< mesh::Element * > & tEdges = tBlock->elements();

            //Cell< mesh::Edge * > & tEdges = tMesh->edges();
            tEdges.set_size( tNumEdges, nullptr );
            for ( uint e=0; e<tNumEdges; ++e)
            {
                id_t edge_id = mMesh->edges()(e)->id();
                id_t node_id_1 = mMesh->edges()(e)->node(0)->id();
                id_t node_id_2 = mMesh->edges()(e)->node(1)->id();

                Element * tEdge = tElementFactory.create_element( ElementType::LINE2, edge_id );
                tEdge->insert_node( tNodeMap(node_id_1), 0  );

                tEdge->insert_node( tNodeMap( node_id_2 ), 1 );
                tEdges( e ) = tEdge ;

            }

            // add the block to the mesh
            tMeshEdge->blocks().push( tBlock );
            tMeshEdge->finalize();

            if (mSuggestHomologies)
            {
                aSuggestedHomology->create_kGeneratorsField(1, tMeshEdge, "Homology");
            }
            mCohomology->create_kGeneratorsField(1, tMeshEdge, "Cohomology");
            tMeshEdge->save( "Homology.exo" );
            delete tMeshEdge ;

        }

//-----------------------------------------------------------------------------

        void
        CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()
        {
            Timer tTimer ;

            // need to backup this number since the CutProcessor will change it
            index_t tNumOriginalSidesets = mMesh->number_of_sidesets() ;

            // create the processor
            CutProcessor tProc( mMesh, mCohomology, mTopology->phi_block_ids(),
                                mTopology->non_phi_block_ids(),
                                mTopology->phi_interface_ids(),
                                mTopology->phi_boundary_ids(),
                                mTopology->phi_periodic_ids() );

            mAbstractNodes.vector_data() = std::move( tProc.abstract_nodes().vector_data() );

            // this must be done before we call the interface processor
            this->link_node_duplicates_and_originals() ;

            // duplicate interface nodes
            InterfaceProcessor tInterfaceProcessor( mMesh, mTopology, mAbstractNodes, tNumOriginalSidesets, tProc.max_node_id() );

            this->collect_orphan_nodes() ;

#if !defined( NDEBUG ) || defined( DEBUG )
            tProc.save_debug_meshes();
#endif

            gLog.message( InfoLevel::Detailed, "    ... done. Creation time : %.3f s\n", tTimer.stop() * 1e-3 );
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::restore_thin_shell_sidesets()
        {
            id_t tID = mMesh->max_block_and_sideset_id() ;

            for ( SideSet * tTemp : mTemporaryThinShellSidesets )
            {
                SideSet * tSideSet = mMesh->sideset( tTemp->id() );

                if ( mUseEnrichment )
                {
                    SideSet * tInterface = new SideSet( ++tID, 0);
                    tInterface->facets() = std::move( tSideSet->facets() );
                    tInterface->set_facet_counter( tInterface->number_of_facets() );
                    tInterface->set_domain_type( DomainType::EnrichedInterface );
                    tInterface->hide( true );
                    mMesh->add_sideset( tInterface );
                    std::cout << "#add sideset " << tInterface->id() << std::endl;
                }
                else
                {
                    // delete temporary facets on mesh
                    for ( Facet * tFacet : tSideSet->facets() )
                    {
                        delete tFacet;
                    }
                }
                tSideSet->reset_facet_container();
                tSideSet->facets() = std::move( tTemp->facets() );
                tSideSet->set_facet_counter( tSideSet->number_of_facets() );

                tTemp->reset_facet_container();
                delete tTemp;
            }
            mTemporaryThinShellSidesets.clear();
        }

//-----------------------------------------------------------------------------

        ElementType
        CutFactory::check_element_types()
        {
            ElementType aType = ElementType::UNDEFINED ;

            for( Element * tElement : mMesh->elements() )
            {
                if( tElement->is_flagged() )
                {
                    aType = tElement->type() ;
                    break ;
                }
            }

            for( Element * tElement : mMesh->elements() )
            {
                if( tElement->is_flagged() )
                {
                    BELFEM_ERROR( tElement->type() == aType, "inconsistent element type for element %lu",
                                  ( long unsigned int ) tElement->id() );
                }
            }

            return aType ;
        }



//-----------------------------------------------------------------------------

        // for debugging
        void
        CutFactory::save_edges( const uint aIndex )
        {
            OrderedMap< id_t, int > & tMap = mCohomology->get_Generators()( 1 )( aIndex )->getSimplicesMap();

            mMesh->unflag_all_nodes();
            mMesh->unflag_all_edges();

            Cell< Node * > & tOldNodes = mMesh->nodes() ;
            Cell< Edge * > & tOldEdges = mMesh->edges() ;

            index_t tEdgeCount = 0 ;
            index_t tNodeCount = 0 ;

            // count edges
            for ( const auto & pair: tMap )
            {
                tOldEdges( pair.first )->flag() ;
                tOldEdges( pair.first )->flag_nodes() ;
                ++tEdgeCount ;
            }

            for( Node * tNode : tOldNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++tNodeCount ;
                }
            }

            Mesh * tMesh = new Mesh( 3, 0, false);


            Cell< Node * > & tNewNodes = tMesh->nodes();
            tNewNodes.set_size( tNodeCount, nullptr );
            Map< id_t, Node * > tNodeMap ;

            tNodeCount = 0 ;

            for( Node * tOldNode : tOldNodes )
            {
                if( tOldNode->is_flagged() )
                {
                    Node * tNewNode = new Node( tOldNode->id(), tOldNode->x(), tOldNode->y(), tOldNode->z() );
                    tNodeMap[ tOldNode->id() ] = tNewNode ;
                    tNewNodes( tNodeCount++ ) = tNewNode ;
                }
            }

            ElementFactory tFactory ;

            Block * tBlock = new Block( 1, tEdgeCount );

            Cell< Element * > & tElements = tBlock->elements();
            tElements.set_size( tEdgeCount, nullptr );
            tEdgeCount = 0 ;

            for ( const auto & pair: tMap )
            {
                Edge * tEdge = tOldEdges( pair.first );

                Element * tElement = tFactory.create_element( ElementType::LINE2, tEdge->id() );
                tElement->insert_node( tNodeMap( tEdge->node( 0 )->id() ), 0 );
                tElement->insert_node( tNodeMap( tEdge->node( 1 )->id() ), 1 );

                tElements( tEdgeCount++ ) = tElement ;
            }

            tMesh->blocks().push( tBlock );
            tMesh->finalize();

            tMesh->save( "edges.vtk");
            delete tMesh ;
        }


//-----------------------------------------------------------------------------


        void
        CutFactory::collect_nodes_and_elements_on_blocks(
            const Vector< id_t > & aBlockIDs,
            Cell< Element * > & aElements,
            Cell< Node * >    & aNodes )
        {
            index_t tCount = 0 ;

            // count elements
            for ( id_t tID : aBlockIDs )
            {
                Cell< Element * > & tElements =  mMesh->block( tID )->elements();
                for ( Element * tElement: tElements )
                {
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        if ( tElement->node( k )->is_flagged() )
                        {
                            tElement->flag();
                            ++tCount;
                            break ;
                        }
                    }
                }
            }

            // collect elements
            aElements.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( id_t tID : aBlockIDs )
            {
                Cell< Element * > & tElements =  mMesh->block( tID )->elements();
                for ( Element * tElement: tElements )
                {
                    if ( tElement->is_flagged() )
                    {
                        aElements( tCount++ ) = tElement ;
                    }
                }
            }

            // count nodes temporarily
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    ++tCount;
                }
            }

            // collect nodes temporarily
            Cell< Node * > tNodes( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNodes( tCount++ ) = tNode ;
                    tNode->unflag();
                }
            }

            // now we extract the nodes that are connected to our elements
            for ( Element * tElement : aElements )
            {
                tElement->flag_nodes();
            }

            // count nodes that are linked to these elements
            tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }

            // collect nodes that are linked to these elements
            aNodes.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                aNodes( tCount++ ) = tNode ;
            }

            // tidy up flags
            for ( Element * tElement : aElements )
            {
                tElement->unflag_nodes();
            }
            for ( Node * tNode : aNodes )
            {
                tNode->flag();
            }

        }


//-----------------------------------------------------------------------------

        SideSet *
        CutFactory::create_cut_sideset_2d(
                const uint aIndex,
                const ElementType aType,
                Cell< Element * > & aElements,
                const Vector< int > & aCases )
        {

            index_t tNumElements = aElements.size() ;

            id_t tSideSetID = mMesh->max_block_and_sideset_id() ;
            id_t tID = mMesh->max_element_id() ;
            ElementFactory tFactory ;


            // flag all elements in this group
            mMesh->unflag_all_elements() ;
            mMesh->unflag_all_edges();

            for( Element * tElement : aElements )
            {
                tElement->flag();
            }

            for( index_t e=0; e<tNumElements; ++e )
            {
                if( 0 < aCases( e ) )
                {
                    aElements( e )->edge( aCases( e ) - 1 )->flag();
                }
            }

            // count edges
            index_t tCount = 0 ;
            for( Edge * tEdge : mMesh->edges() )
            {
                if( tEdge->is_flagged() )
                {
                    ++tCount ;
                }
            }

            Cell< Edge * > tEdges( tCount, 0 );
            index_t tNumEdges = 0 ;
            for( Edge * tEdge : mMesh->edges() )
            {
                if( tEdge->is_flagged() )
                {
                   tEdge->set_index( tNumEdges );
                   tEdges( tNumEdges++ ) = tEdge ;
                }
            }

            // collect master and slave in temporary containers
            Cell< Element * > tMasters( tNumEdges, nullptr );
            Cell< Element * > tSlaves( tNumEdges, nullptr );
            Vector< uint > tMasterIndices( tNumEdges, gNoIndex );
            Vector< uint > tSlaveIndices( tNumEdges, gNoIndex );

            for( index_t e=0; e<tNumElements; ++e )
            {
                if( 0 < aCases( e ) )
                {
                    uint c = aCases( e ) - 1 ;
                    Element * tElement = aElements( e );
                    Edge * tEdge = tElement->edge( c );
                    index_t k = tEdge->index() ;

                    Element * tOther = nullptr ;
                    uint d = BELFEM_UINT_MAX ;

                    for( uint i=0; i<tElement->number_of_elements(); ++i )
                    {
                        bool tFound = false ;
                        for( uint j=0; j<tElement->element( i )->number_of_edges(); ++j )
                        {
                            if( tElement->element( i )->edge( j )->id() == tEdge->id() )
                            {
                                tOther = tElement->element( i ) ;
                                d = j ;
                            }
                        }

                        if( tFound )
                        {
                            break ;
                        }
                    }

                    tMasters( k ) = tElement ;
                    tMasterIndices( k )  = c ;

                    if( tOther != nullptr )
                    {
                        tSlaves( k ) = tOther ;
                        tSlaveIndices( k ) = d ;
                    }
                }
            }

            tCount = 0;
            for( index_t k=0; k<tNumEdges; ++k )
            {
                Element * tMaster = tMasters( k );
                Element * tSlave  = tSlaves( k );

                if( tSlave != nullptr )
                {
                    if( tMaster->is_flagged() xor tSlave->is_flagged() )
                    {
                        ++tCount ;
                    }
                    else
                    {
                        tEdges( k )->unflag() ;
                    }
                }
                else
                {
                    tEdges( k )->unflag() ;
                }
            }

            // create vertices
            Graph tGraph( tCount, nullptr );
            tCount = 0 ;
            for( Edge * tEdge : tEdges )
            {
                if( tEdge->is_flagged() )
                {
                    graph::Vertex * tVertex = new graph::Vertex();
                    tVertex->set_id( tEdge->id() );
                    tVertex->set_index( tEdge->index() );
                    tEdge->set_index( tCount );
                    tGraph( tCount++ ) = tVertex ;
                }
            }

            // connect vertices
            for( graph::Vertex * tVertex : tGraph )
            {
                Edge * tEdge = tEdges( tVertex->index() );

                uint c = 0 ;
                for( uint k=0; k<2; ++k )
                {
                    Node * tNode = tEdge->node( k ) ;

                    for( uint e=0; e<tNode->number_of_edges(); ++e )
                    {
                        if( tNode->edge( e )->is_flagged() && tNode->edge( e )->id() != tEdge->id() )
                        {
                            ++c ;
                        }
                    }
                }

                Cell< Edge * > tNeighbors( c, nullptr );
                c = 0 ;
                for( uint k=0; k<2; ++k )
                {
                    Node * tNode = tEdge->node( k ) ;

                    for( uint e=0; e<tNode->number_of_edges(); ++e )
                    {
                        if( tNode->edge( e )->is_flagged() && tNode->edge( e )->id() != tEdge->id() )
                        {
                            tNeighbors( c++ ) = tNode->edge( e ) ;
                        }
                    }
                }

                unique( tNeighbors );

                tVertex->init_vertex_container( tNeighbors.size() );

                for( Edge * tEdge2 : tNeighbors )
                {
                    tVertex->insert_vertex( tGraph( tEdge2->index() ) );
                }
            }

            index_t tNumFacets = graph::find_connected_partitions( tGraph );

            // correct Vertex indices
            tCount = 0 ;
            for( Edge * tEdge : tEdges )
            {
                if( tEdge->is_flagged() )
                {
                    tGraph( tEdge->index() )->set_index( tCount );
                }
                ++tCount ;
            }

            mMesh->unflag_all_elements() ;
            for( id_t tID2 : mTopology->non_phi_block_ids() )
            {
                mMesh->block( tID2 )->flag_elements() ;
            }

            tCount = 0 ;

            for( index_t f=0; f<tNumFacets; ++f )
            {
                index_t k = tGraph( f )->index() ;

                Element * tMaster = tMasters(  k  ) ;
                Element * tSlave = tSlaves( k );

                if( ! tMaster->is_flagged() && ! tSlave->is_flagged() )
                {
                    ++tCount ;
                }
            }

            // create sideset
            SideSet * aSideSet = new SideSet( ++tSideSetID, tCount );
            string tFormat = "cut_" + format_with_leading_zeros( mAbstractNodes.size() );
            aSideSet->label() =  sprint( tFormat.c_str(), aIndex + 1) ;

            Cell< Facet * > & tFacets = aSideSet->facets() ;
            tCount = 0 ;

            for( index_t f=0; f<tNumFacets; ++f )
            {
                index_t k = tGraph( f )->index() ;

                Element * tMaster = tMasters(  k  ) ;
                Element * tSlave = tSlaves( k );

                uint m = tMasterIndices( k );
                uint s = tSlaveIndices( k );

                if( ! tMaster->is_flagged() && ! tSlave->is_flagged() )
                {
                    Element * tElement = tFactory.create_element( mesh::element_type_of_facet( tMaster->type(), m ), ++tID );

                    Facet * tFacet = new Facet( tElement );
                    tFacet->set_master( tMaster, m );
                    if( tSlave != nullptr )
                    {
                        tFacet->set_slave( tSlave, s );
                    }

                    tFacets( tCount++ ) = tFacet ;
                }
            }

            for( graph::Vertex * tVertex : tGraph )
            {
                delete tVertex ;
            }

            mMesh->sidesets().push( aSideSet );

            return aSideSet ;
        }

//-----------------------------------------------------------------------------

        SideSet *
        CutFactory::create_cut_sideset_3d(
                const uint aIndex,
                const ElementType aType,
                Cell< Element * > & aElements,
                const Vector< int > & aCases )
        {


            id_t tSideSetID = mMesh->max_block_and_sideset_id() ;
            id_t tID = mMesh->max_element_id() ;


            ElementFactory tFactory ;


            // flag all elements in this group
            mMesh->unflag_all_elements() ;
            for( Element * tElement : aElements )
            {
                tElement->flag();
            }
            // unflag all elements on the air block

            // count faces
            index_t tNumElements = aElements.size();
            index_t tCount = 0 ;
            for( index_t e=0; e<tNumElements; ++e )
            {
                int c = aCases( e ) - 1 ;

                if( 0 <= c && c < 4 )
                {
                    Face * tFace = aElements( e )->face( c );
                    if( tFace->slave() != nullptr )
                    {
                        if( tFace->master()->is_flagged() xor tFace->slave()->is_flagged() )
                        {
                            ++tCount ;
                        }
                    }
                }
            }


            Cell< Face * > tFaces( tCount, nullptr );

            // collect faces and make sure that they are unflagged
            tCount = 0 ;
            for( index_t e=0; e<tNumElements; ++e )
            {
                int c = aCases( e ) - 1 ;

                if( 0 <= c && c < 4 )
                {
                    Face * tFace = aElements( e )->face( c );
                    if( tFace->slave() != nullptr )
                    {
                        if( tFace->master()->is_flagged() xor tFace->slave()->is_flagged() )
                        {
                            tFace->unflag() ;
                            tFaces( tCount++ ) = tFace ;
                        }
                    }
                }
            }

            // here we need to crate a graph and identify the number of disconnected subgraphs
            // unflag all faces
            mMesh->unflag_all_faces() ;

            tCount = 0 ;
            for( Face * tFace : tFaces )
            {
                tFace->flag() ;
                tFace->set_index( tCount++ );
            }


            // collect create graph
            Graph tGraph( tCount, nullptr );
            tCount = 0 ;
            for( Face * tFace : tFaces )
            {

                graph::Vertex * tVertex = new graph::Vertex();
                tVertex->set_id( tFace->id() );
                tVertex->set_index( tCount );

                tGraph( tCount++ ) = tVertex ;
            }

            // connect vertices
            for( graph::Vertex * tVertex : tGraph )
            {
                // get face
                Face * tFace = tFaces( tVertex->index() );

                // temporary unflag to exclude face self-connection
                tFace->unflag() ;

                // count memory
                tCount = 0 ;
                for( uint e=0; e<tFace->number_of_edges(); ++e )
                {
                    Edge * tEdge = tFace->edge( e );
                    for( uint f=0; f<tEdge->number_of_faces(); ++f )
                    {
                        if( tEdge->face( f )->is_flagged() )
                        {
                            ++tCount ;
                        }
                    }
                }

                // connect faces
                Graph tVertices( tCount, nullptr );
                tCount = 0 ;
                for( uint e=0; e<tFace->number_of_edges(); ++e )
                {
                    Edge * tEdge = tFace->edge( e );
                    for( uint f=0; f<tEdge->number_of_faces(); ++f )
                    {
                        if( tEdge->face( f )->is_flagged() )
                        {
                            tVertices( tCount++ ) = tGraph( tEdge->face( f )->index() );
                        }
                    }
                }

                unique( tVertices );
                tVertex->init_vertex_container( tVertices.size() );

                for( graph::Vertex * tNeighbor : tVertices )
                {
                    tVertex->insert_vertex( tNeighbor );
                }

                // set face flag back to true
                tFace->flag();
            }

            graph::dfs( tGraph );
            sort( tGraph, opVertexOwner );

            index_t tNumFacets = graph::find_connected_partitions( tGraph );

            mMesh->unflag_all_elements() ;
            for( id_t tID2 : mTopology->non_phi_block_ids() )
            {
                mMesh->block( tID2 )->flag_elements() ;
            }

            tCount = 0 ;

            for( index_t f=0; f<tNumFacets; ++f )
            {
                // grab face
                Face * tFace = tFaces( tGraph( f )->index() ) ;

                if( ! tFace->master()->is_flagged() && ! tFace->slave()->is_flagged() )
                {
                    ++tCount ;
                }
            }

            // create sideset
            SideSet * aSideSet = new SideSet( ++tSideSetID, tCount );
            string tFormat = "cut_" + format_with_leading_zeros( mAbstractNodes.size() );
            aSideSet->label() =  sprint( tFormat.c_str(), aIndex + 1) ;

            Cell< Facet * > & tFacets = aSideSet->facets() ;

            tCount = 0 ;

            for( index_t f=0; f<tNumFacets; ++f )
            {
                // grab face
                Face * tFace = tFaces( tGraph( f )->index() ) ;

                if( ! tFace->master()->is_flagged() && ! tFace->slave()->is_flagged() )
                {
                    // create element from face
                    Element * tElement = tFactory.create_element(  mesh::element_type_of_facet(aType, tFace->index_on_master() ), ++tID );
                    Facet * tFacet = new Facet( tElement );

                    // link master and slave elements
                    tFacet->set_master( tFace->master(), tFace->index_on_master(), true );
                    if( tFace->slave() != nullptr )
                    {
                        tFacet->set_slave( tFace->slave(), tFace->index_on_slave(), tFace->orientation_on_slave() );
                    }

                    // add element to sideset
                    tFacets( tCount++ ) = tFacet ;
                }
            }

            // tidy up memory
            for( graph::Vertex * tVertex : tGraph )
            {
                delete tVertex ;
            }

            mMesh->sidesets().push( aSideSet );

            return aSideSet ;
        }

//-----------------------------------------------------------------------------


        bool
        CutFactory::create_thin_shell_cuts()
        {
            // grab max id from mesh
            mMaxID = mMesh->max_element_id();

            const Vector< id_t > & tThinShellSideSets = mTopology->groups( DomainType::ThinShell );

            if( tThinShellSideSets.length() == 0 )
            {
                // nothing to do here, return false indicating that we don't have thin shells
                return false ;
            }

            this->create_curves_for_thinshells();

            // in 3D, the terminal curves must be oriented BEFORE the node
            // duplication: the relink steps below swap slave-side facet nodes
            // to duplicates, which breaks the original-node facet search in
            // orient_terminal_curves_sub()
            if ( mMesh->number_of_dimensions() == 3 )
            {
                this->orient_terminal_curves();
            }

            this->duplicate_nodes_on_face_sidesets() ;
            this->relink_slave_elements_with_duplicate_nodes();

            // now let's extend the thin shell sidesets
            this->duplicate_and_relink_facets() ;

            // these are the non-thin shell sidesets
            this->relink_non_thinshell_facets();

            this->close_terminal_loops();

            // add duplicate nodes to mesh
            append_move( mMesh->nodes(), mThinShellDuplicates );

            // ensure proper curve orientation
            if ( mMesh->number_of_dimensions() == 2 )
            {
                this->orient_terminal_curves_2D() ;
            }

            mMesh->update_node_indices() ;

            // return true, indicating that we have thin shells
            return true ;
        }

        void
        CutFactory::duplicate_nodes_on_face_sidesets()
        {
            DynamicBitset tBitset( mMesh->nodes().size() );

            // first, we flag all nodes that sit on tapes
            for( id_t tID : mTopology->groups( DomainType::ThinShell ) )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        tBitset.set( tFacet->node( k )->index() );
                    }
                }
            }

            // next, we must unflag all nodes that belong to side curves
            if (mMesh->number_of_dimensions() == 2)
            {
                //in 2D, we must identify side nodes

                // facet flags may hold scratch state from earlier passes
                // (e.g. fix_facet_masters); the tip test below requires that
                // only tape facets are flagged
                mMesh->unflag_all_facets() ;

                //First flag all the edges on thin shells
                for ( Protoshell * tProtoshell : mProtoshells )
                {
                    for ( id_t tID : tProtoshell->sidesets() )
                    {
                        mMesh->sideset( tID )->flag_all_facets() ;
                    }
                }

                //then identify the nodes with only one edge flagged
                for ( Protoshell * tProtoshell : mProtoshells )
                {
                    for ( Curve * tCurve : tProtoshell->terminal_curves() )
                    {
                        uint tNumTips = 0 ;

                        for (Node * tNode : tCurve->nodes() )
                        {
                            uint tCount = 0 ;
                            for (uint i = 0 ; i < tNode->number_of_facets(); ++i)
                            {
                                if (tNode->facet(i)->is_flagged())
                                {
                                    tCount++ ;
                                }
                            }
                            if (tCount == 1)
                            {
                                tBitset.reset( tNode->index() );
                                ++tNumTips ;
                            }
                        }

                        // an open tape must end in exactly two tips; anything else
                        // means the endpoints would be duplicated and the tape
                        // silently torn open
                        BELFEM_ERROR( tNumTips == 2 || ( tCurve->is_closed() && tNumTips == 0 ),
                                      "Terminal curve %lu found %u tip nodes on the tape (expect 2 for an open tape, 0 for a closed loop)",
                                      ( long unsigned int ) tCurve->id(),
                                      ( unsigned int ) tNumTips );
                    }
                }
            }
            else
            {
                //Identify the open and closed side curves
                for ( Protoshell * tProtoshell : mProtoshells )
                {
                    //Flag the nodes on the side curve of this shell
                    for ( Curve * tCurve : tProtoshell->side_curves() )
                    {
                        bool tOpen = false ;
                        for ( Node* tNode : tCurve->nodes() )
                        {
                            tNode->flag();
                        }

                        //Now look if other thin shells have side nodes flagged, if so, this curve must be open
                        for (  Protoshell * tOtherShell : mProtoshells )
                        {
                            if (tProtoshell->id() == tOtherShell->id()) continue ;
                            for ( Curve * tCurve2 : tOtherShell->side_curves() )
                            {
                                for ( Node* tNode : tCurve2->nodes() )
                                {
                                    if (tNode->is_flagged())
                                    {
                                        tOpen = true;
                                        break;
                                    }
                                }
                                if (tOpen) break ;
                            }
                            if (tOpen) break ;
                        }
                        tCurve->set_closed_flag( !tOpen ) ;

                        //Unflag node
                        for ( Node* tNode : tCurve->nodes() )
                        {
                            tNode->unflag();
                        }

                    }

                }

                for ( Protoshell * tProtoshell : mProtoshells )
                {
                    for ( Curve * tCurve : tProtoshell->side_curves() )
                    {
                        //Only applied for closed side curves
                        if (tCurve->is_closed())
                        {
                            for ( Node * tNode : tCurve->nodes() )
                            {
                                tBitset.reset( tNode->index() );
                            }
                        }
                    }
                }
            }

            // grab the nodes that are to be duplicated
            Cell< index_t > tIndices ;
            tBitset.where( tIndices );

            index_t tCount = 0 ;

            mThinShellDuplicates.set_size( tIndices.size(), nullptr );


            // we select the nodes again, this tome for the originals
            for( id_t tID : mTopology->groups( DomainType::ThinShell ) )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        tBitset.set( tFacet->node( k )->index() );
                    }
                }
            }

            // first, we reset all node indices,
            // this will force an error if we do something wrong
            Cell< Node * > & tNodes = mMesh->nodes() ;
            for ( Node * tNode : tNodes )
            {
                tNode->set_index( gNoIndex );
            }

            id_t tID = mMesh->max_node_id() ;

            for ( index_t tIndex : tIndices )
            {
                Node * tOriginal = tNodes( tIndex ) ;
                tOriginal->set_index( tCount );

                Node * tDuplicate = new Node( ++tID, tOriginal->x(), tOriginal->y(), tOriginal->z() );


                //#EYE eye hack for visualization
                /*
                real x = tOriginal->x() ;
                real y = 0.05*(1-x*x);
                real z = tOriginal->z() ;
                tOriginal->set_coords( x, -y, z );
                tDuplicate->set_coords( x, y, z ); */

                mThinShellDuplicates( tCount++ ) = tDuplicate ;
            }

            // with the duplicates created, we must now fix the periodicities
            if ( mMesh->has_periodicity() )
            {
                Periodicity * tPeriodicity = mMesh->periodicity();

                tCount = 0 ;
                for ( index_t tIndex : tIndices )
                {
                    Node * tA = tNodes( tIndex ) ;
                    Node * tC = mThinShellDuplicates( tCount++ ) ;

                    if ( tA->is_periodic() && ! tC->is_flagged() )
                    {
                        Node * tB = tA->periodic();
                        Node * tD = mThinShellDuplicates( tB->index() ) ;

                        tC->set_periodic( tD );
                        tD->set_periodic( tC );
                        tC->flag();
                        tD->flag();

                        tPeriodicity->add_node_pair_to_backup( tC, tD );

                        // flags for master and slave sets
                        if ( tA->is_flagged( 1 ) ) tC->flag( 1 );
                        if ( tA->is_flagged( 2 ) ) tC->flag( 2 );
                        if ( tB->is_flagged( 1 ) ) tD->flag( 1 );
                        if ( tB->is_flagged( 2 ) ) tD->flag( 2 );
                    }
                }
            }

            tBitset.where( tIndices );

            mThinShellMasterNodes.set_size( tIndices.size(), nullptr );
            mThinShellSlaveNodes.set_size( tIndices.size(), nullptr );

            tCount = 0 ;

            for ( index_t tIndex : tIndices )
            {
                Node * tNode = tNodes( tIndex ) ;

                mThinShellMasterNodes( tCount ) = tNode ;

                // check if duplicate exists
                if ( tNode->index() != gNoIndex )
                {
                    mThinShellSlaveNodes( tCount++ ) = mThinShellDuplicates( tNode->index() ) ;
                }
                else
                {
                    mThinShellSlaveNodes( tCount++ ) = tNode ;
                }
            }


        }

        void
        CutFactory::relink_slave_elements_with_duplicate_nodes()
        {
            mMesh->unflag_all_nodes();
            Cell< Element * > & tElements = mMesh->elements() ;

            for ( Element * tElement : tElements )
            {
                tElement->unflag( 0 );
                tElement->unflag( 1 );
            }
            Cell< Block * > & tBlocks = mMesh->blocks() ;

            index_t tCount = 0 ;
            for ( Block * tBlock : tBlocks )
            {
                tBlock->set_index( tCount++ );
            }

            DynamicBitset tElementBitset( tElements.size() );
            DynamicBitset tBlockBitset( tBlocks.size() );

            Cell< index_t > tIndices ;
            for( id_t tID : mTopology->groups( DomainType::ThinShell ) )
            {
                tElementBitset.reset();
                tBlockBitset.reset();
                SideSet * tSideSet = mMesh->sideset( tID ) ;


                Cell<Facet * > & tFacets = tSideSet->facets() ;

                for( Facet * tFacet : tFacets )
                {
                    tBlockBitset.set( mMesh->block( tFacet->slave()->block_id() )->index() );
                }


                tBlockBitset.where( tIndices );
                for ( index_t b : tIndices )
                {
                    tBlocks( b )->flag_elements();
                }

                Cell< Node * > & tNodes = tSideSet->nodes() ;
                for( Node * tNode : tNodes )
                {
                    for ( uint e=0; e<tNode->number_of_elements(); ++e )
                    {
                        Element * tElement = tNode->element( e ) ;
                        if ( tElement->is_flagged() )
                        {
                            tElementBitset.set( tElement->index() );
                        }
                    }
                }
                for ( index_t b : tIndices )
                {
                    tBlocks( b )->unflag_elements();
                }

                tSideSet->flag_all_nodes();
                tElementBitset.where( tIndices );

                for ( index_t e : tIndices )
                {
                    Element * tElement = tElements( e ) ;
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        Node * tNode = tElement->node( k ) ;
                        if ( tNode->is_flagged() && tNode->index() != gNoIndex )
                        {
                            // tag that this element has been altered
                            tElement->flag( 1 );
                            tElement->insert_node( mThinShellDuplicates( tNode->index() ), k );
                        }
                    }
                }
                tSideSet->unflag_all_nodes();
            }
        }


        void
        CutFactory::duplicate_and_relink_facets()
        {
            mMesh->unflag_all_facets( 1 );

            for( id_t tID : mTopology->groups( DomainType::ThinShell ) )
            {
                // first we create a temporary sideset
                SideSet * tSideSet = new SideSet( tID, mMesh->sideset( tID )->number_of_facets() );

                // remember sideset for later
                mTemporaryThinShellSidesets.push( tSideSet );

                // move the original facets
                tSideSet->facets().vector_data() = std::move( mMesh->sideset( tID )->facets().vector_data() );

                // get the originals
                Cell< Facet * > & tOriginals = tSideSet->facets() ;
                Cell< Facet * > & tDuplicates = mMesh->sideset( tID )->facets() ;

                index_t tNumFacets = tOriginals.size() ;
                tDuplicates.set_size( 2 * tNumFacets, nullptr );

                index_t tCountA = 0 ;
                index_t tCountB = tNumFacets ;
                id_t tFacetIdA = mMaxID ;
                id_t tFacetIdB = mMaxID + tNumFacets ;

                ElementFactory tFactory ;

                for ( Facet * tOrg : tOriginals )
                {
                    // create a new facet
                    Facet * tA = new Facet( tFactory.create_element( tOrg->element()->type(), ++tFacetIdA ) );
                    Facet * tB = new Facet( tFactory.create_element( tOrg->element()->type(), ++tFacetIdB ) );

                    tA->set_master( tOrg->master(), tOrg->index_on_master() );
                    tB->set_master( tOrg->slave(), tOrg->index_on_slave() );

                    // tag needed for facet relinking
                    tOrg->master()->flag( 1 );
                    tOrg->slave()->flag( 1 );

                    // add duplicates to container
                    tDuplicates( tCountA++ ) = tA ;
                    tDuplicates( tCountB++ ) = tB ;

                    tA->flag( 1 );
                    tB->flag( 1 );
                }

                mMaxID = tFacetIdB ;
            }

            mMesh->update_facet_indices();
        }

        void
        CutFactory::relink_non_thinshell_facets()
        {
            Cell< Facet * > & tFacets = mMesh->facets() ;

            Cell< Node * > tNodes ;
            for ( Facet * tFacet : tFacets )
            {
                // facets without master shouldn't exist,
                // but we check anyways
                if ( ! tFacet->has_master() ) continue ;

                Element * tMaster = tFacet->master() ;

                if ( ! tMaster->is_flagged( 1 ) ) continue;

                // we can skip facets that have been processed by
                // duplicate_and_relink_facets
                if ( tFacet->is_flagged( 1 ) )
                {
                    tFacet->unflag( 1 );
                    continue;
                }

                tMaster->get_nodes_of_facet( tFacet->index_on_master(), tNodes );

                // relink the facet
                uint k = 0 ;
                Element * tElement = tFacet->element();
                for ( Node * tNode : tNodes )
                {
                    tElement->insert_node( tNode, k++ );
                }
            }

            mMesh->unflag_all_elements( 1 );
        }

        void
        CutFactory::close_terminal_loops()
        {
            ElementFactory tFactory ;

            for ( Protoshell * tProtoshell : mProtoshells )
            {
                for ( Curve * tCurve : tProtoshell->terminal_curves() )
                {
                    // container for originals
                    Cell< Segment * > & tSegments = tCurve->segments() ;

                    if ( tSegments.size() == 0 ) continue ;

                    // container for duplicate segments
                    Cell< Segment * > tDuplicates( tSegments.size(), nullptr );

                    // get first element type
                    ElementType tType = tSegments( 0 )->element()->type() ;

                    index_t tCount = tSegments.size() ;

                    // create new segments
                    for ( Segment * tOrg : tSegments )
                    {
                        // create a new segment
                        Segment * tDup = new Segment( tFactory.create_element( tType, ++mMaxID ) );

                        // collect nodes
                        uint n = tOrg->number_of_nodes();
                        for ( uint k=0; k<n; ++k )
                        {
                            Node * tNode = tOrg->node( k )->index() == gNoIndex ? tOrg->node( k ) : mThinShellDuplicates( tOrg->node( k )->index() );
                            tDup->element()->insert_node( tNode, n-1-k );
                        }

                        // add duplicate to temporary container
                        tDuplicates( --tCount ) = tDup ;
                    }

                    // fix the arc length coordinates
                    switch ( tType )
                    {
                        case( ElementType::LINE2 ) :
                        {
                            // Get initial arclengths and segment count
                            index_t n = tCurve->segments().size();
                            Vector<real >& tS = tCurve->arclength();
                            Vector<real> tS0( tS ); // Copy initial arclengths

                            // Resize arclength array
                            tS.set_size(2 * n + 1);

                            // Copy initial arclengths
                            index_t c = 0;
                            for (real s : tS0)
                            {
                                tS(c++) = s;
                            }

                            // Append cumulative arclengths
                            real s = tS(c - 1);
                            for (index_t k = 0; k < n; ++k)
                            {
                                s += tS0(k + 1) - tS0(k);
                                tS(c++) = s;
                            }

                            break;
                        }
                        case( ElementType::LINE3 ) :
                        {
                            // Get initial arclengths and segment count
                            index_t n = tCurve->segments().size();
                            Vector<real>& tS = tCurve->arclength();
                            Vector<real> tS0( tS ); // Copy initial arclengths

                            // Resize arclength array
                            tS.set_size(4 * n + 1);

                            // Copy initial arclengths
                            index_t c = 0;
                            for (real s : tS0) {
                                tS(c++) = s;
                            }

                            // Append cumulative arclengths with midpoints
                            real s = tS(c - 1);
                            for (index_t k = 0; k < n; ++k)
                            {
                                real l = tS0(2 * k + 2) - tS0(2 * k); // Segment length
                                tS(c++) = s + 0.5 * l; // Midpoint
                                s += l;
                                tS(c++) = s; // Endpoint
                            }
                            break;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Invalid element type at terminal loop %lu",
                                          ( long unsigned int ) tCurve->id() );
                        }
                    }

                    // append duplicates to segmnent list
                    append_move( tSegments, tDuplicates );

                    // fix the node container
                    Cell< Node * > & tNodes = tCurve->nodes() ;

                    switch ( tType )
                    {
                        case ElementType::LINE2 :
                        {
                            // we only copy the first node because this is a closed loop
                            tNodes.set_size( tSegments.size(), nullptr );
                            tCount = 0 ;
                            for ( Segment * tSegment : tSegments )
                            {
                                tNodes( tCount++ ) = tSegment->node( 0 );
                            }
                            break;
                        }
                        case ElementType::LINE3 :
                        {
                            // we only copy the first node because this is a closed loop
                            tNodes.set_size( 2 * tSegments.size(), nullptr );
                            tCount = 0 ;
                            for ( Segment * tSegment : tSegments )
                            {
                                tNodes( tCount++ ) = tSegment->node( 0 );
                                tNodes( tCount++ ) = tSegment->node( 2 ); // node 2 is the center node
                            }
                            break;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Invalid element type at terminal loop %lu",
                                          ( long unsigned int ) tCurve->id() );
                        }
                    }

                    tCurve->set_closed_flag( true );
                }
            }
        }

        void
        CutFactory::orient_terminal_curves()
        {
            // To identify the orientation of a curve, we must compare the
            // tangential vector of each edge with the normal vectors of the boundary
            // and the thin shell sideset. In practice, however, it is sufficient
            // to just check the first segment of each curve
            mMesh->unflag_all_nodes();

            Vector< real > tP( 3 );
            Vector< real > tQ( 3 );

            Vector< real > tT( 3 );
            Vector< real > tS( 3 );
            Vector< real > tB( 3 );
            Vector< real > tR( 3 );

            for ( Protoshell * tProtoshell : mProtoshells )
            {
                Cell< Curve * > & tCurves = tProtoshell->terminal_curves() ;
                for ( Curve * tCurve : tCurves )
                {
                    // grab the first segment of the curve
                    Segment * tSegment = tCurve->segments()( 0 );

                    // compute the tangential vector
                    tT( 0 ) = tSegment->node( 1 )->x() - tSegment->node( 0 )->x();
                    tT( 1 ) = tSegment->node( 1 )->y() - tSegment->node( 0 )->y();
                    tT( 2 ) = tSegment->node( 1 )->z() - tSegment->node( 0 )->z();
                    tT /= norm( tT );

                    tSegment->node( 0 )->flag();
                    tSegment->node( 1 )->flag();

                    // next, we look for the surface that sits on the shell and compute its normal
                    this->orient_terminal_curves_sub( tCurve->sideset_a(), tP, tQ, tB );

                    // this is the surface for the boundary
                    this->orient_terminal_curves_sub( tCurve->sideset_b(), tP, tQ, tS );

                    // tidy up
                    tSegment->node( 0 )->unflag();
                    tSegment->node( 1 )->unflag();

                    // this is the direction vector
                    tR = cross( tS, tB );
                    tR /= norm( tR );

                    // if the difference between T and R is close 0, the curve runs counter clockwise, if it is close to 2,
                    // it runs clockwise and must be flipped. We chose 1.0 as criterion to allow for imprecisions
                    // for higher order elements
                    if ( norm( tT - tR ) > 1.0 )
                    {
                        tCurve->reverse();
                    }

                }
            }
        }

        void
        CutFactory::orient_terminal_curves_sub(
            SideSet * aSideSet,
            Vector< real > & aWorkA,
            Vector< real > & aWorkB,
            Vector< real > & aN )
        {
            Facet * tFacet = nullptr ;
            for ( Facet * tF : aSideSet->facets() )
            {
                // count flagged nodes
                uint tCount = 0 ;
                for ( uint k=0; k<tF->number_of_corner_nodes(); ++k )
                {
                    if ( tF->node( k )->is_flagged() )
                    {
                        ++tCount ;
                    }
                }

                if ( tCount > 1 )
                {
                    // we have a surface
                    tFacet = tF ;
                    break ;
                }
            }

            BELFEM_ASSERT( tFacet != nullptr, "No surface found on tape/boundary" );

            // next, we compute the normal of the tape
            aWorkA( 0 ) = tFacet->node( 1 )->x() - tFacet->node( 0 )->x();
            aWorkA( 1 ) = tFacet->node( 1 )->y() - tFacet->node( 0 )->y();
            aWorkA( 2 ) = tFacet->node( 1 )->z() - tFacet->node( 0 )->z();
            aWorkB( 0 ) = tFacet->node( 2 )->x() - tFacet->node( 0 )->x();
            aWorkB( 1 ) = tFacet->node( 2 )->y() - tFacet->node( 0 )->y();
            aWorkB( 2 ) = tFacet->node( 2 )->z() - tFacet->node( 0 )->z();
            aN = cross( aWorkA, aWorkB );
            aN /= norm( aN );
        }

        void
        CutFactory::orient_terminal_curves_2D()
        {

            mMesh->unflag_all_nodes();

            Vector< real > tP( 3 );

            Vector< real > tT( 3 );
            Vector< real > tS( 3 );
            Vector< real > tB( 3 );
            Vector< real > tR( 3 );

            //in 2-D the surface is always in +z
            tS(2) = 1.0 ;

            for ( Protoshell * tProtoshell : mProtoshells )
            {
                Cell< Curve * > & tCurves = tProtoshell->terminal_curves() ;
                for ( Curve * tCurve : tCurves )
                {
                    // grab the first segment of the curve
                    Segment * tSegment = tCurve->segments()( 0 );

                    // compute the tangential vector
                    tT( 0 ) = tSegment->node( 1 )->x() - tSegment->node( 0 )->x();
                    tT( 1 ) = tSegment->node( 1 )->y() - tSegment->node( 0 )->y();
                    tT( 2 ) = tSegment->node( 1 )->z() - tSegment->node( 0 )->z();
                    tT /= norm( tT );

                    tSegment->node( 0 )->flag();
                    tSegment->node( 1 )->flag();

                    Facet * tFacet = nullptr ;

                    for ( Facet * tF : tCurve->sideset_a()->facets() )
                    {
                        // count flagged nodes
                        uint tCount = 0 ;
                        for ( uint k=0; k<tF->number_of_corner_nodes(); ++k )
                        {
                            if ( tF->node( k )->is_flagged() )
                            {
                                ++tCount ;
                            }
                        }

                        if ( tCount > 1 )
                        {
                            // we have a surface
                            tFacet = tF ;
                            break ;
                        }
                    }

                    BELFEM_ASSERT( tFacet != nullptr, "No surface found on tape/boundary" );

                    // next, we compute the normal of the tape
                    tP( 0 ) = tFacet->node( 1 )->x() - tFacet->node( 0 )->x();
                    tP( 1 ) = tFacet->node( 1 )->y() - tFacet->node( 0 )->y();
                    tB(0) = -tP(1);
                    tB(1) = tP(0) ;
                    tB /= norm( tB );

                    // tidy up
                    tSegment->node( 0 )->unflag();
                    tSegment->node( 1 )->unflag();

                    // this is the direction vector
                    tR = cross( tS, tB );
                    tR /= norm( tR );

                    // if the difference between T and R is close 0, the curve runs counter clockwise, if it is close to 2,
                    // it runs clockwise and must be flipped. We chose 1.0 as criterion to allow for imprecisions
                    // for higher order elements
                    if ( norm( tT - tR ) > 1.0 )
                    {
                        tCurve->reverse();
                    }
                }
            }
        }

        void
        CutFactory::create_curves_for_thinshells()
        {
            if ( mMesh->number_of_dimensions()  == 3 )
            {
                this->create_side_curves_for_thinshells_3d() ;
            }
            mMesh->unflag_all_nodes();
            mMesh->create_curve_map();
        }


        void
        CutFactory::collect_boundary_sidesets( Vector< id_t > & aIDs )
        {
            // count sidesets
            index_t tCount = 0 ;
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {

                if ( tSideSet->number_of_facets() > 0 )
                {
                    // within one sideset, we know that either all facets have slaves or none.
                    if ( ! tSideSet->facets()(0)->has_slave() )
                    {
                        ++tCount ;
                    }
                }
            }

            // allocate memory
            aIDs.set_size( tCount );
            tCount = 0 ;

            // collect sidesets
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {

                if ( tSideSet->number_of_facets() > 0 )
                {
                    // within one sideset, we know that either all facets have slaves or none.
                    if ( ! tSideSet->facets()(0)->has_slave() )
                    {
                        aIDs( tCount++ ) = tSideSet->id() ;
                    }
                }
            }
        }


        void
        CutFactory::create_side_curves_for_thinshells_3d()
        {
            CurveFactory tFactory( mMesh );

            Vector< id_t > tBoundaries ;
            this->collect_boundary_sidesets( tBoundaries );

            index_t tCount = 0 ;

            for ( Protoshell * tProtoshell : mProtoshells )
            {
                // create the side curves for the sidesets
                 tProtoshell->side_curves() = tFactory.thin_shell_side_curves( tProtoshell->sidesets(), tBoundaries ) ;
            }

            // set labels for side curves
            tCount = 0 ;

            for ( Protoshell * tProtoshell : mProtoshells )
            {
                for ( Curve * tCurve : tProtoshell->side_curves() )
                {
                    tCurve->label() = sprint( "side_%u", ( unsigned int ) ++tCount );
                }
            }

            //Add side curves to the mesh
            Cell< Curve * > & tCurves = mMesh->curves() ;
            for ( Protoshell * tProtoshell : mProtoshells )
            {
                for ( Curve * tCurve : tProtoshell->side_curves() )
                {
                    tCurves.push(tCurve) ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::flag_nodes_and_facets_of_tape_sidesets( Cell< Node * > & aNodes )
        {
            mMesh->unflag_all_nodes() ;
            mMesh->unflag_all_facets() ;
            mMesh->unflag_all_elements() ;

            // first we flag all nodes on the thin shell sidesets
            for( id_t tID : mTopology->groups( DomainType::ThinShell ) )
            {
                mMesh->sideset( tID )->flag_all_nodes() ;
                mMesh->sideset( tID )->flag_all_facets() ;
            }

            Cell< Node * > & tNodes = mMesh->nodes();

            index_t aCount = 0 ;

            if ( mMesh->number_of_dimensions() == 2 )
            {
                // reindex nodes

                for( Node * tNode : tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        tNode->set_index( aCount++ );
                    }
                }

                // count facets per node
                Vector< uint > tNumFacets( aCount, 0 );
                for( Facet * tFacet : mMesh->facets() )
                {
                    if( tFacet->is_flagged() )
                    {
                        for( uint k=0; k<tFacet->number_of_nodes(); ++k )
                        {
                            ++tNumFacets( tFacet->node( k )->index() );
                        }
                    }
                }

                aCount = 0 ;

                for( Node * tNode : mMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        if( tNumFacets( tNode->index() ) > 1 )
                        {
                            ++aCount ;
                        }
                        else
                        {
                            tNode->set_index( gNoIndex );
                            tNode->unflag();
                        }
                    }
                }
            }
            else
            {
                // remove nodes on tape edges
                for ( Protoshell * tProtoshell : mProtoshells )
                {
                    for ( Curve * tCurve : tProtoshell->side_curves() )
                    {
                        for ( Node * tNode : tCurve->nodes() )
                        {
                            tNode->unflag();
                        }
                    }
                }

                // count flagged nodes
                for( Node * tNode : tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        tNode->set_index( aCount++ );
                    }
                    else
                    {
                        tNode->set_index( gNoIndex );
                    }
                }
            }

            aNodes.set_size( aCount, nullptr );
            aCount = 0 ;
            for( Node * tNode : tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    aNodes( aCount++ ) = tNode ;
                }
            }
        }

//-----------------------------------------------------------------------------

            Facet *
            CutFactory::create_facet( Element * aElement, Cell< Element * > & aCandidates )
            {
                Cell< Node * > tNodes ;

                // reset nodes of this child
                for ( Element * tCandidate: aCandidates )
                {
                    tCandidate->unflag_nodes();
                }

                aElement->flag_corner_nodes() ;

                for ( Element * tCandidate: aCandidates )
                {
                    uint tNumFacets = tCandidate->number_of_facets() ;

                    for( uint f=0; f<tNumFacets; ++f )
                    {
                        tCandidate->get_corner_nodes_of_facet( f, tNodes );

                        uint tCount = 0 ;

                        for( Node * tNode : tNodes )
                        {
                            if( tNode->is_flagged() )
                            {
                                ++tCount ;
                            }
                        }

                        if( tCount == aElement->number_of_corner_nodes() )
                        {
                            Facet * aFacet = new Facet( aElement );
                            aFacet->set_master( tCandidate, f );

                            for ( Element * tElement: aCandidates )
                            {
                                tElement->unflag_nodes();
                            }

                            return aFacet ;
                        }
                    }
                }

                BELFEM_ERROR( false, "something went wrong while creating cut facets");
                return nullptr ;
            }

//-----------------------------------------------------------------------------

        bool
        CutFactory::connect_facet_to_slave( Facet * aFacet )
        {
            // get number of nodes from this facet
            uint tNumNodes = aFacet->element()->number_of_corner_nodes();

            // count max size of candidates
            uint tCount = 0;
            for ( uint k = 0; k < tNumNodes; ++k )
            {
                mesh::Node * tNode = aFacet->node( k );

                tCount += tNode->number_of_elements();
            }

            // create list of Element candidates
            Cell<mesh::Element *> tCandidates( tCount, nullptr );

            // reset counter
            tCount = 0;

            // populate candidates
            for ( uint k = 0; k < tNumNodes; ++k )
            {
                mesh::Node * tNode = aFacet->node( k );

                for ( uint e = 0; e < tNode->number_of_elements(); ++e )
                {
                    tCandidates( tCount++ ) = tNode->element( e );
                }
            }

            // make result unique
            unique( tCandidates );

            Cell<mesh::Node *> tNodes;

            mesh::Element * tMaster = aFacet->master() ;

            for ( mesh::Element * tElement : tCandidates )
            {
                uint tIndex = compute_facet_index( aFacet, tElement, tNodes );

                if ( tIndex < BELFEM_UINT_MAX )
                {
                    if ( tElement->id() != tMaster->id() )
                    {;
                        aFacet->set_slave( tElement, tIndex );
                        return true ;
                    }
                }
            } // end loop over candidates

            return false ;
        }

//-----------------------------------------------------------------------------

        Cell< Node * > &
        CutFactory::abstract_nodes()
        {
            return mAbstractNodes ;
        }

//-----------------------------------------------------------------------------

        Cell< Node * > &
        CutFactory::orphaned_nodes()
        {
            return mOrphanedNodes ;
        }

//-----------------------------------------------------------------------------

        Cell< SideSet * > &
        CutFactory::cuts()
        {
            return mCuts ;
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::link_node_duplicates_and_originals()
        {
            mMesh->unflag_all_nodes() ;
            Cell< Node * > & tNodes = mMesh->nodes();

            // count hanging nodes
            index_t tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_hanging() )
                {
                    ++tCount ;
                    for ( uint s=0; s<tNode->number_of_sources(); ++s )
                    {
                        // check if this source is a node
                        if ( tNode->source( s )->entity_type() == EntityType::NODE )
                        {
                            tNode->source( s )->flag();
                        }
                    }
                }
            }

            // unflag the abstract nodes
            for ( Node * tNode : mAbstractNodes )
            {
                tNode->unflag();
            }

            // collect hanging nodes
            Cell< Node * > tDuplicates( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_hanging() )
                {
                    tDuplicates( tCount++ ) = tNode ;
                }
            }

            // count source nodes
            tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }

            // collect source nodes
            Cell< Node * > tOriginals( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    tOriginals( tCount++ ) = tNode ;
                }
            }

            // count duplicates per node
            Vector< uint > tNumDuplicates( tCount, 0 );
            for ( Node * tNode : tDuplicates )
            {
                for ( uint s=0; s<tNode->number_of_sources(); ++s )
                {
                    // check if this source is a node
                    if ( tNode->source( s )->entity_type() == EntityType::NODE && tNode->source( s )->is_flagged() )
                    {
                        ++tNumDuplicates( tNode->source( s )->index() );
                    }
                }
            }


            // allocate duplicate containers
            for ( Node * tNode : tOriginals )
            {
                tNode->allocate_duplicate_container( tNumDuplicates( tNode->index() ) );
            }

            for ( Node * tDup : tDuplicates )
            {
                for ( uint s=0; s<tDup->number_of_sources(); ++s )
                {
                    // check if this source is a node
                    if ( tDup->source( s )->entity_type() == EntityType::NODE && tDup->source( s )->is_flagged() )
                    {
                        // grab node
                        Node * tOrg = reinterpret_cast< Node * >( tDup->source( s ) ) ;

                        // link node with original
                        tOrg->original()->add_duplicate( tDup );
                        tDup->set_original( tOrg->original() );
                    }
                }
            }
        }

        void
        CutFactory::collect_orphan_nodes()
        {
            mMesh->unflag_all_nodes() ;
            for ( Node * tNode : mAbstractNodes )
            {
                tNode->flag();
            }
            for ( Block * tBlock : mMesh->blocks() )
            {
                tBlock->flag_nodes() ;
            }

            index_t tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( ! tNode->is_flagged() )
                {
                    ++tCount ;
                }
            }
            if ( tCount == 0 ) return ;

            mOrphanedNodes.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( ! tNode->is_flagged() )
                {
                    mOrphanedNodes( tCount++ ) = tNode ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::unflag_symmetry_sidesets()
        {
            if (mAlgorithm == CutAlgorithm::PellikkaGeneralized ||
                mAlgorithm == CutAlgorithm::Pellikka ||
                mAlgorithm == CutAlgorithm::CCR)
            {
                for(SideSet * tSideset : mMesh->sidesets())
                {
                    DomainType tType = tSideset->domain_type();
                    if ( tType == DomainType::AirSymmetry ||
                        tType == DomainType::Symmetry ||
                        tType == DomainType::FerroSymmetry ||
                        tType == DomainType::ConductorSymmetry )
                    {
                        for(Facet * tFacet : tSideset->facets())
                        {
                            for (uint i = 0 ; i < tFacet->number_of_edges(); ++i)
                            {
                                //Unflagging the simplices along the symmetry plane just makes sure that they don't interact with the cut
                                //This is just forcing the first reduction algorithm steps to be on the cut
                                tFacet->edge(i)->unflag();
                                for(uint j = 0 ; j < tFacet->edge(i)->number_of_nodes(); ++j)
                                {
                                    tFacet->edge(i)->node(j)->unflag() ;
                                }

                            }
                            for (uint i = 0 ; i < tFacet->number_of_faces(); ++i)
                            {
                                tFacet->face(i)->unflag();
                            }
                        }
                        for (Face * tFace : mMesh->faces())
                        {
                            if (tFace->is_flagged())
                            {
                                uint tCountEdges = 0;
                                for (uint i = 0; i < tFace->number_of_edges();++i)
                                {
                                    if (!tFace->edge(i)->is_flagged())
                                    {
                                        tCountEdges++;
                                    }
                                }
                                if (tCountEdges == 3)
                                {
                                    tFace->unflag();
                                }
                            }
                        }
                    }

                }
            }

        }

//-----------------------------------------------------------------------------

        void
        CutFactory::save_curve_debug_meshes()
        {
            for ( Curve * tCurve : mMesh->curves() )
            {
                std::cout << "#SAVING CURVE " << tCurve->id() << " " << tCurve->is_terminal() << " " << tCurve->is_closed() << std::endl ;
                string tLabel = sprint( "curve_%u.vtk", tCurve->id() ) ;
                tCurve->save( tLabel.c_str() );
            }
        }

        void
        CutFactory::compute_element_adjacencies()
        {
            mMesh->update_element_indices();

            if ( mMesh->number_of_dimensions() == 2 )
            {
                // there are no faces in 2D, so we cannot rebuild the
                // adjacencies here. The cut surgery has invalidated any
                // cached element-to-element data, so we drop it and let
                // the next kernel recompute it from scratch.
                for ( Element * tElement : mMesh->elements() )
                {
                    tElement->reset_element_container();
                }
                mMesh->reset_connectivity( Connectivity::ElementToElement );
                return ;
            }

            Cell< Face * > & tFaces = mMesh->faces();

            Vector< index_t > tCount( mMesh->number_of_elements(), 0 );

            for ( Face * tFace : tFaces )
            {
                if ( tFace->master() == nullptr || tFace->slave() == nullptr ) continue ;
                ++tCount( tFace->master()->index() );
                ++tCount( tFace->slave()->index() );
            }

            Cell< Element * > & tElements = mMesh->elements();
            for ( Element * tElement : tElements )
            {
                tElement->allocate_element_container( tCount( tElement->index() ));
            }
            for ( Face * tFace : tFaces )
            {
                if ( tFace->master() == nullptr || tFace->slave() == nullptr ) continue ;
                tFace->master()->insert_element( tFace->slave() );
                tFace->slave()->insert_element( tFace->master() );
            }

            // caution: not valid for conducting elements and thin shells
            // we need these data for the thin shell factory
            mMesh->set_connectivity( Connectivity::ElementToElement );
        }

//-----------------------------------------------------------------------------

        void
        CutFactory::compute_poisson_problem()
        {

                fem::KernelParameters * tParams = new fem::KernelParameters( mMesh );
                fem::Kernel tKernel( tParams );

                tKernel.create_field(  tKernel.create_equation( IwgType::Poisson  ) );
                fem::DofManager * tField = tKernel.dofmgr();

                SolverParameters tSolverParams( gDefaultSolver );
                tField->set_solver(  tSolverParams );
                if ( mCommRank == 0 )
                {

                    mMesh->unflag_all_nodes( 0 );
                    mMesh->unflag_all_nodes( 1 );

                    for ( Block * tBlock : mMesh->blocks() )
                    {
                        if ( tBlock->domain_type() == DomainType::Conductor )
                        {
                            for ( Element * tElement : tBlock->elements() )
                            {
                                tElement->flag_nodes( 1 );
                            }
                        }
                    }



                    for ( SideSet * tSideSet : mMesh->sidesets() )
                    {
                        if ( tSideSet->domain_type() == DomainType::ThinShell || tSideSet->domain_type() == DomainType::AirPeriodic )
                        {
                            for ( Facet * tFacet : tSideSet->facets() )
                            {
                                tFacet->flag_nodes( 1 );
                            }
                        }
                        else if ( tSideSet->domain_type() == DomainType::AirAntiSymmetry || tSideSet->domain_type() == DomainType::AirSymmetry )
                        {
                            for ( Facet * tFacet : tSideSet->facets() )
                            {
                                tFacet->flag_nodes( 0 );
                            }
                        }
                    }
                    for ( fem::Dof * tDof : tField->dofs() )
                    {
                        if ( tDof->mesh_basis()->is_flagged( 1 ) )
                        {
                            tDof->fix( 1.0 );
                        }
                        else if ( tDof->mesh_basis()->is_flagged( 0 ) )
                        {
                            tDof->fix( 0.0 );
                        }
                    }
                }

                tField->initialize();
                tField->compute_jacobian();
                tField->solve();

                if ( mCommRank == 0 )
                {
                    mMesh->update_node_indices();
                    Cell< Node * > & tNodes = mMesh->nodes();

                    const Vector< real > & tPhi = mMesh->field_data( "phi" );

                    Cell< std::pair< Node *, real > > tPairs ;
                    tPairs.reserve( tNodes.size() );
                    for ( Node * tNode : tNodes )
                    {
                        tPairs.push( std::pair( tNode, tPhi( tNode->index() ) ) );
                    }

                    std::sort( tPairs.begin(), tPairs.end(),
                        []( std::pair< Node *, real > tA, std::pair< Node *, real > tB ) { return tA.second < tB.second; } );

                    index_t tCount = 0 ;
                    for ( auto tPair : tPairs )
                    {
                        tNodes( tCount++ ) = tPair.first ;
                    }
                    mMesh->update_node_indices();
                }

        }

    }
}
