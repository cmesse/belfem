//
// Created by christian on 9/15/21.
//

#include "cl_Mesh_Scissors.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "cl_Element.hpp"
#include "cl_Facet.hpp"
#include "cl_Matrix.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace mesh
    {
        namespace scissors
        {

//------------------------------------------------------------------------------

            CutData::CutData(
                     const id_t aCutID,
                     const id_t aSideSetID,
                     const id_t aPlusBlock,
                     const id_t aMinusBlock ) :
                         mID( aCutID )
            {
                Vector< id_t > tData( 3 );
                tData( 0 ) = aSideSetID ;
                tData( 1 ) = aPlusBlock ;
                tData( 2 ) = aMinusBlock ;
                mData.set_size( 1, tData );
            }

//------------------------------------------------------------------------------

            CutData::CutData(
                    const id_t aCutID,
                    const Vector< id_t > & aSideSetIDs,
                    const id_t aPlusBlock,
                    const id_t aMinusBlock ) :
                         mID( aCutID )
            {
                Vector< id_t > tData( 3, 0 );
                tData( 1 ) = aPlusBlock ;
                tData( 2 ) = aMinusBlock ;

                index_t tNumSidesets = aSideSetIDs.length();
                mData.set_size( tNumSidesets, tData );

                for( index_t k=0; k<tNumSidesets; ++k )
                {
                    mData( k )( 0 ) = aSideSetIDs( k );
                }
            }

//------------------------------------------------------------------------------

            CutData::CutData( const id_t aCutID,
                     const Vector< id_t > & aSideSetIDs,
                     const Vector< id_t > & aPlusBlocks,
                     const Vector< id_t > & aMinusBlocks ):
                     mID( aCutID )
            {
                BELFEM_ASSERT( aSideSetIDs.length() == aPlusBlocks.length(),
                              "Length of PlusBlocks does not match ( is %u, expect %u )",
                              ( unsigned int ) aPlusBlocks.length(),
                              ( unsigned int ) aSideSetIDs.length() );

                BELFEM_ASSERT( aSideSetIDs.length() == aMinusBlocks.length(),
                              "Length of MinusBlocks does not match ( is %u, expect %u )",
                              ( unsigned int ) aMinusBlocks.length(),
                              ( unsigned int ) aSideSetIDs.length() );

                index_t tNumSidesets = aSideSetIDs.length();
                Vector< id_t > tData( 3, 0 );
                mData.set_size( tNumSidesets, tData );

                for( index_t k=0; k<tNumSidesets; ++k )
                {
                    mData( k )( 0 ) = aSideSetIDs( k );
                    mData( k )( 1 ) = aPlusBlocks( k );
                    mData( k )( 2 ) = aMinusBlocks( k );
                }

            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        Scissors::Scissors( Mesh * aMesh, const Vector< id_t > & aAirBlockIDs  ) :
            mMesh( aMesh ),
            mMyRank( comm_rank() ),
            mNumberOfOriginalNodes( aMesh->number_of_nodes() ),
            mNodes( aMesh->nodes() ),
            mFacets( aMesh->facets() ),
            mElements( aMesh->elements() ),
            mAirBlockIDs( aAirBlockIDs )
        {
            if( mMyRank == mMesh->master() )
            {

                // create the edges if they don't exist already
                if ( !aMesh->edges_exist() )
                {
                    BELFEM_ERROR( false, "Edges for this mesh have not been created");
                    // aMesh->create_edges();
                }
                this->compute_max_sideset_id();

            }
        }

//------------------------------------------------------------------------------

        Scissors::~Scissors()
        {
            for( scissors::CutData * tCut : mCuts )
            {
                delete tCut ;
            }
            for( scissors::CutData * tTape : mTapes )
            {
                delete tTape ;
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::cut(
                const id_t aSidesetID,
                const id_t aPlusBlock,
                const id_t aMinusBlock,
                const bool aIsSheet )
        {
            if( mMyRank == mMesh->master() )
            {
                if( aIsSheet )
                {
                    mTapes.push( new scissors::CutData(
                            ++mMaxSideSetID,
                            aSidesetID,
                            aPlusBlock,
                            aMinusBlock ));
                    ++mNumberOfTapes;
                }
                else
                {
                    mCuts.push( new scissors::CutData(
                            ++mMaxSideSetID,
                            aSidesetID,
                            aPlusBlock,
                            aMinusBlock ));
                    ++mNumberOfCuts;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::cut( const Vector< id_t > & aSidesetIDs,
             const id_t aPlusBlock,
             const id_t aMinusBlock,
             const bool aIsSheet )
        {
            if( mMyRank == mMesh->master() )
            {
                if( aIsSheet )
                {
                    mTapes.push( new scissors::CutData(
                            ++mMaxSideSetID,
                            aSidesetIDs,
                            aPlusBlock,
                            aMinusBlock ));
                    ++mNumberOfTapes;
                }
                else
                {
                    mCuts.push( new scissors::CutData(
                            ++mMaxSideSetID,
                            aSidesetIDs,
                            aPlusBlock,
                            aMinusBlock ));
                    ++mNumberOfCuts;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::cut( const Vector< id_t > & aSidesetIDs,
             const Vector< id_t > & aMinusBlocks,
             const Vector< id_t > & aPlusBlocks,
             const bool aIsSheet )
        {
            if( mMyRank == mMesh->master() )
            {
                if( aIsSheet )
                {
                    mTapes.push( new scissors::CutData(
                            ++mMaxSideSetID,
                            aSidesetIDs,
                            aPlusBlocks,
                            aMinusBlocks ));
                    ++mNumberOfTapes;
                }
                else
                {
                    mCuts.push( new scissors::CutData(
                            ++mMaxSideSetID,
                            aSidesetIDs,
                            aPlusBlocks,
                            aMinusBlocks ));
                    ++mNumberOfCuts;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        Scissors::compute_max_sideset_id()
        {
            mMaxSideSetID = 0 ;
            for( SideSet * tSideSet : mMesh->sidesets() )
            {
                mMaxSideSetID = tSideSet->id() > mMaxSideSetID ?
                                tSideSet->id() : mMaxSideSetID ;

            }
            ++mMaxSideSetID ;
        }

//------------------------------------------------------------------------------

        void
        Scissors::compute_max_node_id()
        {
            mMaxNodeID = 0 ;
            for( mesh::Node * tNode : mNodes )
            {
                mMaxNodeID = tNode->id() > mMaxNodeID ?
                        tNode->id() : mMaxNodeID ;
            }
            ++mMaxNodeID;
        }

//------------------------------------------------------------------------------

        void
        Scissors::compute_max_element_id()
        {
            mMaxElementID = 0 ;

            for( mesh::Element * tElement : mElements )
            {
                mMaxElementID = tElement->id() > mMaxElementID ?
                        tElement->id() : mMaxElementID ;
            }

            for( mesh::Facet * tFacet : mFacets )
            {
                mMaxElementID = tFacet->id() > mMaxElementID ?
                        tFacet->id() : mMaxElementID ;
            }
            ++mMaxElementID;
        }

//------------------------------------------------------------------------------

        void
        Scissors::finalize()
        {
            // only do something on master
            if( mMyRank == mMesh->master() )
            {
                mMesh->unflag_all_nodes() ;
                mMesh->unflag_all_elements() ;

                // first we create the new sheets
                if( mTapes.size() > 0 )
                {
                    this->check_compatibility( true );

                    this->collect_cut_data( true );
                    this->count_nodes( true );
                    this->compute_max_node_id();
                    this->duplicate_nodes( true );
                    this->create_cut_table( true ) ;
                    this->count_facets( true );
                    this->compute_max_element_id();
                    // this->duplicate_facets();

                    this->relink_elements( true );
                    this->relink_facets();
                }

                // second, we create the cuts
                if( mCuts.size() > 0 )
                {
                    this->check_compatibility( false );

                    this->collect_cut_data( false );
                    this->count_nodes( false );
                    this->compute_max_node_id();
                    this->duplicate_nodes( false );
                    this->create_cut_table( false ) ;

                    this->count_facets( false );
                    this->compute_max_element_id();

                    this->relink_elements( false );
                    this->create_connectors();
                }

                mMesh->unfinalize();
                mMesh->finalize() ;

            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::collect_cut_data( const bool aTapeMode )
        {
            // select which cell to use
            Cell< scissors::CutData * > & tCuts = aTapeMode ? mTapes : mCuts ;

            // count total cuts
            index_t tCount = 0 ;
            for( scissors::CutData * tCut : tCuts )
            {
                tCount += tCut->num_sidesets();
            }

            // flatten table
            Matrix< id_t > & tData = aTapeMode ? mTapeData : mCutData ;

            tData.set_size( 4, tCount );
            tCount = 0 ;

            for( scissors::CutData * tCut : tCuts )
            {
                index_t tNumSideSets = tCut->num_sidesets();
                for( index_t k=0; k<tNumSideSets; ++k )
                {
                    const Vector< id_t > & tCutData = tCut->data( k );

                    tData( 0, tCount ) = tCut->id();
                    tData( 1, tCount ) = tCutData( 0 );  // sideset
                    tData( 2, tCount ) = tCutData( 1 );  // plus
                    tData( 3, tCount ) = tCutData( 2 );  // minus

                    ++tCount ;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Scissors::count_nodes( const bool aTapeMode )
        {

            Matrix< id_t > & tData = aTapeMode ? mTapeData : mCutData ;
            index_t tNumSideSets = tData.n_cols() ;

            // fix facet links
            if( ! aTapeMode )
            {
                for( index_t s=0 ; s<tNumSideSets; ++s )
                {
                    Cell< Facet * > & tFacets = mMesh->sideset( tData( 1, s ) )->facets() ;

                    uint tNumNodes = mesh::number_of_nodes( mMesh->sideset( tData( 1, s ) )->element_type() );

                    for( Facet * tFacet : tFacets )
                    {
                        if( tFacet->master()->is_flagged() && tFacet->slave()->is_flagged() )
                        {
                            for ( uint k = 0; k < tNumNodes; ++k )
                            {
                                if ( tFacet->node( k )->is_flagged())
                                {
                                    tFacet->element()->insert_node( mDuplicateNodes( tFacet->node( k )->index()), k );
                                }
                            }
                        }
                    }
                }

            }

            mMesh->unflag_all_nodes() ;

            for( index_t s=0 ; s<tNumSideSets; ++s )
            {
                Cell< Facet * > & tFacets = mMesh->sideset( tData( 1, s ) )->facets() ;

                for( Facet * tFacet : tFacets )
                {

                    tFacet->flag_nodes() ;
                }
            }

            mNumberOfNodes = 0 ;

            for( mesh::Node * tNode : mNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++mNumberOfNodes ;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Scissors::count_facets( const bool aTapeMode  )
        {
            Matrix< id_t > & tData = aTapeMode ? mTapeData : mCutData ;

            mMesh->unflag_all_facets();
            index_t tNumSideSets = tData.n_cols() ;

            for( index_t k=0 ; k<tNumSideSets; ++k )
            {
                mMesh->sideset( tData( 1, k ) )->flag_all_facets() ;
            }

            mNumberOfFacets = 0 ;

            for( mesh::Facet * tFacet : mFacets )
            {
                if( tFacet->is_flagged() )
                {
                    ++mNumberOfFacets ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::duplicate_nodes( const bool aTapeMode )
        {
            index_t tNumNodes = mMesh->number_of_nodes() ;

            // backup original nodes on mesh
            Cell< mesh::Node * > tNodes ;
            tNodes.vector_data() = std::move( mMesh->nodes().vector_data() );

            // resize array
            mNodes.clear() ;
            mNodes.set_size( tNumNodes + mNumberOfNodes, nullptr );

            // copy original nodes back
            index_t tAllCount = 0 ;
            for( mesh::Node * tNode : tNodes )
            {
                tNode->set_index( gNoIndex );
                mNodes( tAllCount++ ) = tNode ;
            }

            mOriginalNodes.set_size( mNumberOfNodes, nullptr );
            mDuplicateNodes.set_size( mNumberOfNodes, nullptr );

            // crate the node duplicates
            index_t tCount = 0 ;

            // crate the node duplicates
            for ( mesh::Node * tNode: tNodes )
            {
                if ( tNode->is_flagged())
                {
                    // create duplicate
                    mesh::Node * tDuplicate = new Node(
                            mMaxNodeID++,
                            tNode->x(),
                            tNode->y(),
                            tNode->z());

                    // set index of original node
                    tNode->set_index( tCount );

                    // set index of duplicate ( needed for connector creation )
                    tDuplicate->set_index( tCount + mNumberOfNodes );

                    // store old and new nodes in special containers
                    mOriginalNodes( tCount ) = tNode ;
                    mDuplicateNodes( tCount++ ) = tDuplicate ;

                    // also add new node to mesh
                    mNodes( tAllCount++ ) = tDuplicate;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::create_cut_table( const bool aTapeMode )
        {
            Matrix< id_t > & tTable = aTapeMode ?
                    mMesh->node_tape_table() : mMesh->node_cut_table() ;

            index_t tNumNodes = mOriginalNodes.size() ;

            tTable.set_size( 2, tNumNodes );

            for( index_t k=0; k<tNumNodes; ++k )
            {
                tTable( 0 , k ) = mOriginalNodes( k )->id() ;
                tTable( 1, k )  = mDuplicateNodes( k )->id() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::relink_facets()
        {
            // each facet has a master and a slave element
            // for interface elements, the master is always solid (superconductor, coil or iron)
            // and the slave is always air

            // loop over al facets
            for( mesh::Facet * tFacet : mFacets )
            {
                if( tFacet->has_slave() )
                {
                    // check if the facet needs to be relinked
                    if ( tFacet->master()->is_flagged() && !tFacet->slave()->is_flagged())
                    {
                        // get element of facet
                        Element * tElement = tFacet->element();

                        // number of nodes on the element
                        uint tNumNodes = tElement->number_of_nodes();

                        // relink facet
                        for ( uint k = 0; k < tNumNodes; ++k )
                        {
                            // grab node
                            Node * tNode = tElement->node( k );

                            // check if node is to be replaced
                            if ( tNode->is_flagged() )
                            {
                                tElement->insert_node( mDuplicateNodes( tNode->index()), k );
                            }
                        } // end node loop

                        // relink edges of facet
                        if( tElement->has_edges() )
                        {
                            for( uint e=0; e<tElement->number_of_edges(); ++e )
                            {
                                // grab edge
                                Edge * tEdge = tElement->edge( e );

                                // loop over all nodes of this edge
                                for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                                {
                                    // grab node
                                    Node * tNode = tEdge->node( k );

                                    // check if node is to be replaced
                                    if ( tNode->is_flagged() )
                                    {
                                        tEdge->insert_node( mDuplicateNodes( tNode->index()), k );
                                    }
                                } // end node loop
                            } // end edge loop
                        } // end edges exist
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::relink_elements( const bool aTapeMode )
        {
            Matrix< id_t > & tData = aTapeMode ? mTapeData : mCutData ;

            // unflag all nodes on the mesh
            mMesh->unflag_all_nodes() ;

            // unflag all elements on the mesh
            mMesh->unflag_all_elements() ;

            // flag the nodes that are to be duplicated
            for( Node * tNode : mOriginalNodes )
            {
                tNode->flag() ;
            }

            uint tNumCuts = tData.n_cols() ;

            for( uint c=0; c<tNumCuts; ++c )
            {
                Block * tBlock = mMesh->block( tData( 2, c ) );

                // grab elements from block
                Cell< Element * > & tElements = tBlock->elements() ;

                uint tNumNodes = mesh::number_of_nodes( tBlock->element_type() );

                // loop over all elements on plus block
                for( Element * tElement  : tElements )
                {
                    for( uint k=0; k<tNumNodes; ++k )
                    {
                        if( tElement->node( k )->is_flagged() )
                        {
                            // relink the node
                            tElement->insert_node( mDuplicateNodes( tElement->node( k )->index() ), k );

                            // make sure that the element is flagged
                            tElement->flag() ;
                        }
                    }
                }
            } // end loop over all blocks

        }

//------------------------------------------------------------------------------
#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif
        void
        Scissors::create_connectors()
        {
            // the factory that we use to create the new elements
            ElementFactory tFactory ;

            index_t tNumSideSets = mCutData.n_cols() ;

            mConnectorSetIDs.set_size( tNumSideSets );

            index_t tNodeCount = 0 ;

            mMesh->unflag_all_nodes() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 1: count how many connectors are created in total per sideset
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            Vector< index_t > tNumConnectors( tNumSideSets, 0 );

            for( index_t s=0 ; s<tNumSideSets; ++s )
            {
                SideSet * tSideSet = mMesh->sideset( mCutData( 1, s )) ;

                uint tNumNodesPerFacet = number_of_nodes( tSideSet->element_type() ) ;

                // get facets
                Cell< Facet * > & tFacets = tSideSet->facets() ;

                // loop over all facets
                for( Facet * tFacet : tFacets )
                {
                    for( uint k=0; k<tNumNodesPerFacet; ++k )
                    {
                        Node * tNode = tFacet->node( k );
                        if( ! tNode->is_flagged() )
                        {
                            tNode->flag() ;
                            ++tNumConnectors( s );
                        }
                    }
                }
            }

            // sanity check
            BELFEM_ERROR( sum( tNumConnectors ) == mNumberOfNodes,
                          "Something went wrong while creating connectors for cuts" ) ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 2: create the connectors
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for( index_t s=0 ; s<tNumSideSets; ++s )
            {
                SideSet * tSideSet = mMesh->sideset( mCutData( 1, s ));

                // get facets from sideset
                Cell < Facet * > & tFacets = tSideSet->facets();

                // create a new sideset
                SideSet * tCut = new SideSet( mCutData( 0, s ), tNumConnectors( s ) );

                // get the connectors for this sideset
                Cell< Facet * > & tConnectors = tCut->facets() ;

                uint tNumNodesPerFacet = number_of_nodes( tSideSet->element_type());

                // initialize a counter
                index_t tCount = 0 ;

                // id of plus block
                id_t tPlus  = mCutData( 2, s );

                // loop over all facets
                for ( Facet * tFacet: tFacets )
                {
                    for( uint k=0; k<tNumNodesPerFacet; ++k )
                    {
                        Node * tNode = tFacet->node( k );
                        if( tNode->is_flagged() )
                        {
                            // create a new element
                            Element * tElement = tFactory.create_lagrange_element( ElementType::LINE2, ++mMaxElementID );

                            // link to negative side
                            tElement->insert_node( tNode, 0 );

                            // link to positive side
                            tElement->insert_node( mDuplicateNodes( tNode->index() ), 1 );

                            // we need this information for the dof linking (todo, can be done prettier)
                            tElement->set_block_id( tPlus );

                            // add new facet to  sideset
                            tConnectors( tCount++ ) = new Facet( tElement );

                            // unflag this node
                            tNode->unflag() ;
                        }
                    }
                }

                mSideSetToCutIDs[ tSideSet->id() ] = tCut->id() ;

                // remember ID of new block
                mConnectorSetIDs( s ) = tCut->id() ;

                // add cut to mesh
                mMesh->cuts().push( tCut );

            }
        }
#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
//------------------------------------------------------------------------------

        void
        Scissors::get_cut_ids( const Vector< id_t > & aSideSetIDs,
                                     Vector< id_t > & aCutIDs )
        {
            proc_t tCommSize = comm_size() ;
            proc_t tMyRank   = comm_rank() ;

            if( tMyRank == mMesh->master() )
            {
                index_t tCount = 0 ;

                // populate output
                if( aSideSetIDs.length() > 0 )
                {
                    aCutIDs.set_size( aSideSetIDs.length() );

                    for( id_t tSideSetID : aSideSetIDs )
                    {
                        aCutIDs( tCount++ ) = mSideSetToCutIDs( tSideSetID );
                    }
                }

                // communicate to other procs
                if( tCommSize > 1 )
                {
                    // create comm table
                    Vector< proc_t > tCommTable( tCommSize );

                    for( proc_t tRank = 0; tRank<tCommSize; ++tRank )
                    {

                        tCommTable( tRank ) = tRank ;
                    }

                    send_same( tCommTable, aCutIDs );
                }
            }
            else
            {
                receive( mMesh->master(), aCutIDs );
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::check_compatibility( const bool aTapeMode )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - -
            // step 1: collect all minus blocks from the cuts
            // - - - - - - - - - - - - - - - - - - - - - - - -

            Cell< scissors::CutData * > & tCuts = aTapeMode ? mTapes : mCuts ;

            if( tCuts.size() > 0 )
            {
                // count all minus blocks from cuts
                index_t tCount = 0;
                for ( scissors::CutData * tCut: tCuts )
                {
                    tCount += tCut->num_sidesets();
                }

                // collect all minus blocks from cuts
                Vector< id_t > tM( tCount );
                Vector< id_t > tP( tCount );
                tCount = 0;
                for ( scissors::CutData * tCut: tCuts )
                {
                    for ( index_t k = 0; k < tCut->num_sidesets(); ++k )
                    {
                        tP( tCount ) = tCut->data( k )( 1 );
                        tM( tCount ) = tCut->data( k )( 2 );
                        ++tCount;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - -
                // step 2: create unique maps
                // - - - - - - - - - - - - - - - - - - - - - - - -

                unique( tM );
                unique( tP );

                Map< id_t, index_t > tPlusMap ;
                Map< id_t, index_t > tMinusMap ;
                Map< id_t, index_t > & tPlus = aTapeMode ? mPlus : tPlusMap ;
                Map< id_t, index_t > & tMinus = aTapeMode ? mMinus : tMinusMap ;


                tCount = 0;
                for ( id_t tID: tM )
                {
                    tMinus[ tID ] = tCount++;
                }
                tCount = 0;
                for ( id_t tID: tP )
                {
                    tPlus[ tID ] = tCount++;
                }


                string tString = aTapeMode ? "tapes" : "cuts";

                // - - - - - - - - - - - - - - - - - - - - - - - -
                // step 3: check all tapa data
                // - - - - - - - - - - - - - - - - - - - - - - - -
                for ( scissors::CutData * tTape: tCuts )
                {
                    for ( index_t k = 0; k < tTape->num_sidesets(); ++k )
                    {
                        Vector< id_t > & tData = tTape->data( k );

                        id_t tA = tData( 1 );
                        id_t tB = tData( 2 );

                        BELFEM_ERROR( tPlus.key_exists( tA ) ^ tMinus.key_exists( tA ),
                                      "Incompatible cutting scheme: Block %lu can not be positive and negative at the same time when creating %s",
                                      ( long unsigned int ) tA, tString.c_str() );


                        BELFEM_ERROR( tPlus.key_exists( tB ) ^ tMinus.key_exists( tB ),
                                      "Incompatible cutting scheme: Block %lu can not be positive and negative at the same time when creating %s",
                                      ( long unsigned int ) tB, tString.c_str() );

                    }
                }
            }
        }

    }
}