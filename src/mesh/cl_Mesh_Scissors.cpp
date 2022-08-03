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
                    aMesh->create_edges();
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
                // first we create the new sheets
                if( mTapes.size() > 0 )
                {
                    this->check_compatibility( true );

                    this->collect_cut_data( true );
                    this->count_nodes();
                    this->compute_max_node_id();
                    this->duplicate_nodes( true );

                    this->count_facets();
                    this->compute_max_element_id();
                    // this->duplicate_facets();

                    this->relink_elements( true );
                    this->relink_original_facets( true );

                    // this->create_connectors( true ); // < -- no connector elements for tapes!

                    mMesh->unfinalize();

                    mMesh->finalize();
                }

                // second, we create the cuts
                if( mCuts.size() > 0 )
                {
                    this->check_compatibility( false );

                    this->collect_cut_data( false );
                    this->count_nodes();
                    this->compute_max_node_id();
                    this->duplicate_nodes( false );

                    this->count_facets();
                    this->compute_max_element_id();
                    // this->duplicate_facets();

                    this->relink_elements( false );
                    this->relink_original_facets( false );

                    this->create_connectors( false );

                    mMesh->unfinalize();

                    mMesh->finalize();
                }
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
            mData.set_size( 4, tCount );
            tCount = 0 ;

            for( scissors::CutData * tCut : tCuts )
            {
                index_t tNumSideSets = tCut->num_sidesets();
                for( index_t k=0; k<tNumSideSets; ++k )
                {
                    const Vector< id_t > & tData = tCut->data( k );

                    mData( 0, tCount ) = tCut->id();
                    mData( 1, tCount ) = tData( 0 );  // sideset
                    mData( 2, tCount ) = tData( 1 );  // plus
                    mData( 3, tCount ) = tData( 2 );  // minus

                    ++tCount ;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Scissors::count_nodes()
        {
            mMesh->unflag_all_nodes() ;
            index_t tNumSideSets = mData.n_cols() ;

            for( index_t k=0 ; k<tNumSideSets; ++k )
            {
                mMesh->sideset( mData( 1, k ) )->flag_all_nodes() ;
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
        Scissors::count_facets()
        {
            mMesh->unflag_all_facets();
            index_t tNumSideSets = mData.n_cols() ;

            for( index_t k=0 ; k<tNumSideSets; ++k )
            {
                mMesh->sideset( mData( 1, k ) )->flag_all_facets() ;
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
            Map< id_t, mesh::Node * > & tNodeMap = aTapeMode ? mNodeMapTapes : mNodeMapCuts ;

            index_t tNumNodes = mMesh->number_of_nodes() ;

            // make backup copy of original node list
            Cell< mesh::Node * > tOldNodes(
                    tNumNodes,
                    nullptr );

            index_t tCount = 0 ;

            // populate temporary container
            for( mesh::Node * tNode : mNodes )
            {
                tOldNodes( tCount++ ) = tNode ;
            }

            // resize array
            mNodes.clear() ;
            mNodes.set_size( tNumNodes + mNumberOfNodes, nullptr );

            // copy original nodes back
            tCount = 0 ;
            for( mesh::Node * tNode : tOldNodes )
            {
                mNodes( tCount++ ) = tNode ;
            }

            if( ! aTapeMode )
            {
                Matrix< id_t > & tNodeTable = mMesh->node_cut_table();
                tNodeTable.set_size( 2, mNumberOfNodes );

                // create a new counter
                index_t tDuplicateCount = 0;

                // crate the node duplicates
                for ( mesh::Node * tNode: tOldNodes )
                {
                    if ( tNode->is_flagged())
                    {
                        mesh::Node * tDuplicate = new Node(
                                mMaxNodeID++,
                                tNode->x(),
                                tNode->y(),
                                tNode->z());

                        tNodeMap[ tNode->id() ] = tDuplicate;

                        // if tapes exist, this node might already have been duplicated
                        // we use this trick to correct the cut
                        if ( mNodeMapTapes.key_exists( tNode->id()))
                        {
                            mNodeMapCuts[ mNodeMapTapes( tNode->id())->id() ] = tDuplicate;
                        }

                        // remember relation between original and duplicate
                        tNodeTable( 0, tDuplicateCount ) = tNode->id();
                        tNodeTable( 1, tDuplicateCount ) = tDuplicate->id();
                        ++tDuplicateCount;

                        // add new node to memory
                        mNodes( tCount++ ) = tDuplicate;
                    }
                }

                BELFEM_ASSERT( tDuplicateCount == mNumberOfNodes, "Internal Error in Cut routine: number of nodes does not match" );


            }
            else
            {
                // crate the node duplicates
                for ( mesh::Node * tNode: tOldNodes )
                {
                    if ( tNode->is_flagged())
                    {
                        mesh::Node * tDuplicate = new Node(
                                mMaxNodeID++,
                                tNode->x(),
                                tNode->y(),
                                tNode->z() );

                        tNodeMap[ tNode->id() ] = tDuplicate;

                        // add new node to memory
                        mNodes( tCount++ ) = tDuplicate;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::duplicate_facets( const bool aTapeMode )
        {
            Map< id_t, mesh::Node * > & tNodeMap = aTapeMode ? mNodeMapTapes : mNodeMapCuts ;

            index_t tNumFacets = mMesh->number_of_facets() ;

            // make backup copy of original node list
            Cell< mesh::Facet * > tOldFacets(
                    mMesh->number_of_facets(),
                    nullptr );

            index_t tCount = 0 ;
            for( mesh::Facet * tFacet : mFacets )
            {
                tOldFacets( tCount++ ) = tFacet ;
            }

            mFacets.clear() ;
            tCount = tNumFacets + mNumberOfFacets ;

            mFacets.set_size( tCount, nullptr );

            tCount = 0 ;
            for( mesh::Facet * tFacet : tOldFacets )
            {
                mFacets( tCount++ ) = tFacet ;
            }

            // the factory for elements
            ElementFactory tFactory ;

            for( mesh::Facet * tFacet : tOldFacets )
            {
                if( tFacet->is_flagged() )
                {
                    BELFEM_ASSERT( tFacet->master() != nullptr, "Facet %lu does not have a master element",
                                  ( long unsigned int ) tFacet->id() );

                    BELFEM_ASSERT( tFacet->slave() != nullptr, "Facet %lu does not have a slave element",
                                  ( long unsigned int ) tFacet->id() );

                    // create duplicate element
                    mesh::Element * tElement  = tFactory.create_lagrange_element(
                            tFacet->element()->type() , mMaxElementID++ );

                    // link element with duplicated nodes
                    uint tNumNodes = tFacet->element()->number_of_nodes() ;

                    for( uint k=0; k<tNumNodes; ++k )
                    {
                        tElement->insert_node( tNodeMap( tFacet->element()->node( k )->id() ), k );
                    }

                    // create new facet
                    mesh::Facet * tDuplicate = new mesh::Facet( tElement );

                    tDuplicate->set_master( tFacet->master(), tFacet->master_index(), false );
                    tDuplicate->set_slave( tFacet->slave(), tFacet->slave_index() );

                    mFacetMap[ tFacet->id() ] = tDuplicate ;
                    mFacets( tCount++ ) = tDuplicate ;
                }
            }

            BELFEM_ERROR( false, "TODO: New facets must be stored in new sidesets. Otherwise: Memory leak");
        }

//------------------------------------------------------------------------------

        void
        Scissors::relink_original_facets( const bool aTapeMode )
        {
            Map< id_t, mesh::Node * > & tNodeMap = aTapeMode ? mNodeMapTapes : mNodeMapCuts ;

            mMesh->unflag_all_facets() ;
            mMesh->unflag_all_elements() ;

            // now we flag all elements that are connected to nodes that are to be replaced
            for( Element * tElement : mElements )
            {
                uint tNumNodes = tElement->number_of_nodes() ;

                for( uint k=0; k<tNumNodes; ++k )
                {
                    if( tElement->node( k )->is_flagged() )
                    {
                        tElement->flag() ;
                        break ;
                    }
                }
            }

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
                                tElement->insert_node( tNodeMap( tNode->id()), k );
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
                                        tEdge->insert_node( tNodeMap( tNode->id()), k );
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
            Map< id_t, mesh::Node * > & tNodeMap = aTapeMode ? mNodeMapTapes : mNodeMapCuts ;

            index_t tNumSideSets = mData.n_cols() ;


            uint tDimension = mMesh->number_of_dimensions() ;
            for( index_t s=0 ; s<tNumSideSets; ++s )
            {
                // unflag all nodes on the mesh
                mMesh->unflag_all_nodes() ;

                // unflag all elements on the mesh
                mMesh->unflag_all_elements() ;

                uint tNumBlocks   = mMesh->number_of_blocks() ;

                id_t tPlus  = mData( 2, s );
                id_t tMinus = mData( 3, s );


                // get the sideset
                SideSet * tSideSet =  mMesh->sideset( mData( 1, s ) );

                // find out which blocks the plus block touches
                if( ! aTapeMode )
                {
                    Cell< Node * > & tNodes = tSideSet->nodes() ;

                    for( Node * tNode : tNodes )
                    {
                        tNode->flag();
                        uint tNumElems = tNode->number_of_elements();
                        for ( uint e = 0; e < tNumElems; ++e )
                        {
                            tNode->element( e )->flag();
                        }
                    }
                }
                else
                {
                    // get elements on plus block
                    Cell< Element * > & tPlusElements
                        = mMesh->block( tPlus )->elements() ;

                    // flag all corner nodes
                    for( Element * tElement : tPlusElements )
                    {
                        tElement->flag_corner_nodes() ;
                    }

                    Vector< index_t > tBlockCount( mMesh->number_of_blocks(), 0 );

                    // loop over all blocks
                    for( uint b=0; b<tNumBlocks; ++b )
                    {
                        // get block
                        Block * tBlock = mMesh->blocks()( b );

                        // check if this is neither the minus nor the plus block
                        if( tBlock->id() != tMinus )
                        {
                            // get element container
                            Cell< Element * > & tElements = tBlock->elements() ;

                            uint tCount  ;

                            // loop over all elements
                            for( Element * tElement : tElements )
                            {
                                tCount = 0 ;

                                // count flagged nodes
                                for( uint k=0; k<tElement->number_of_corner_nodes(); ++k )
                                {
                                    if( tElement->node( k )->is_flagged() )
                                    {
                                        // increment counter
                                        ++tCount ;
                                    }
                                }

                                // check if element touches
                                if( tCount >= tDimension )
                                {
                                    // set flag for this block
                                    tBlockCount( b ) = 1 ;
                                    break ;
                                }
                            }
                        } // end if block is not minus
                    } // end loop over all blocks

                    // create id list with touching blocks

                    Vector< id_t > tTouchBlocks( sum( tBlockCount ), 0 );
                    uint tCount = 0 ;
                    for( uint b=0; b<tNumBlocks; ++b )
                    {
                        if( tBlockCount( b ) == 1 )
                        {
                            tTouchBlocks( tCount++ ) = mMesh->blocks()(b)->id() ;
                        }
                    }
                    tNumBlocks = tCount ;

                    // unflag all nodes on the mesh
                    mMesh->unflag_all_nodes() ;

                    // flag all nodes on this sideset
                    Cell< Node * > & tNodes = tSideSet->nodes() ;
                    for( Node * tNode : tNodes )
                    {
                        tNode->flag();
                    }

                    // loop over all blocks of interest and flag the elements that are connected
                    // to this sideset
                    for( uint b=0; b<tNumBlocks; ++b )
                    {
                        Cell< Element * > & tElements = mMesh->block( tTouchBlocks( b ) )->elements() ;

                        // loop over all elements
                        for( Element * tElement : tElements )
                        {
                            for( uint k=0; k<tElement->number_of_nodes(); ++k )
                            {
                                if( tElement->node( k )->is_flagged() )
                                {
                                    // flag this element
                                    tElement->flag() ;

                                    // continue with next element
                                    break ;
                                }
                            }
                        }
                    }
                } // end tape mode

                // get the element container
                Cell< Element * > & tElements = aTapeMode ? mMesh->elements() : mMesh->block( tPlus )->elements() ;

                // loop over all elements in plus block
                for( Element * tElement : tElements )
                {
                    // check if element has been selected
                    if( tElement->is_flagged() )
                    {
                        // loop over all nodes on element
                        for( uint i=0; i<tElement->number_of_nodes(); ++i )
                        {
                            // grab node
                            Node * tNode = tElement->node( i );

                            // check if node sits on sideset
                            if( tNode->is_flagged() )
                            {
                                // replace node
                                tElement->insert_node( tNodeMap( tNode->id() ), i );
                            }
                        }

                        // relink edges of element
                        if( tElement->has_edges() )
                        {
                            for( uint e=0; e<tElement->number_of_edges(); ++e )
                            {
                                // grab edge
                                Edge * tEdge = tElement->edge( e );

                                // loop over all nodes of this edge
                                for( uint i=0; i<tEdge->number_of_nodes(); ++i )
                                {
                                    // grab node
                                    Node * tNode = tEdge->node( i );

                                    // check if node is to be replaced
                                    if ( tNode->is_flagged() )
                                    {
                                        tEdge->insert_node( tNodeMap( tNode->id()), i );
                                    }
                                } // end node loop
                            } // end edge loop
                        } // end edges exist
                    } // end if element is flagged
                } // end loop over all elements
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::create_connectors( const bool aTapeMode )
        {
            Map< id_t, mesh::Node * > & tNodeMap = aTapeMode ? mNodeMapTapes : mNodeMapCuts ;

            ElementFactory tFactory ;

            index_t tNumSideSets = mData.n_cols() ;

            mConnectorSetIDs.set_size( tNumSideSets );

            mSideSetToCutIDs.clear() ;

            for( index_t k=0 ; k<tNumSideSets; ++k )
            {
                // get sideset
                SideSet * tSideSet =  mMesh->sideset( mData( 1, k ) );

                // number of nodes on this sideset
                Cell< Node * > & tNodes = tSideSet->nodes();

                // create a new block
                SideSet * tCut = new SideSet( mData( 0, k ), tNodes.size() );

                Cell< Facet * > & tFacets = tCut->facets() ;

                index_t tCount = 0 ;

                // get id of master block
                id_t tPlusID = mData( 2, k );

                // needed for a bugfix that the connectors sit correct on cuts that intersect tapes
                bool tSwitch = ! aTapeMode && mPlus.key_exists( tPlusID );

                for( Node * tNode : tNodes )
                {
                    // create a new element
                    Element * tElement = tFactory.create_lagrange_element( ElementType::LINE2, ++mMaxElementID );

                    if( tSwitch && mNodeMapTapes.key_exists( tNode->id())  )
                    {
                        // insert negative side
                        tElement->insert_node( mNodeMapTapes( tNode->id() ), 0 );

                    }
                    else
                    {
                        // insert negative side
                        tElement->insert_node( tNode, 0 );
                    }


                    // insert positive side
                    tElement->insert_node( tNodeMap( tNode->id() ), 1 );
                    tElement->set_block_id( tPlusID );

                    // add facet to block
                    tFacets( tCount++ ) = new Facet( tElement );
                }

                mSideSetToCutIDs[ tSideSet->id() ] = tCut->id() ;

                // remember ID of new block
                mConnectorSetIDs( k ) = tCut->id() ;

                // add sideset to mesh
                mMesh->cuts().push( tCut );
            }
        }

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

//------------------------------------------------------------------------------
    }
}