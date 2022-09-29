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
#include "fn_max.hpp"

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
                     const id_t aMinusBlock,
                     const bool aIsTape  ) :
                         mID( aCutID ),
                         mIsTape( aIsTape )
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
                    const id_t aMinusBlock,
                    const bool aIsTape  ) :
                         mID( aCutID ),
                         mIsTape( aIsTape )
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
                     const Vector< id_t > & aMinusBlocks,
                     const bool aIsTape ):
                     mID( aCutID ),
                     mIsTape( aIsTape )
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
            for( scissors::CutData * tCutData : mCutsAndTapes )
            {
                delete tCutData ;
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::cut(
                const id_t aSidesetID,
                const id_t aPlusBlock,
                const id_t aMinusBlock,
                const bool aIsTape )
        {
            if( mMyRank == mMesh->master() )
            {
                mCutsAndTapes.push( new scissors::CutData(
                        ++mMaxSideSetID,
                        aSidesetID,
                        aPlusBlock,
                        aMinusBlock,
                        aIsTape ) );

                if( aIsTape )
                {
                    ++mNumberOfTapes;
                }
                else
                {
                    ++mNumberOfCuts;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::cut( const Vector< id_t > & aSidesetIDs,
             const id_t aPlusBlock,
             const id_t aMinusBlock,
             const bool aIsTape )
        {
            if( mMyRank == mMesh->master() )
            {
                mCutsAndTapes.push( new scissors::CutData(
                        ++mMaxSideSetID,
                        aSidesetIDs,
                        aPlusBlock,
                        aMinusBlock,
                        aIsTape ) );

                if( aIsTape )
                {
                    ++mNumberOfTapes;
                }
                else
                {
                    ++mNumberOfCuts;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::cut( const Vector< id_t > & aSidesetIDs,
             const Vector< id_t > & aMinusBlocks,
             const Vector< id_t > & aPlusBlocks,
             const bool aIsTape )
        {
            if( mMyRank == mMesh->master() )
            {
                mCutsAndTapes.push( new scissors::CutData(
                        ++mMaxSideSetID,
                        aSidesetIDs,
                        aPlusBlocks,
                        aMinusBlocks,
                        aIsTape ) );

                if( aIsTape )
                {
                    ++mNumberOfTapes;
                }
                else
                {
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
                this->compute_max_node_id() ;
                this->compute_max_element_id() ;
                this->compute_max_sideset_id() ;

                mMesh->unflag_all_nodes();
                mMesh->unflag_all_elements();
                mMesh->unflag_all_facets() ;
                this->collect_sidesets() ;
                this->collect_original_nodes() ;
                this->count_duplicate_nodes() ;
                this->create_node_duplicates();
                this->relink_elements();
                this->relink_facets() ;
                this->create_connectors() ;
                this->add_nodes_to_mesh() ;

                mMesh->unfinalize();
                mMesh->finalize() ;

            }

            // visualization hack
            /*mMesh->unflag_all_nodes() ;
            for( uint b=2; b<6; ++b )
            {
                for( Element * tElement : mMesh->block( b )->elements() )
                {
                    uint tNumNodes = tElement->number_of_nodes() ;
                    for( uint k=0; k<tNumNodes; ++k )
                    {
                        Node * tNode = tElement->node( k );

                        tElement->node( k )->set_coords( tNode->x(), tNode->y(), b * 0.001 );
                    }
                }
            }

            mMesh->save( "test.vtk");*/
            // don't forget to add new entities to the mesh!
            //exit( 0  );
        }

//------------------------------------------------------------------------------

        void
        Scissors::collect_cut_data( const bool aTapeMode )
        {

            // count total cuts
            index_t tCount = 0 ;
            for( scissors::CutData * tCut : mCutsAndTapes )
            {
                if ( tCut->is_tape() == aTapeMode )
                {
                    tCount += tCut->num_sidesets();
                }
            }

            // flatten table
            Matrix< id_t > & tData = aTapeMode ? mTapeData : mCutData ;

            tData.set_size( 4, tCount );
            tCount = 0 ;

            for( scissors::CutData * tCut : mCutsAndTapes )
            {
                if ( tCut->is_tape() == aTapeMode )
                {
                    index_t tNumSideSets = tCut->num_sidesets();
                    for ( index_t s = 0; s < tNumSideSets; ++s )
                    {
                        const Vector< id_t > & tCutData = tCut->data( s );

                        tData( 0, tCount ) = tCut->id();
                        tData( 1, tCount ) = tCutData( 0 );  // sideset
                        tData( 2, tCount ) = tCutData( 1 );  // plus
                        tData( 3, tCount ) = tCutData( 2 );  // minus

                        ++tCount;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::get_cut_ids( const Vector< id_t > & aSideSetIDs,
                                     Vector< id_t > & aCutIDs )
        {
            if( comm_rank() == mMesh->master() )
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
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::check_compatibility( const bool aTapeMode )
        {
            /*
            // - - - - - - - - - - - - - - - - - - - - - - - -
            // step 1: collect all minus blocks from the cuts
            // - - - - - - - - - - - - - - - - - - - - - - - -

            Cell< scissors::CutData * > & tCuts = aTapeMode ? mTapes : mCuts ;

            for ( scissors::CutData * tCut: tCuts )
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
            }*/
        }

//------------------------------------------------------------------------------

        void
        Scissors::collect_sidesets()
        {
            // counters
            index_t tS = 0 ;
            index_t tC = 0 ;
            index_t tT = 0 ;

            // count sidesets
            for ( scissors::CutData * tCutData : mCutsAndTapes )
            {
                tS += tCutData->num_sidesets() ;

                if( tCutData->is_tape() )
                {
                    tT += tCutData->num_sidesets() ;
                }
                else
                {
                    tC += tCutData->num_sidesets() ;
                }
            }

            // allocate memory
            mAllSideSets.set_size( tS, nullptr );
            mTapeSideSets.set_size( tT, nullptr );
            mCutSideSets.set_size( tC, nullptr );

            // reset counters
            tS = 0 ;
            tT = 0 ;
            tC = 0 ;
            for ( scissors::CutData * tCutData : mCutsAndTapes )
            {
                for( uint s = 0 ; s<tCutData->num_sidesets() ; ++s )
                {
                    // get sideset
                    mAllSideSets( tS++ ) = mMesh->sideset( tCutData->data( s )( 0 ) );
                }
                if( tCutData->is_tape() )
                {
                    for( uint s = 0 ; s<tCutData->num_sidesets() ; ++s )
                    {
                        // get sideset
                        mTapeSideSets( tT++ ) = mMesh->sideset( tCutData->data( s )( 0 ) );
                    }
                }
                else
                {
                    for( uint s = 0 ; s<tCutData->num_sidesets() ; ++s )
                    {
                        // get sideset
                        mCutSideSets( tC++ ) = mMesh->sideset( tCutData->data( s )( 0 ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::collect_original_nodes()
        {
            mMesh->unflag_all_nodes() ;
            mMesh->unflag_all_facets() ;

            // flag all nodes and facets relevant for this mesh
            for( SideSet * tSideSet : mAllSideSets )
            {
                for ( Facet * tFacet: tSideSet->facets() )
                {
                    tFacet->element()->flag_nodes();
                    tFacet->flag() ;
                }
            }

            // count nodes
            index_t tCount = 0 ;
            for( Node * tNode : mMesh->nodes() )
            {
                if( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ ) ;
                }
            }

            // collect nodes
            Cell< Node * > tNodes( tCount, nullptr );
            tCount = 0 ;
            for( Node * tNode : mMesh->nodes() )
            {
                if( tNode->is_flagged() )
                {
                    tNodes( tCount++ ) = tNode ;
                }
            }

            // count facets per node
            Vector< uint > tTapeCount( tCount, 0 );
            for( SideSet * tTape : mTapeSideSets )
            {
                for ( Facet * tFacet: tTape->facets() )
                {
                    for( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        ++tTapeCount( tFacet->node( k )->index() );
                    }
                }
            }

            // count tapes per node
            Vector< uint > tCutCount( tCount, 0 );
            for( SideSet * tCut : mCutSideSets )
            {
                for ( Facet * tFacet: tCut->facets() )
                {
                    for( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        ++tCutCount( tFacet->node( k )->index() );
                    }
                }
            }


            // next, we unflag side nodes
            for ( Node * tNode: tNodes )
            {
                if (    tTapeCount( tNode->index()) == 1      // all nodes that have one tape facet
                     && tCutCount( tNode->index())  == 0 )     // but no cuts facet
                {
                    tNode->unflag();
                }
            }

            // but: center nodes are actualy needed!
            for( SideSet * tSideSet : mAllSideSets )
            {
                uint tNumNodes = mesh::number_of_nodes( tSideSet->element_type() );
                uint tFirstInnerNode = this->inside_node( tSideSet->element_type() );

                if( tFirstInnerNode < tNumNodes )
                {
                    for ( Facet * tFacet: tSideSet->facets())
                    {
                        for( uint k=tFirstInnerNode; k<tNumNodes; ++k )
                        {
                            tFacet->node( k )->flag() ;
                        }
                    }
                }
            }


            // count
            tCount = 0 ;
            for( Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }

            mNumberOfNodes = tCount ;
            mOriginalNodes.set_size( tCount, nullptr );
            tCount = 0 ;
            for( Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    mOriginalNodes( tCount++ ) = tNode ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::count_duplicate_nodes()
        {
            mNumTapesPerNode.set_size( mNumberOfNodes, 0 );
            mNumCutsPerNode.set_size( mNumberOfNodes, 0 );
            mNodeCounter.set_size( mNumberOfNodes );

            // count cuts
            for( SideSet * tSideSet : mCutSideSets )
            {
                mNodeCounter.fill( 0 );
                for( Facet * tFacet : tSideSet->facets() )
                {
                    uint tNumNodes = tFacet->number_of_nodes() ;
                    for( uint k=0; k< tNumNodes; ++k )
                    {
                        if( tFacet->node( k )->is_flagged() )
                        {
                            mNodeCounter( tFacet->node( k )->index() ) = 1;
                        }
                    }
                }
                mNumCutsPerNode += mNodeCounter ;
            }

            // count tapes
            for( SideSet * tSideSet : mTapeSideSets )
            {
                mNodeCounter.fill( 0 );
                for( Facet * tFacet : tSideSet->facets() )
                {
                    uint tNumNodes = tFacet->number_of_nodes() ;
                    for( uint k=0; k< tNumNodes; ++k )
                    {
                        if( tFacet->node( k )->is_flagged() )
                        {
                            mNodeCounter( tFacet->node( k )->index()) = 1;
                        }
                    }
                }
                mNumTapesPerNode += mNodeCounter ;
            }

            // list of plus blocks per node
            mPlusBlocksTable.set_size( mNumberOfNodes, {} );

            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                mPlusBlocksTable( k ).set_size( mNumTapesPerNode( k ) + mNumCutsPerNode( k ) );
            }

            // additional counter
            Vector< uint > tCount( mNumberOfNodes, 0 );

            // populate container
            for ( scissors::CutData * tCut : mCutsAndTapes )
            {
                for ( index_t s = 0; s < tCut->num_sidesets(); ++s )
                {
                    mNodeCounter.fill( 0 );
                    Cell< Facet * > & tFacets = mMesh->sideset( tCut->data( s )( 0 )  )->facets();

                    for( Facet * tFacet : tFacets )
                    {
                        uint tNumNodes = tFacet->number_of_nodes() ;
                        for( uint k=0; k< tNumNodes; ++k )
                        {
                            if( tFacet->node( k )->is_flagged() )
                            {
                                mNodeCounter( tFacet->node( k )->index()) = 1;
                            }
                        }
                    }

                    uint tPlusBlockID =  tCut->data( s )( 1 );

                    for( uint k=0; k<mNumberOfNodes; ++k )
                    {
                        if( mNodeCounter( k ) == 1 )
                        {
                            mPlusBlocksTable( k )( tCount( k )++ ) = tPlusBlockID ;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::create_node_duplicates()
        {
            // counter for new nodes
            index_t tCount = 0 ;
            for( Node * tNode : mOriginalNodes )
            {
                // identify the number of duplicated that are to be crated
                uint tN = mLookup( mNumTapesPerNode( tNode->index() ), mNumCutsPerNode( tNode->index() ) );

                // remember the number of duplicates
                mNodeCounter( tNode->index() ) = tN ;

                tCount += tN ;

                // allocate the memory for the duplicate container
                tNode->allocate_duplicate_container( tN );
            }

            mDuplicateNodes.set_size( tCount, nullptr );
            tCount = 0 ;

            Cell< Node * > tWork( 8, nullptr );

            for( Node * tNode : mOriginalNodes )
            {
                // identify the number of duplicated that are to be created
                uint tN = mNodeCounter( tNode->index() ) ;

                for ( uint d = 0; d < tN; ++d )
                {
                    Node * tDuplicate = new Node( ++mMaxNodeID,
                                                  tNode->x(),
                                                  tNode->y(),
                                                  tNode->z());

                    // the linking will be needed for parallelization
                    tNode->add_duplicate( tDuplicate );

                    tWork( d ) = tDuplicate ;

                    // save the duplicate into array
                    mDuplicateNodes( tCount++ ) = tDuplicate;
                }

                tWork( tN  ) = tNode ;

                // link duplicates with original and others
                for ( uint d = 0; d < tN; ++d )
                {
                   Node * tDuplicate = tNode->duplicate(  d );

                    tDuplicate->allocate_duplicate_container( tN );

                    for( uint k=0; k<= tN; ++k )
                    {
                        if( tWork( k )->id() != tDuplicate->id() )
                        {
                            tDuplicate->add_duplicate(  tWork( k ) );
                        }
                    }
                }
            }

            BELFEM_ASSERT( tCount = mDuplicateNodes.size(),
                           "Number of duplicated notes does not match ( is %lu, but expect %lu )",
                           ( long unsigned int ) tCount ,
                           ( long unsigned int ) mDuplicateNodes.size()  );

        }

//------------------------------------------------------------------------------

        void
        Scissors::relink_elements()
        {
            for( Node * tNode : mOriginalNodes )
            {
                if( tNode->number_of_vertices() == 1 )  // simple case
                {
                    id_t tPlus = mPlusBlocksTable( tNode->index() )( 0 );

                    // loop over  all elements connected to node
                    for ( uint e = 0; e < tNode->number_of_elements(); ++e )
                    {
                        // get elements
                        Element * tElement = tNode->element( e );

                        // check if element is positive
                        if( tElement->geometry_tag() == tPlus )
                        {
                            // loop over all nodes
                            for ( uint k = 0; k < tElement->number_of_nodes(); ++k )
                            {
                                // check if this is the right node
                                if ( tElement->node( k )->id() == tNode->id() )
                                {
                                    // replace node
                                    tElement->insert_node( tNode->duplicate( 0 ), k );

                                    // cancel this loop
                                    break;
                                }
                            }

                        }
                    }
                }
                else // complicated case
                {
                    // create duplicate map
                    Map< id_t, uint > tPlus;
                    uint tCount = 0;
                    for ( id_t tID : mPlusBlocksTable( tNode->index() ) )
                    {
                        if ( !tPlus.key_exists( tID ))
                        {
                            tPlus[ tID ] = tCount++;
                        }
                    }

                    uint tNumElems = tNode->number_of_elements();

                    // loop over all elements
                    for ( uint e = 0; e < tNumElems; ++e )
                    {
                        Element * tElement = tNode->element( e );

                        // check if we have to replace the node
                        if ( tPlus.key_exists( tElement->geometry_tag()))
                        {
                            // find node_id
                            uint tNumNodes = tElement->number_of_nodes();
                            for ( uint k = 0; k < tNumNodes; ++k )
                            {
                                if ( tElement->node( k )->id() == tNode->id())
                                {
                                    // replace node
                                    tElement->insert_node(  tNode->duplicate( tPlus( tElement->geometry_tag() ) ), k );

                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::create_connectors()
        {

            ElementFactory tFactory ;

            index_t tSideSetCount = 0 ;
            for ( scissors::CutData * tCutData : mCutsAndTapes )
            {
                if ( !tCutData->is_tape() )
                {
                    tSideSetCount += tCutData->num_sidesets() ;
                }
            }
            mConnectorSetIDs.set_size( tSideSetCount, 0 );
            tSideSetCount = 0 ;

            for ( scissors::CutData * tCutData : mCutsAndTapes )
            {
                if( ! tCutData->is_tape() )
                {
                    index_t tNumSideSets = tCutData->num_sidesets();

                    // count connectors
                    this->unflag_nodes() ;

                    // count connectors for this cut
                    index_t tCount = 0 ;

                    for ( index_t s = 0; s < tNumSideSets; ++s )
                    {
                        // get the sideset
                        SideSet * tSideSet = mMesh->sideset( tCutData->data( s )( 0 ) ) ;

                        // get the facets of the sideset
                        Cell< Facet * > & tFacets = tSideSet->facets() ;

                        for( Facet * tFacet : tFacets )
                        {
                            uint tNumNodes = tFacet->number_of_nodes() ;
                            for( uint k=0; k<tNumNodes; ++k )
                            {
                                Node * tNode = tFacet->node( k );
                                if( ! tNode->is_flagged() )
                                {
                                    tNode->flag() ;
                                    ++tCount ;
                                }
                            }
                        }
                    }

                    // create new cut
                    SideSet * tCut = new SideSet( tCutData->id(), tCount );

                    //SideSet * tCut = new SideSet( 10, tCount );
                    Cell< Facet * > & tConnectors = tCut->facets() ;

                    // reset counter
                    tCount = 0 ;

                    for ( index_t s = 0; s < tNumSideSets; ++s )
                    {

                        // get the sideset
                        SideSet * tSideSet = mMesh->sideset( tCutData->data( s )( 0 ) ) ;

                        // get the facets of the sideset
                        Cell< Facet * > & tFacets = tSideSet->facets() ;

                        // get id of plus block
                        id_t tPlus = tCutData->data( s )( 1 ) ;

                        uint tNumNodesPerFacet = mesh::number_of_nodes( tSideSet->element_type() );
                        Cell< Node * > tMasterNodes( tNumNodesPerFacet, nullptr );
                        Cell< Node * > tSlaveNodes( tNumNodesPerFacet, nullptr );
                        Cell< Node * > tSlaveNodesOrg( tNumNodesPerFacet, nullptr );
                        Vector< uint > tWork( tNumNodesPerFacet, 0 );

                        for( Facet * tFacet : tFacets )
                        {
                            tFacet->master()->get_nodes_of_facet( tFacet->master_index(), tMasterNodes );
                            tFacet->slave()->get_nodes_of_facet( tFacet->slave_index(), tSlaveNodesOrg );
                            this->rearrange_nodes( tMasterNodes, tSlaveNodesOrg, tSlaveNodes, tWork );
                            for( uint k=0; k<tNumNodesPerFacet; ++k )
                            {
                                Node * tNode = tFacet->node( k );
                                if( tNode->is_flagged() )
                                {
                                    // create a new element
                                    Element * tElement = tFactory.create_lagrange_element( ElementType::LINE2, ++mMaxElementID );

                                    if( tFacet->master()->geometry_tag() == tPlus )
                                    {
                                        tElement->insert_node( tSlaveNodes( k ), 0 );
                                        tElement->insert_node( tMasterNodes( k ), 1 );
                                    }
                                    else if ( tFacet->slave()->geometry_tag() == tPlus )
                                    {
                                        tElement->insert_node( tMasterNodes( k ), 0 );
                                        tElement->insert_node( tSlaveNodes( k ), 1 );
                                    }
                                    else
                                    {
                                        BELFEM_ERROR( false, "Error while creating cut on Sideset %lu", ( long unsigned int ) tSideSet->id() );
                                    }

                                    // we need this information for the dof linking (todo, can be done prettier)
                                    tElement->set_block_id( tPlus );

                                    // add new facet to  sideset
                                    tConnectors( tCount++ ) = new Facet( tElement );
                                    tNode->unflag() ;
                                }
                            }
                        }

                        mSideSetToCutIDs[ tSideSet->id() ] = tCut->id() ;

                        // remember ID of new block
                        mConnectorSetIDs( tSideSetCount++ ) = tCut->id() ;

                        // add cut to mesh
                        mMesh->cuts().push( tCut );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::relink_facets()
        {
            // identify facets that need to be relinked
            mMesh->unflag_all_nodes() ;
            this->flag_nodes() ;

            for( mesh::SideSet * tSideSet : mMesh->sidesets() )
            {
                Cell< Facet * > & tFacets = tSideSet->facets() ;
                uint tNumNodes = mesh::number_of_nodes( tSideSet->element_type() );
                Cell< Node * > tNodes( tNumNodes, nullptr );

                for( Facet * tFacet : tFacets )
                {
                    bool tFlag = false ;
                    for( uint k=0; k<tNumNodes; ++k )
                    {
                        if( tFacet->node( k )->is_flagged() )
                        {
                            tFlag = true ;
                            break ;
                        }
                    }
                    if( tFlag )
                    {
                        // grab nodes from master
                        tFacet->master()->get_nodes_of_facet( tFacet->master_index(), tNodes );

                        // relink the nodes on the element
                        for( uint k=0; k<tNumNodes; ++k )
                        {
                            tFacet->element()->insert_node( tNodes( k ), k );
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        uint
        Scissors::inside_node( const ElementType aElementType )
        {
            switch( aElementType )
            {
                case( ElementType::LINE2 ) :
                case( ElementType::LINE3 ) :
                case( ElementType::LINE4 ) :
                case( ElementType::LINE5 ) :
                {
                    return 2 ;
                }
                case( ElementType::TRI3 ) :
                {
                    return 3 ;
                }
                case( ElementType::TRI6 ) :
                {
                    return 6 ;
                }
                case( ElementType::TRI10 ) :
                {
                    return 9 ;
                }
                case( ElementType::TRI15 ) :
                {
                    return 12 ;
                }
                case( ElementType::QUAD8 ) :
                case( ElementType::QUAD9 ) :
                {
                    return 8 ;
                }
                case( ElementType::QUAD16 ) :
                {
                    return 12 ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Element Type" );
                    return  BELFEM_UINT_MAX ;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        Scissors::unflag_nodes()
        {
            for( Node * tNode : mOriginalNodes )
            {
                tNode->unflag() ;
            }
            for( Node * tNode : mDuplicateNodes )
            {
                tNode->unflag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::flag_nodes()
        {
            for( Node * tNode : mOriginalNodes )
            {
                tNode->flag() ;
            }
            for( Node * tNode : mDuplicateNodes )
            {
                tNode->flag() ;
            }
        }

//----------------------------------------------------------------------------

        void
        Scissors::add_nodes_to_mesh()
        {
            Cell< Node * > tOldNodes ;
            tOldNodes.vector_data() = std::move( mNodes.vector_data() );
            mNodes.set_size( tOldNodes.size() + mDuplicateNodes.size(), nullptr );

            uint tCount = 0 ;
            for( Node * tNode : tOldNodes )
            {
                mNodes( tCount++ ) = tNode ;
            }
            for( Node * tNode : mDuplicateNodes )
            {
                mNodes( tCount++ ) = tNode ;
            }
        }

//------------------------------------------------------------------------------

        void
        Scissors::rearrange_nodes(
                Cell< Node * > & aMaster,
                Cell< Node * > & aSlaveIn,
                Cell< Node * > & aSlaveOut,
                Vector< uint > & aFlags )
        {
            // backup flags and flag all nodes
            uint tCount = 0 ;
            for( Node * tNode : aSlaveIn )
            {
                aFlags( tCount++ ) = tNode->is_flagged() ? 1 : 0 ;
                tNode->flag() ;
            }
            tCount = 0 ;
            for( Node * tMaster : aMaster )
            {
                for( Node * tSlave : aSlaveIn )
                {
                    if( tSlave->is_flagged() )
                    {
                        if( tSlave->x() == tMaster->x() )
                        {
                            if( tSlave->y() == tMaster->y() )
                            {
                                if( tSlave->z() == tMaster->z() )
                                {
                                    aSlaveOut( tCount++ ) = tSlave ;
                                    tSlave->unflag() ;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // restore flags
            tCount = 0 ;
            for( Node * tNode : aSlaveIn )
            {
                aFlags( tCount++ ) == 1 ? tNode->flag() : tNode->unflag() ;
            }
        }

//------------------------------------------------------------------------------
    }
}