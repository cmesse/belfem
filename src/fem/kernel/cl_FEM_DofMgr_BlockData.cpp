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
#include "cl_FEM_DofMgr_BlockData.hpp"

#include "cl_DynamicBitset.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//------------------------------------------------------------------------------

            BlockData::BlockData(  DofManager * aParent ) :
                mParent( aParent ),
                mKernel( aParent->parent() ),
                mMesh( aParent->parent()->mesh() ),
                mCommRank( comm_rank() ),
                mCommSize( comm_size() )
            {

            }

//------------------------------------------------------------------------------

            BlockData::~BlockData()
            {
               this->reset() ;
            }

//------------------------------------------------------------------------------

            void
            BlockData::reset()
            {
                if( mEmptyBlock != nullptr )
                {
                    delete mEmptyBlock;
                    mEmptyBlock = nullptr ;
                }


                for( Block * tBlock : mBlocks )
                {
                    delete tBlock ;
                }

                mBlocks.clear() ;
                mBlockMap.clear() ;
            }

//------------------------------------------------------------------------------

            void
            BlockData::create_blocks()
            {
                // restore factory settings
                this->reset() ;

                // for good measure, we also reset element indices
                mMesh->update_element_indices() ;
                mMesh->unflag_all_elements() ;
                mMesh->unflag_all_facets() ;

                // check which blocks have been selected
                const Vector< id_t > & tBlockIDs = mParent->iwg()->selected_blocks() ;

                // We use key_exists to check if a block is selected
                index_t tCount = 0 ;
                Map< id_t, index_t > tBlockIndices ;
                for ( id_t tID : tBlockIDs )
                {
                    if ( mMesh->block_exists( tID ) )
                    {
                        tBlockIndices[ tID ] = tCount++ ;
                    }
                }

                Vector< id_t> tAllSelectedThinShellFacets ;
                this->collect_thin_shell_facet_ids( tBlockIndices, tAllSelectedThinShellFacets );


                const Vector< id_t > & tSideSetIDs = mParent->iwg()->selected_sidesets() ;

                Map< id_t, index_t > tSideSetIndices ;
                tCount = 0 ;
                for ( id_t tID : tSideSetIDs )
                {
                    if ( mMesh->sideset_exists( tID ) )
                    {
                        if (   mMesh->sideset( tID )->domain_type() == DomainType::ThinShell
                            || mMesh->sideset( tID )->domain_type() == DomainType::Inactive ) continue;
                        tSideSetIndices[ tID ] = tCount++ ;
                    }
                }

                Cell< Vector< index_t > > tOwnedElementIndices, tAuraElementIndices ;
                this->collect_element_indices(
                    tBlockIndices,
                    tSideSetIndices,
                    tAllSelectedThinShellFacets,
                    tOwnedElementIndices,
                    tAuraElementIndices );

                for ( id_t b : tBlockIDs )
                {
                    // skip of block doesn't exist on this proc
                    if ( ! tBlockIndices.key_exists( b  ) ) continue ;

                    const Vector< index_t > & tOwnedIndices = tOwnedElementIndices( tBlockIndices( b ) ) ;
                    const Vector< index_t > & tAuraIndices = tAuraElementIndices( tBlockIndices( b ) ) ;

                    if ( tOwnedIndices.length() > 0 || tAuraIndices.length() > 0 )
                    {
                         Block * tBlock = new Block(
                            mParent,
                            mMesh->block( b ),
                            tOwnedIndices,
                            tAuraIndices );
                        mBlocks.push( tBlock );
                        mBlockMap[ b ] = tBlock ;
                    }
                }

                // empty block class. Might not be necessary
                mEmptyBlock = new Block( mParent ) ;

                if ( tAllSelectedThinShellFacets.length() > 0 )
                {
                    this->link_thin_shell_facets() ;
                }

                // might not be needed but can't hurt
                comm_barrier();
            }

            void
            BlockData::collect_thin_shell_facet_ids(
                const Map< id_t, index_t > & aBlockIndexMap,
                Vector< id_t > & aAllSelectedThinShellFacets )
            {
                if ( mCommRank == 0 )
                {
                    // counter for selected thin shell facets
                    index_t tCount = 0 ;
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        for ( mesh::Block * tBlock : tShell->blocks() )
                        {
                            if ( aBlockIndexMap.key_exists( tBlock->id() ) )
                            {
                                tCount += tShell->facets().size() ;
                                break;
                            }
                        }
                    }

                    // now we collect the IDs for all thin shell facets
                    aAllSelectedThinShellFacets.set_size( tCount );
                    tCount = 0 ;
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        for ( mesh::Facet * tFacet : tShell->facets() )
                        {
                            aAllSelectedThinShellFacets( tCount++ ) = tFacet->id() ;
                        }
                    }

                    comm_barrier() ;
                    share( aAllSelectedThinShellFacets );
                }
                else
                {
                    comm_barrier() ;
                    receive( aAllSelectedThinShellFacets );
                }
            }


             void
             BlockData::collect_element_indices(
                   const Map< id_t, index_t > & aBlockIndexMap,
                   const Map< id_t, index_t > & aSideSetIndexMap,
                   const Vector< id_t >       & aAllSelectedThinShellFacets,
                   Cell< Vector< index_t > >  & aOwnedElementIndices,
                   Cell< Vector< index_t > >  & aAuraElementIndices )
            {
                DynamicBitset tElementBitset( mMesh->number_of_elements() );
                DynamicBitset tFacetBitset( mMesh->number_of_facets() );
                Cell< index_t > tElementIndices ;
                Cell< index_t > tFacetIndices ;
                for ( id_t f : aAllSelectedThinShellFacets )
                {
                    // check if the facet exists on the current mesh
                    if ( mMesh->facet_exists( f ) )
                    {
                        mesh::Facet * tFacet = mMesh->facet( f ) ;

                        // check if we own it
                        if ( tFacet->owner() == mCommRank )
                        {
                            tFacetBitset.set( tFacet->index() );
                        }
                        else // flag for aura processing
                        {
                            tFacet->flag() ;
                        }
                    }
                }


                // first we select the elements that we own
                for ( auto tPair : aBlockIndexMap )
                {
                    // grab the block
                    mesh::Block * tBlock = mMesh->block( tPair.first ) ;

                    for ( mesh::Element * tElement : tBlock->elements() )
                    {
                        if ( tElement->owner() == mCommRank )
                        {
                            tElementBitset.set( tElement->index() );
                        }
                        else // tag this element for aura checking
                        {
                            tElement->flag() ;
                        }
                    }
                }

                // next, we check the elements that belong to the facets
                for ( auto tPair : aSideSetIndexMap )
                {
                    // grab the sideset
                    mesh::SideSet * tSideSet = mMesh->sideset( tPair.first ) ;

                    for ( mesh::Facet * tFacet : tSideSet->facets() )
                    {
                        if ( tFacet->owner() == mCommRank )
                        {
                            if ( tFacet->has_master() )
                            {
                                mesh::Element * tElement = tFacet->master() ;
                                if ( ! aBlockIndexMap.key_exists( tElement->block_id() ) ) continue ;
                                if ( tElement->owner() == mCommRank )
                                {
                                    tElementBitset.set( tElement->index() );
                                }
                                else // tag this element for aura checking
                                {
                                    tElement->flag() ;
                                }
                            }
                            if ( tFacet->has_slave() )
                            {
                                mesh::Element * tElement = tFacet->slave() ;
                                if ( ! aBlockIndexMap.key_exists( tElement->block_id() ) ) continue ;
                                if ( tElement->owner() == mCommRank )
                                {
                                    tElementBitset.set( tElement->index() );
                                }
                                else // tag this element for aura checking
                                {
                                    tElement->flag() ;
                                }
                            }
                        }
                    }
                }

                // this logic selects the aura elements. the aura must span the
                // full node disc of the owned elements, not just their facet
                // neighbors: the patch recovery in the postprocessors needs
                // every element that touches an owned node, including
                // vertex-only neighbors ( the distributor ghosts the full
                // disc, so all flagged candidates exist on this mesh )
                tElementBitset.where( tElementIndices );
                tElementBitset.reset();

                Cell< mesh::Element * > & tElements = mMesh->elements() ;
                for ( index_t e : tElementIndices )
                {
                    mesh::Element * tElement = tElements(e);

                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        mesh::Node * tNode = tElement->node( k )->original() ;

                        for ( uint i=0; i<tNode->number_of_elements(); ++i )
                        {
                            if ( tNode->element( i )->is_flagged() )
                            {
                                tElementBitset.set( tNode->element( i )->index() );
                            }
                        }
                        for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                        {
                            mesh::Node * tDup = tNode->duplicate( d );
                            for ( uint i=0; i<tDup->number_of_elements(); ++i )
                            {
                                if ( tDup->element( i )->is_flagged() )
                                {
                                    tElementBitset.set( tDup->element( i )->index() );
                                }
                            }
                        }
                    }
                }

                // the walk above spans the discs of the owned elements. a proc
                // can also own a node without owning any element around it
                // ( reentrant corners of the partition boundary ). the recovery
                // patches of such nodes need the full disc as well, so we also
                // expand from the owned nodes ( the distributor ghosts their
                // complete discs, so all flagged candidates exist locally )
                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    if ( tNode->owner() != mCommRank ) continue ;

                    mesh::Node * tOrg = tNode->original() ;

                    for ( uint i=0; i<tOrg->number_of_elements(); ++i )
                    {
                        if ( tOrg->element( i )->is_flagged() )
                        {
                            tElementBitset.set( tOrg->element( i )->index() );
                        }
                    }
                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        mesh::Node * tDup = tOrg->duplicate( d );
                        for ( uint i=0; i<tDup->number_of_elements(); ++i )
                        {
                            if ( tDup->element( i )->is_flagged() )
                            {
                                tElementBitset.set( tDup->element( i )->index() );
                            }
                        }
                    }
                }

                // we also tag the facets that should belong to the aura
                tFacetBitset.where( tFacetIndices );
                tFacetBitset.reset();
                Cell< mesh::Facet * > & tFacets = mMesh->facets() ;

                for ( index_t f : tFacetIndices )
                {
                    mesh::Facet * tFacet = tFacets( f ) ;

                    for ( uint k=0; k<tFacet->number_of_facets(); ++k )
                    {
                        mesh::Facet * tOtherFacet = tFacet->facet( k ) ;
                        if ( tOtherFacet->is_flagged() ) // check if this facet sits on the aura
                        {
                            // only elements of SELECTED blocks may enter the
                            // aura: kernels that do not select all blocks
                            // ( e.g. the thermal kernel excludes air ) have no
                            // container for foreign elements, and the counting
                            // loop below would die on the block-index lookup
                            if ( tOtherFacet->has_master() )
                            {
                                mesh::Element * tElement = tOtherFacet->master() ;
                                if ( tElement->owner() != mCommRank // check if master should be on aura
                                     && aBlockIndexMap.key_exists( tElement->block_id() ) )
                                {
                                    tElementBitset.set( tElement->index() );
                                }
                            }
                            if ( tOtherFacet->has_slave() )
                            {
                                mesh::Element * tElement = tOtherFacet->slave() ;
                                if ( tElement->owner() != mCommRank // check if facet should be on aura
                                     && aBlockIndexMap.key_exists( tElement->block_id() ) )
                                {
                                    tElementBitset.set( tElement->index() );
                                }
                            }
                        }
                    }
                }
                uint tNumBlocks = aBlockIndexMap.size();

                aOwnedElementIndices.set_size( tNumBlocks, {} );
                aAuraElementIndices.set_size( tNumBlocks, {} );


                // count elements
                Vector< index_t > tElementCounters( tNumBlocks );

                for ( uint k=0; k<2; ++k )
                {
                    Cell< Vector< index_t > > & tIndices = ( k==0 ) ? aOwnedElementIndices : aAuraElementIndices ;

                    // reset counters
                    tElementCounters.fill( 0 );

                    // count elements
                    for ( index_t e : tElementIndices )
                    {
                        ++tElementCounters( aBlockIndexMap( tElements( e )->block_id() ) );
                    }

                    // allocate memory
                    for ( uint b=0; b<tNumBlocks; ++b )
                    {
                        tIndices( b ).set_size( tElementCounters( b ) );
                    }

                    // reset counters
                    tElementCounters.fill( 0 );

                    // populate indices
                    for ( index_t e : tElementIndices )
                    {
                        // get block index
                        uint b = aBlockIndexMap( tElements( e )->block_id() ) ;

                        // add element index to table
                        tIndices( b )( tElementCounters( b )++ ) = e ;
                    }

                    if ( mCommSize  < 2 || k == 1 ) break ;

                   tElementBitset.where( tElementIndices );
                }

            }

//------------------------------------------------------------------------------

            void
            BlockData::link_thin_shell_facets()
            {
                index_t tSize = 0 ;
                if ( mCommRank == 0 ) tSize = mMesh->thin_shells().size() ;

                comm_barrier() ;
                broadcast( tSize );

                if ( tSize == 0 ) return ;

                if ( comm_size() < 2 )
                {
                    this->link_thin_shell_facets_serial() ;
                }
                else
                {
                    this->link_thin_shell_facets_parallel() ;
                }
            }

            void
            BlockData::link_thin_shell_facets_serial()
            {
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    // grab the facets
                    Cell< mesh::Facet * > & tFacets = tShell->facets();

                    // loop over all blocks
                    for ( mesh::Block * tMeshBlock : tShell->blocks() )
                    {
                        // get the fem block
                        Block * tBlock = mBlockMap( tMeshBlock->id() ) ;

                        index_t tCount = 0 ;
                        for ( mesh::Element * tElement : tMeshBlock->elements() )
                        {
                            tBlock->element( tElement->id() )->set_facet( tFacets( tCount++ ) );
                        }
                    }

                    // we don't have the connectors in the thermal problem
                    if ( mParent->iwg()->type() != IwgType::Maxwell ) continue ;

                    for ( mesh::Block * tMeshBlock : tShell->side_connector_blocks() )
                    {
                        // get the fem block
                        Block * tBlock = mBlockMap( tMeshBlock->id() ) ;

                        for ( mesh::Element * tElement : tMeshBlock->elements() )
                        {
                            // the recovery facet sits at wall id + 1 by construction
                            mesh::Facet * tFacet = mMesh->facet( tElement->id() + 1 );

                            // setup runs once per run: always-active check,
                            // the deref below is unconditional
                            BELFEM_ERROR( tFacet->master() != nullptr,
                                "recovery facet %lu has no master",
                                ( long unsigned int ) tFacet->id() );

                            Element * tWall = tBlock->element( tElement->id() );

                            // the kernel reads element()->facet() and facet()->index_on_master()
                            tWall->set_facet( tFacet );

                            // the facet's master is the layer block element the wall spans
                            tWall->set_reference(
                                mBlockMap( tFacet->master()->block_id() )
                                    ->element( tFacet->master()->id() ) );
                        }
                    }
                }
            }

            void
            BlockData::link_thin_shell_facets_parallel()
            {
                Vector< id_t > tData ;
                Vector< id_t > tSideData ;

                if ( mCommRank == 0 )
                {
                    // all element-facet pairs go to everyone: besides its owner,
                    // an element may sit on the aura of other procs, and those
                    // copies need the facet link as well
                    index_t tCount = 0 ;
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        // loop over all blocks
                        for ( mesh::Block * tMeshBlock : tShell->blocks() )
                        {
                            tCount += tMeshBlock->elements().size() ;
                        }
                    }

                    tData.set_size( 2 * tCount );
                    tCount = 0 ;

                    // populate data
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        // grab the facets
                        Cell< mesh::Facet * > & tFacets = tShell->facets();

                        // loop over all blocks
                        for ( mesh::Block * tBlock : tShell->blocks() )
                        {
                            index_t k = 0 ;
                            for ( mesh::Element * tElement : tBlock->elements() )
                            {
                                // write element id
                                tData( tCount++ ) = tElement->id();

                                // write facet id
                                tData( tCount++ ) = tFacets( k++ )->id();
                            }
                        }
                    }

                    // now the side connectors (if they exist)

                    tCount = 0 ;
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        // loop over all blocks
                        for ( mesh::Block * tMeshBlock : tShell->side_connector_blocks() )
                        {
                            tCount += tMeshBlock->elements().size() ;
                        }
                    }
                    tSideData.set_size( tCount );
                    tCount = 0 ;
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        // loop over all blocks
                        for ( mesh::Block * tMeshBlock : tShell->side_connector_blocks() )
                        {
                            for ( mesh::Element * tMeshElement : tMeshBlock->elements() )
                            {
                                tSideData( tCount++ ) = tMeshElement->id() ;
                            }
                        }
                    }

                    comm_barrier() ;
                    share( tData );
                    share( tSideData );
                }
                else
                {
                    comm_barrier() ;
                    receive( tData );
                    receive( tSideData );
                }

                index_t tNumElements = tData.length() / 2 ;
                index_t tCount = 0 ;

                for ( index_t e=0; e<tNumElements; ++e )
                {
                    id_t tElementID = tData( tCount++ );
                    id_t tFacetID   = tData( tCount++ );

                    // element or facet may not be part of this submesh
                    if ( ! mMesh->element_exists( tElementID ) ) continue ;
                    if ( ! mMesh->facet_exists( tFacetID ) ) continue ;

                    mesh::Element * tElement = mMesh->element( tElementID );

                    if ( ! mBlockMap.key_exists( tElement->block_id() ) ) continue ;

                    Block * tBlock = mBlockMap( tElement->block_id() );

                    BELFEM_ASSERT( tElement->owner() != mCommRank
                        || tBlock->element_exists( tElementID ),
                        "FEM Element %lu on block %lu was not created.",
                        ( long unsigned int ) tElementID,
                        ( long unsigned int ) tElement->block_id() );

                    // an aura copy exists only if the element touches this proc
                    if ( ! tBlock->element_exists( tElementID ) ) continue ;

                    tBlock->element( tElementID )->set_facet( mMesh->facet( tFacetID ) );
                }

                for ( id_t e : tSideData )
                {
                    if ( ! mMesh->element_exists( e ) ) continue ;

                    mesh::Element * tElement = mMesh->element( e );

                    if ( ! mBlockMap.key_exists( tElement->block_id() ) ) continue ;
                    Block * tBlock = mBlockMap( tElement->block_id() );

                    BELFEM_ASSERT( tElement->owner() != mCommRank
                       || tBlock->element_exists( e ),
                       "FEM Element %lu on block %lu was not created.",
                       ( long unsigned int ) e,
                       ( long unsigned int ) tElement->block_id() );

                    // an aura copy exists only if the element touches this proc
                    if ( ! tBlock->element_exists( e ) ) continue ;

                    // owned wall + missing recovery facet = distribution bug.
                    // setup runs once per run: always-active check - a silent
                    // skip would leave the wall without a facet and crash
                    // later in assembly
                    BELFEM_ERROR( tElement->owner() != mCommRank
                        || mMesh->facet_exists( e + 1 ),
                        "recovery facet %lu missing for owned side connector %lu",
                        ( long unsigned int )( e + 1 ),
                        ( long unsigned int ) e );

                    // on a rank that holds the wall only as an aura copy the
                    // facet may not have shipped - skip, assembly only runs
                    // owned elements
                    if ( ! mMesh->facet_exists( e + 1 ) ) continue ;

                    // the recovery facet sits at wall id + 1 by construction
                    mesh::Facet * tFacet = mMesh->facet( e + 1 );

                    Element * tWall = tBlock->element( e );

                    // the kernel reads element()->facet() and facet()->index_on_master()
                    tWall->set_facet( tFacet );

                    // the facet's master is the layer block element the wall
                    // spans; always-active check before the deref below
                    BELFEM_ERROR( tFacet->master() != nullptr,
                        "recovery facet %lu has no master",
                        ( long unsigned int ) tFacet->id() );

                    mesh::Element * tMaster = tFacet->master() ;

                    if ( ! mBlockMap.key_exists( tMaster->block_id() ) ) continue ;
                    if ( ! mBlockMap( tMaster->block_id() )->element_exists( tMaster->id() ) ) continue ;

                    tWall->set_reference(
                        mBlockMap( tMaster->block_id() )->element( tMaster->id() ) );
                }


            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */