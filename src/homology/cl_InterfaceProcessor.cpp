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

#include "cl_InterfaceProcessor.hpp"
#include "fn_unique.hpp"

namespace belfem
{
    namespace mesh
    {

        InterfaceSet::InterfaceSet( Mesh * aMesh, const string aBitsetHex ) :
            mMesh( aMesh ),
            mBitset( new DynamicBitset( aMesh->number_of_blocks() ) )
        {
            mBitset->set_from_hex( aBitsetHex );

            // a coil on either side decouples the interface from the solve
            for ( uint b=0; b<aMesh->number_of_blocks(); ++b )
            {
                if ( mBitset->test( b ) and aMesh->blocks()( b )->domain_type() == DomainType::Coil )
                {
                    mTreatment = InterfaceTreatment::Decouple ;
                    break ;
                }
            }
        }

        InterfaceSet::~InterfaceSet()
        {
            delete mBitset;
        }

        void
        InterfaceSet::increment_sideset_counter()
        {
            ++mNumberOfSideSets ;
        }

        void
        InterfaceSet::allocate_sideset_container()
        {
            mSideSets.set_size( mNumberOfSideSets, nullptr );
            mNumberOfSideSets = 0;
        }

        void
        InterfaceSet::add_sideset( SideSet * aSideSet )
        {
            mSideSets( mNumberOfSideSets++ ) = aSideSet ;
        }

        void
        InterfaceSet::collect_nodes_and_elements()

        {
            Cell< Node * > tNodes ;
            for ( SideSet * tSideSet : mSideSets )
            {
                for ( Facet * tFacet : tSideSet->facets() )
                {
                    // note that we can't use tFacet->flag_nodes() here
                    // because the mesh is not finalized
                    tFacet->master()->get_nodes_of_facet( tFacet->index_on_master(), tNodes ) ;
                    for ( Node * tNode : tNodes )
                    {
                        // get the original of the node (might be the node itself)
                        Node * tOrg = tNode->original();
                        tOrg->flag();
                        for ( uint k=0; k<tOrg->number_of_duplicates(); ++k )
                        {
                            tOrg->duplicate( k )->flag();
                        }
                    }
                }
            }

            // collect nodes (note that for cuts, we also flag duplicates at the surface)
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    mOriginals[ tNode->id() ]= tNode ;
                }
            }

            // collect the blocks
            Vector< id_t > tBlockIDs( mSideSets.size(), 0 );
            index_t tCount = 0 ;
            for ( SideSet * tSideSet : mSideSets )
            {
                // get first facet
                Facet * tFacet = tSideSet->facets()( 0 );

                // add blockid of master block to container
                // ( we must duplicate  the masters since we don't want nodes on the air side hanging
                //     unless for the cuts )
                tBlockIDs( tCount++ ) = tFacet->master()->block_id() ;
            }

            unique( tBlockIDs );

            // count the elements that need to be relinked
            tCount = 0 ;
            for ( id_t tBlockID : tBlockIDs )
            {
                Cell< Element * > & tElements = mMesh->block( tBlockID )->elements() ;
                for ( Element * tElement : tElements )
                {
                    // check if any node of this element is flagged
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        if ( tElement->node( k )->is_flagged() )
                        {
                            tElement->flag();
                            ++tCount ;
                            break ;
                        }
                    }
                }
            }

            // collect elements
            mElements.set_size( tCount, nullptr ) ;
            tCount = 0 ;
            for ( id_t tBlockID : tBlockIDs )
            {
                Cell< Element * > & tElements = mMesh->block( tBlockID )->elements() ;
                for ( Element * tElement : tElements )
                {
                    if ( tElement->is_flagged() )
                    {
                        mElements( tCount++ ) = tElement;
                        tElement->unflag() ;
                    }
                }
            }

            // finally, we unflag these nodes for now
            for ( auto tPair : mOriginals )
            {
                tPair.second->unflag();
            }

            // now all nodes and all elements should be unflagged
        }

        void
        InterfaceSet::duplicate_nodes( id_t & aMaxNodeID )
        {
            Cell< Basis * > tSources ;
            Vector< real >  tWeights ;

            for ( auto tPair : mOriginals )
            {
                Node * tOriginal = tPair.second ;
                Node * tDuplicate = new Node( ++aMaxNodeID, tOriginal->x(), tOriginal->y(), tOriginal->z() ) ;
                if ( tOriginal->is_hanging() )
                {
                    uint n = tOriginal->number_of_sources();

                    tSources.set_size( n, nullptr ) ;
                    tWeights.set_size( n, 0.0 );
                    for ( uint k=0; k<n; ++k )
                    {
                        tSources( k ) = tOriginal->source( k ) ;
                        tWeights( k ) = tOriginal->weight( k ) ;
                    }

                    tDuplicate->set_sources( tSources, tWeights );
                }
                mDuplicates[ tOriginal->id() ] = tDuplicate ;
            }

            if ( mMesh->has_periodicity() )
            {
                Periodicity * tPeriodicity = mMesh->periodicity();

                for ( auto tPair : mOriginals )
                {
                    Node * A = tPair.second ;

                    if ( ! A->is_periodic() )
                    {
                        A->flag( 6 );
                    }
                    else
                    {
                        A->unflag( 6 );
                    }
                }
                for ( auto tPair : mOriginals )
                {
                    Node * A = tPair.second ;
                    if ( A->is_flagged( 6 ) ) continue;
                    Node * B = A->periodic();

                    // the periodic map may cross block pairs ( e.g. twisted helix:
                    // front of conductor k maps to back of conductor k-1 ), so the
                    // partner can live in a sibling interface set. those pairs are
                    // tied by InterfaceProcessor::pair_cross_set_periodic_duplicates()
                    if ( ! mDuplicates.key_exists( B->id() ) ) continue;

                    Node * C = mDuplicates( A->id() );
                    Node * D = mDuplicates( B->id() );

                    A->flag( 6 );
                    B->flag( 6 );

                    C->set_periodic( D );
                    D->set_periodic( C );

                    tPeriodicity->add_node_pair_to_backup( C, D );
                }
                for ( auto tPair : mOriginals )
                {
                    tPair.second->unflag( 6 );
                }
            }
        }

        void
        InterfaceSet::relink_elements()
        {
            // select all original nodes
            for ( auto tPair : mOriginals )
            {
                tPair.second->flag();
            }

            for ( Element * tElement : mElements )
            {
                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    if ( tElement->node( k )->is_flagged() )
                    {
                        tElement->insert_node( mDuplicates( tElement->node( k )->id() ), k );
                    }
                }
            }

            for ( auto tPair : mOriginals )
            {
                tPair.second->unflag();
            }
        }


        index_t
        InterfaceSet::number_of_nodes() const
        {
            return mOriginals.size();
        }

        void
        InterfaceSet::add_duplicates( index_t & aCount, Cell< Node * > & aDuplicates, Map< id_t, Node * > & aOriginalMap )
        {
            for ( auto tPair : mDuplicates )
            {
                aDuplicates( aCount++ ) = tPair.second;

                // decoupled duplicates stay unsourced and unregistered
                if ( mTreatment != InterfaceTreatment::Decouple )
                {
                    aOriginalMap[ tPair.second->id() ] = mOriginals( tPair.first );
                }
            }
        }

        InterfaceProcessor::InterfaceProcessor(
                Mesh * aMesh,
                Topology * aTopology,
                Cell< Node * > & aAbstractNodes,
                const uint aNumOriginalSideSets,
                const id_t aMaxNodeID ):
            mMesh( aMesh ),
            mTopology( aTopology ),
            mAbstractNodes( aAbstractNodes ),
            mNumOriginalSidesets( aNumOriginalSideSets ),
            mMaxNodeID( aMaxNodeID )
        {
            // create the block indices
            index_t tCount = 0 ;
            for ( Block * tBlock : aMesh->blocks() )
            {
                mBlockIndices[ tBlock->id() ] = tCount++ ;
            }

            // create the bitset
            mIsAirBlock = new DynamicBitset( aMesh->number_of_blocks() );
            mIsFerroBlock = new DynamicBitset( aMesh->number_of_blocks() );

            // populate the bitset
            for ( id_t tID : mTopology->groups( DomainType::Air ) )
            {
                mIsAirBlock->set( mBlockIndices( tID ) );
            }
            for ( id_t tID : mTopology->groups( DomainType::Ferro ) )
            {
                mIsFerroBlock->set( mBlockIndices( tID ) );
            }

            this->create_interface_sets();
            this->connect_sidesets_to_interface_sets();
            this->populate_interface_sets();
            this->duplicate_nodes();
            this->relink_elements();
            this->add_duplicate_nodes_to_mesh();
        }

        InterfaceProcessor::~InterfaceProcessor()
        {
            delete mIsAirBlock ;
            delete mIsFerroBlock ;
        }

        void
        InterfaceProcessor::create_interface_sets()
        {
            // count interface classes
            DynamicBitset tBitset( mMesh->number_of_blocks() );

            // collect strings
            Cell< string > tStrings ;

            // temporary map to link sideset ids with strings
            Map< id_t, string > tMapA ;

            // first we need to identify how many interface sets we have
            for ( index_t s=0 ; s < mNumOriginalSidesets ; ++s )
            {
                // get first facet
                Facet * tFacet = mMesh->sidesets()( s )->facets()(0);

                // skip if this is a boundary
                if ( ! tFacet->has_slave() ) continue;

                // get block index of master
                index_t tMasterIndex = mBlockIndices( tFacet->master()->block_id() );

                // get block index of slave
                index_t tSlaveIndex = mBlockIndices( tFacet->slave()->block_id() ) ;

                // check if we want to duplicate at  this interface
                if ( mIsAirBlock->test( tMasterIndex ) and mIsAirBlock->test( tSlaveIndex ) ) continue ;

                // no duplication if neither side is ferro or air
                if ( ! ( ( mIsAirBlock->test( tMasterIndex ) or mIsFerroBlock->test( tMasterIndex ) )
                       or ( mIsAirBlock->test( tSlaveIndex ) or mIsFerroBlock->test( tSlaveIndex ) ) ) ) continue ;

                // now we create the bit string
                tBitset.reset();

                tBitset.set( tMasterIndex );
                tBitset.set( tSlaveIndex );

                // add bitset string to cell
                tStrings.push( tBitset.to_hex() );
                tMapA[  mMesh->sidesets()( s )->id() ] = tBitset.to_hex()  ;
            }

            unique( tStrings );

            Map< string, InterfaceSet * > tMapB ;


            // create the interface sets
            mSets.set_size( tStrings.size(), nullptr );
            uint tCount = 0 ;
            for ( string tString : tStrings )
            {

                InterfaceSet * tSet = new InterfaceSet( mMesh, tString );

                // link set to temporary map
                tMapB[ tString ] = tSet ;

                // add set to container
                mSets( tCount++ ) = tSet ;
            }

            // link sidesets with cutsets
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( tMapA.key_exists( tSideSet->id() ) )
                {
                    mSetsMap[ tSideSet->id() ] = tMapB( tMapA( tSideSet->id() ) );
                }
            }
        }

        void
        InterfaceProcessor::connect_sidesets_to_interface_sets()
        {
            // first we need to allocate the containers
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( mSetsMap.key_exists( tSideSet->id() ) )
                {
                    mSetsMap[ tSideSet->id() ]->increment_sideset_counter() ;
                }
            }

            // allocate the memory for the containers
            for ( InterfaceSet * tSet : mSets )
            {
                tSet->allocate_sideset_container();
            }

            // connect sidesets with interface sets
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( mSetsMap.key_exists( tSideSet->id() ) )
                {
                    mSetsMap[ tSideSet->id() ]->add_sideset( tSideSet );
                }
            }
        }

        void
        InterfaceProcessor::populate_interface_sets()
        {
            mMesh->unflag_all_elements();
            mMesh->unflag_all_nodes() ;
            for ( InterfaceSet * tSet : mSets )
            {
                tSet->collect_nodes_and_elements();
            }
        }

        void
        InterfaceProcessor::duplicate_nodes()
        {
            for ( InterfaceSet * tSet : mSets )
            {
                tSet->duplicate_nodes( mMaxNodeID );
            }

            if ( mMesh->has_periodicity() )
            {
                this->pair_cross_set_periodic_duplicates();
            }
        }

        void
        InterfaceProcessor::pair_cross_set_periodic_duplicates()
        {
            Periodicity * tPeriodicity = mMesh->periodicity();

            for ( InterfaceSet * tSet : mSets )
            {
                Map< id_t, Node * > & tDuplicates = tSet->duplicate_map() ;

                for ( auto tPair : tSet->original_map() )
                {
                    Node * A = tPair.second ;
                    if ( ! A->is_periodic() ) continue;

                    Node * C = tDuplicates( A->id() );

                    // already tied, either within its own set or by this pass
                    // ( from the partner's side )
                    if ( C->is_periodic() ) continue;

                    Node * B = A->periodic();

                    // find the partner's duplicate in the sibling sets
                    Node * D = nullptr ;
                    for ( InterfaceSet * tOther : mSets )
                    {
                        if ( tOther == tSet ) continue;
                        if ( tOther->duplicate_map().key_exists( B->id() ) )
                        {
                            BELFEM_ERROR( D == nullptr,
                                "periodic partner %lu of interface node %lu was duplicated in more than one interface set, pairing is ambiguous",
                                ( long unsigned int ) B->id(),
                                ( long unsigned int ) A->id() );

                            D = tOther->duplicate_map()( B->id() );
                        }
                    }

                    BELFEM_ERROR( D != nullptr,
                        "periodic partner %lu of interface node %lu was not duplicated in any interface set",
                        ( long unsigned int ) B->id(),
                        ( long unsigned int ) A->id() );

                    BELFEM_ERROR( ! D->is_periodic(),
                        "duplicate of periodic partner %lu of interface node %lu is already tied to another node",
                        ( long unsigned int ) B->id(),
                        ( long unsigned int ) A->id() );

                    C->set_periodic( D );
                    D->set_periodic( C );

                    tPeriodicity->add_node_pair_to_backup( C, D );
                }
            }
        }

        void InterfaceProcessor::relink_elements()
        {
            for ( InterfaceSet * tSet : mSets )
            {
                tSet->relink_elements();
            }
        }

        void
        InterfaceProcessor::add_duplicate_nodes_to_mesh()
        {
            // count duplicates
            index_t tCount = 0 ;
            for ( InterfaceSet * tSet : mSets )
            {
                tCount+= tSet->number_of_nodes() ;
            }
            Cell< Node * > tDuplicates( tCount, nullptr );

            // temporaty map with originals
            Map< id_t, Node * > tOriginalMap ;
            // collect duplicates
            tCount = 0 ;
            for ( InterfaceSet * tSet : mSets )
            {
               tSet->add_duplicates( tCount, tDuplicates, tOriginalMap );
            }

            for ( Node * tDup : tDuplicates )
            {
                // decoupled (coil) duplicates have no original -> no sources
                auto tIt = tOriginalMap.find( tDup->id() );
                if ( tIt == tOriginalMap.end() ) continue ;

                // get the original
                Node * tOrg = tIt->second ;

                Cell< Node * > tOneOrg( 1 , nullptr );
                Vector< real > tOne( 1.0, 1 );

                // check if the original is hanging
                if ( ! tOrg->is_hanging() )
                {
                    tOneOrg( 0 ) = tOrg ;
                    tDup->set_sources(  tOneOrg, tOne );
                }
                else if ( tDup->number_of_sources() == 0 )
                {
                    // grab sources and coefficients form original
                    Cell< Basis * > tSources( tOrg->number_of_sources(), nullptr );

                    Vector< real > tWeights( tOrg->number_of_sources() );
                    for ( uint k=0; k<tOrg->number_of_sources(); ++k )
                    {
                        tSources( k ) = tOrg->source( k );
                        tWeights( k ) = tOrg->weight( k );
                    }

                    tDup->set_sources( tSources, tWeights );
                }
            }

            // add the new duplicates to the existing originals
            this->unify_duplicates( tDuplicates );

            // add duplicates to mesh
            append( mMesh->nodes(), tDuplicates );
        }

        void
        InterfaceProcessor::unify_duplicates(  Cell< Node * > & aDuplicates )
        {
            mMesh->unflag_all_nodes();

            // flag all originals
            for ( Node * tNode : aDuplicates )
            {
                for ( uint k=0; k<tNode->number_of_sources(); ++k )
                {
                    tNode->source( k )->flag();
                }
            }

            // exclude abstract nodes
            for ( Node * tNode : mAbstractNodes )
            {
                tNode->unflag();
            }

            // collect all originals
            index_t tCount = 0 ;

            // count originals
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }

            // collect originals
            Cell< Node * > tOriginals( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tOriginals( tCount++ ) = tNode ;
                }
            }

            // count duplicates that need to be added
            Vector< uint > tNumDuplicates( tCount, 0 );

            // in this scenario, we know that the only node that is flagged is our original
            // so we link it accordingly
            for ( Node * tNode : aDuplicates )
            {
                for ( uint k=0; k<tNode->number_of_sources(); ++k )
                {
                    Node * tOrg = reinterpret_cast< Node * >( tNode->source( k ) );
                    if ( tOrg->is_flagged() )
                    {
                        tNode->set_original( tOrg );
                        ++tNumDuplicates( tOrg->index() );
                        break ;
                    }
                }
            }

            Cell< Node * > tDups ;

            // next we need to make sure that we enlarge the containers for the originals
            for ( Node * tNode : tOriginals )
            {
                if ( tNode->number_of_duplicates() > 0 )
                {
                    // save duplicates
                    tDups.set_size( tNode->number_of_duplicates(), nullptr );
                    for ( uint k=0; k<tNode->number_of_duplicates(); ++k )
                    {
                        tDups( k ) = tNode->duplicate( k );
                    }

                    // reset container
                    tNode->reset_duplicate_container();

                    // enlarge the container
                    tNode->allocate_duplicate_container( tDups.size() + tNumDuplicates( tNode->index() ) );

                    // repopulate duplicates
                    for ( Node * tDup : tDups )
                    {
                        tNode->add_duplicate( tDup );
                    }
                }
                else
                {
                    tNode->allocate_duplicate_container( tNumDuplicates( tNode->index() ) );
                }
            }

            // next we add the duplicates to the originals
            // ( decoupled duplicates never got set_original() -> skip them )
            for ( Node * tNode : aDuplicates )
            {
                if ( tNode->is_duplicate() )
                {
                    tNode->original()->add_duplicate( tNode );
                }
            }

        }

    }
}
