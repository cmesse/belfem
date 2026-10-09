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
#include <string>
#include "fn_unique.hpp"

namespace belfem
{
    namespace mesh
    {

        namespace
        {
            const char *
            treatment_label( const InterfaceTreatment aTreatment )
            {
                switch ( aTreatment )
                {
                    case InterfaceTreatment::TieWeight1 :
                    {
                        return "tie-weight" ;
                    }
                    case InterfaceTreatment::Decouple :
                    {
                        return "decoupled" ;
                    }
                }
                BELFEM_ERROR( false, "unknown interface treatment %d", static_cast< int >( aTreatment ) );
                return "" ;
            }

            // the name with the enumerator value, since to_string maps
            // several types to one name ( both coatings print as Coating )
            string
            domain_type_label( const DomainType aType )
            {
                return to_string( aType ) + " ( " + std::to_string( static_cast< int >( aType ) ) + " )" ;
            }

            // the keys of a node map in ascending order: Map is a hash map
            // and its iteration order must not reach the mesh
            void
            sorted_keys( const Map< id_t, Node * > & aMap, Cell< id_t > & aKeys )
            {
                aKeys.set_size( aMap.size(), 0 );
                index_t tCount = 0 ;
                for ( const auto & tPair : aMap )
                {
                    aKeys( tCount++ ) = tPair.first ;
                }
                sort( aKeys );
            }
        }

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
                        // the original may be the node itself
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
                Facet * tFacet = tSideSet->facets()( 0 );

                // the master blocks are duplicated too; otherwise nodes on the
                // air side would be left hanging, except at the cuts
                tBlockIDs( tCount++ ) = tFacet->master()->block_id() ;
            }

            mMasterType = mMesh->block( tBlockIDs( 0 ) )->domain_type() ;

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
        InterfaceSet::duplicate_nodes( id_t & aMaxNodeID, Map< id_t, Node * > & aSharedDuplicates )
        {
            Cell< Basis * > tSources ;
            Vector< real >  tWeights ;

            // copies are numbered in ascending id of their originals
            Cell< id_t > tIDs ;
            sorted_keys( mOriginals, tIDs );

            for ( id_t tID : tIDs )
            {
                if ( aSharedDuplicates.key_exists( tID ) )
                {
                    mDuplicates[ tID ] = aSharedDuplicates( tID );
                    continue ;
                }

                Node * tOriginal = mOriginals( tID );
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
                mDuplicates[ tID ] = tDuplicate ;
                aSharedDuplicates[ tID ] = tDuplicate ;
            }

            if ( mMesh->has_periodicity() )
            {
                Periodicity * tPeriodicity = mMesh->periodicity();

                for ( id_t tID : tIDs )
                {
                    Node * A = mOriginals( tID );
                    if ( ! A->is_periodic() ) continue;
                    Node * B = A->periodic();

                    // the periodic map may cross block pairs ( e.g. twisted helix:
                    // front of conductor k maps to back of conductor k-1 ), so the
                    // partner can live in a sibling interface set. those pairs are
                    // tied by InterfaceProcessor::pair_cross_set_periodic_duplicates()
                    if ( ! mDuplicates.key_exists( B->id() ) ) continue;

                    Node * C = mDuplicates( A->id() );
                    Node * D = mDuplicates( B->id() );

                    // shared copies may have been tied by an earlier set, and
                    // the partner ties the pair again from its side
                    if ( C->is_periodic() || D->is_periodic() )
                    {
                        BELFEM_ERROR( C->periodic() == D && D->periodic() == C,
                            "duplicates of periodic interface nodes %lu and %lu are already tied to other nodes",
                            ( long unsigned int ) A->id(),
                            ( long unsigned int ) B->id() );
                        continue;
                    }

                    C->set_periodic( D );
                    D->set_periodic( C );

                    tPeriodicity->add_node_pair_to_backup( C, D );
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
                if ( ! this->shares_duplicates() )
                {
                    aDuplicates( aCount++ ) = tPair.second;
                }

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
            mMaxNodeID( aMaxNodeID ),
            mIsAirBlock( aMesh->number_of_blocks() ),
            mIsFerroBlock( aMesh->number_of_blocks() )
        {
            // create the block indices
            index_t tCount = 0 ;
            for ( Block * tBlock : aMesh->blocks() )
            {
                mBlockIndices[ tBlock->id() ] = tCount++ ;
            }

            // populate the bitsets
            for ( id_t tID : mTopology->groups( DomainType::Air ) )
            {
                mIsAirBlock.set( mBlockIndices( tID ) );
            }
            for ( id_t tID : mTopology->groups( DomainType::Ferro ) )
            {
                mIsFerroBlock.set( mBlockIndices( tID ) );
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
        }

        void
        InterfaceProcessor::create_interface_sets()
        {
            this->refuse_mixed_master_types();

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
                if ( ! this->duplicates_at( tFacet ) ) continue;

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

        bool
        InterfaceProcessor::duplicates_at( Facet * aFacet ) const
        {
            if ( ! aFacet->has_slave() ) return false ;

            index_t tMasterIndex = mBlockIndices( aFacet->master()->block_id() );
            index_t tSlaveIndex  = mBlockIndices( aFacet->slave()->block_id() );

            // nothing to duplicate inside the air
            if ( mIsAirBlock.test( tMasterIndex ) and mIsAirBlock.test( tSlaveIndex ) ) return false ;

            // nor where neither side is ferro or air
            return ( mIsAirBlock.test( tMasterIndex ) or mIsFerroBlock.test( tMasterIndex ) )
                or ( mIsAirBlock.test( tSlaveIndex )  or mIsFerroBlock.test( tSlaveIndex ) ) ;
        }

        void
        InterfaceProcessor::refuse_mixed_master_types() const
        {
            for ( index_t s = 0; s < mNumOriginalSidesets; ++s )
            {
                SideSet * tSideSet = mMesh->sidesets()( s );

                // a sideset is an interface sideset when one of its facets
                // is duplicated at; every other sideset is left alone
                bool tInterface = false ;
                for ( Facet * tFacet : tSideSet->facets() )
                {
                    if ( this->duplicates_at( tFacet ) )
                    {
                        tInterface = true ;
                        break;
                    }
                }
                if ( ! tInterface ) continue;

                Facet * tFirst = nullptr ;
                DomainType tType = DomainType::UNDEFINED ;

                for ( Facet * tFacet : tSideSet->facets() )
                {
                    // boundary facets take no part
                    if ( ! tFacet->has_slave() ) continue;

                    DomainType tMaster = mMesh->block( tFacet->master()->block_id() )->domain_type() ;

                    if ( tFirst == nullptr )
                    {
                        tFirst = tFacet ;
                        tType  = tMaster ;
                        continue;
                    }

                    BELFEM_ERROR( tMaster == tType,
                        "interface sideset %lu mixes master domain types: facet %lu ( element %lu ) has a %s master, facet %lu ( element %lu ) a %s master; split the sideset so that each sideset has one master domain type",
                        ( long unsigned int ) tSideSet->id(),
                        ( long unsigned int ) tFirst->id(),
                        ( long unsigned int ) tFirst->master()->id(),
                        domain_type_label( tType ).c_str(),
                        ( long unsigned int ) tFacet->id(),
                        ( long unsigned int ) tFacet->master()->id(),
                        domain_type_label( tMaster ).c_str() );
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
                if ( tSet->shares_duplicates() )
                {
                    tSet->duplicate_nodes( mMaxNodeID, mSharedDuplicates );
                }
                else
                {
                    // every other set keeps private copies
                    Map< id_t, Node * > tPrivate ;
                    tSet->duplicate_nodes( mMaxNodeID, tPrivate );
                }
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

            Cell< id_t > tIDs ;

            for ( InterfaceSet * tSet : mSets )
            {
                Map< id_t, Node * > & tDuplicates = tSet->duplicate_map() ;
                Map< id_t, Node * > & tOriginals  = tSet->original_map() ;

                sorted_keys( tOriginals, tIDs );

                for ( id_t tID : tIDs )
                {
                    Node * A = tOriginals( tID );
                    if ( ! A->is_periodic() ) continue;

                    Node * C = tDuplicates( A->id() );
                    Node * B = A->periodic();

                    // the partner's duplicate, in this set or in a sibling set
                    // of the same treatment and class; a sibling that holds
                    // one but is excluded is remembered for the refusal
                    Node * D = nullptr ;
                    InterfaceSet * tOtherType = nullptr ;
                    InterfaceSet * tOtherTreatment = nullptr ;
                    if ( tDuplicates.key_exists( B->id() ) )
                    {
                        D = tDuplicates( B->id() );
                    }
                    else
                    {
                        for ( InterfaceSet * tOther : mSets )
                        {
                            if ( tOther == tSet || ! tOther->duplicate_map().key_exists( B->id() ) ) continue;

                            if ( ! this->pair_within( tSet, tOther ) )
                            {
                                if ( tSet->treatment() == tOther->treatment() )
                                {
                                    tOtherType = tOther ;
                                }
                                else
                                {
                                    tOtherTreatment = tOther ;
                                }
                                continue;
                            }

                            Node * tCandidate = tOther->duplicate_map()( B->id() );

                            BELFEM_ERROR( D == nullptr || D == tCandidate,
                                "periodic partner %lu of interface node %lu was duplicated in more than one interface set, pairing is ambiguous",
                                ( long unsigned int ) B->id(),
                                ( long unsigned int ) A->id() );

                            D = tCandidate ;
                        }
                    }

                    if ( D == nullptr )
                    {
                        this->refuse_unpaired_periodic_copy( tSet, A, B, tOtherType, tOtherTreatment );
                    }

                    // already tied, within a set or by this pass from the
                    // partner's side
                    if ( C->is_periodic() )
                    {
                        BELFEM_ERROR( C->periodic() == D,
                            "duplicate of interface node %lu is tied to another node than the duplicate of its periodic partner %lu",
                            ( long unsigned int ) A->id(),
                            ( long unsigned int ) B->id() );
                        continue;
                    }

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

            // every element of a conductor block takes the shared copy, also in
            // blocks that touch an interface node without owning an interface
            // facet
            if ( mSharedDuplicates.empty() ) return;

            for ( Block * tBlock : mMesh->blocks() )
            {
                if ( tBlock->domain_type() != DomainType::Conductor ) continue;

                for ( Element * tElement : tBlock->elements() )
                {
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        auto tIt = mSharedDuplicates.find( tElement->node( k )->id() );
                        if ( tIt != mSharedDuplicates.end() )
                        {
                            tElement->insert_node( tIt->second, k );
                        }
                    }
                }
            }
        }

        void
        InterfaceProcessor::add_duplicate_nodes_to_mesh()
        {
            // count duplicates: the shared copies once, plus the private
            // copies of every other set
            index_t tCount = mSharedDuplicates.size() ;
            for ( InterfaceSet * tSet : mSets )
            {
                if ( ! tSet->shares_duplicates() )
                {
                    tCount += tSet->number_of_nodes() ;
                }
            }
            Cell< Node * > tDuplicates( tCount, nullptr );

            // collect duplicates
            tCount = 0 ;
            for ( auto tPair : mSharedDuplicates )
            {
                tDuplicates( tCount++ ) = tPair.second ;
            }

            // temporaty map with originals
            Map< id_t, Node * > tOriginalMap ;
            for ( InterfaceSet * tSet : mSets )
            {
               tSet->add_duplicates( tCount, tDuplicates, tOriginalMap );
            }

            BELFEM_ASSERT( tCount == tDuplicates.size(),
                "collected %lu interface duplicates, expected %lu",
                ( long unsigned int ) tCount,
                ( long unsigned int ) tDuplicates.size() );

            // the maps were walked in hash order; the mesh takes the copies
            // in ascending id, which is their allocation order
            auto tByID = []( const Node * aA, const Node * aB ) { return aA->id() < aB->id(); };
            sort( tDuplicates, tByID );

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

        bool
        InterfaceProcessor::pair_within( InterfaceSet * aSet, InterfaceSet * aOther ) const
        {
            // a decoupled copy pairs with a decoupled copy only, a tie-weight
            // copy with a copy of a set with the same master type
            return aSet->treatment() == aOther->treatment()
                && ( aSet->treatment() == InterfaceTreatment::Decouple
                     || aSet->master_type() == aOther->master_type() );
        }

        void
        InterfaceProcessor::refuse_unpaired_periodic_copy( InterfaceSet * aSet,
                                                            Node * aNode,
                                                            Node * aPartner,
                                                            InterfaceSet * aOtherType,
                                                            InterfaceSet * aOtherTreatment ) const
        {
            BELFEM_ERROR( aOtherType == nullptr,
                "periodic partner %lu of interface node %lu is duplicated on an interface of another domain type ( node %lu: %s, node %lu: %s ): the periodic plane coincides with a change of domain type, which the periodic tie cannot represent; move the periodic plane into a section of one domain type",
                ( long unsigned int ) aPartner->id(),
                ( long unsigned int ) aNode->id(),
                ( long unsigned int ) aNode->id(),
                to_string( aSet->master_type() ).c_str(),
                ( long unsigned int ) aPartner->id(),
                to_string( aOtherType->master_type() ).c_str() );

            BELFEM_ERROR( aOtherTreatment == nullptr,
                "interface node %lu and its periodic partner %lu are duplicated on interfaces of different treatment ( node %lu: %s, node %lu: %s ): the decoupled copy touches a coil and carries no unknown, so the tie-weight copy has nothing to tie to; move the periodic plane away from the coil",
                ( long unsigned int ) aNode->id(),
                ( long unsigned int ) aPartner->id(),
                ( long unsigned int ) aNode->id(),
                treatment_label( aSet->treatment() ),
                ( long unsigned int ) aPartner->id(),
                treatment_label( aOtherTreatment->treatment() ) );

            BELFEM_ERROR( false,
                "periodic partner %lu of interface node %lu was not duplicated in any interface set",
                ( long unsigned int ) aPartner->id(),
                ( long unsigned int ) aNode->id() );
        }

    }
}
