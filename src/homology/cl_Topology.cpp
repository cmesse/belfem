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

#include "cl_Topology.hpp"

#include "cl_FEM_DofMgr_SideSetData.hpp"
#include "fn_combine.hpp"
#include "fn_append.hpp"

namespace belfem
{
    namespace mesh
    {
        Topology::Topology( Mesh * aMesh ) :
            mCommRank( comm_rank() ),
            mMesh( aMesh )
        {

        }

        Topology::~Topology()
        {
            for ( auto tPair : mTypeMap )
            {
                delete tPair.second ;
            }
        }

        void
        Topology::run()
        {
            if ( mCommRank == 0 )
            {
                this->collect_block_and_sideset_types();
                this->detect_sideset_types();
                this->select_blocks();
                this->select_sidesets();
            }
        }

        void
        Topology::run_on_enriched_mesh()
        {
            if ( mCommRank == 0 )
            {
                this->collect_enrichment_ids();
                this->collect_block_and_sideset_types();

                // no detect_sideset_types() here: the reloaded mesh already
                // carries the final sideset classification, and the cut
                // sidesets must keep DomainType::Default like on the fresh
                // path, where they are created after this map is built

                this->select_blocks();
                this->select_sidesets();
            }
        }

        void
        Topology::collect_enrichment_ids()
        {
            index_t tCount = 0 ;

            for ( ThinShell * tShell : mMesh->thin_shells() )
            {
                for ( Block * tBlock : tShell->blocks() )
                {
                    mEnrichedBlocks[ tBlock->id() ] = tCount++ ;
                }

                mEnrichedSideSets[ tShell->id() ] = tCount++ ;

                if ( tShell->ghost_id() != gNoID )
                {
                    mEnrichedSideSets[ tShell->ghost_id() ] = tCount++ ;
                }
            }
        }

        void
        Topology::synchronize_maps()
        {
            comm_barrier() ;

            Vector< id_t > tData ;
            if ( mCommRank == 0 )
            {
                index_t n = static_cast< index_t >( DomainType::UNDEFINED ) + 1 ;

                index_t tCount = 2 * n ;

                for ( auto tPair : mTypeMap )
                {
                    tCount += tPair.second->length() ;
                }

                tData.set_size( tCount );
                tCount = 0 ;
                for ( auto tPair : mTypeMap )
                {
                    tData( tCount++ ) = static_cast< id_t >( tPair.first ) ;
                    tData( tCount++ ) = tPair.second->length() ;
                    for ( id_t tId : *tPair.second )
                    {
                        tData( tCount++ ) = tId ;
                    }
                }
            }
            else
            {
                for ( auto tPair : mTypeMap )
                {
                    delete tPair.second ;
                }
                mTypeMap.clear() ;
            }

            broadcast( tData );

            if ( mCommRank != 0 )
            {
                index_t tCount = 0 ;

                index_t n = static_cast< index_t >( DomainType::UNDEFINED ) + 1 ;

                for ( index_t i=0; i<n; ++i )
                {
                    // get the domain type
                    DomainType tDomainType = static_cast< DomainType >( tData( tCount++ ) ) ;

                    // get the number of entries
                    index_t tLength = tData( tCount++ ) ;

                    // allocate vector
                    mTypeMap[ tDomainType ] = new Vector< id_t >( tLength );

                    Vector< id_t > & tVector = *mTypeMap[ tDomainType ] ;

                    // copy the data
                    for ( index_t j=0; j<tLength; ++j )
                    {
                        tVector( j ) = tData( tCount++ ) ;
                    }
                }
            }

            broadcast( mPhiBlockIDs );
            broadcast( mNonPhiBlockIDs );
            broadcast( mPhiBoundaryIDs );
            broadcast( mPhiInterfaceIDs );

            this->update_block_map();
        }

        void
        Topology::collect_block_and_sideset_types()
        {
            BELFEM_ASSERT( mCommRank == 0, "this function must be called on rank 0" );

            index_t n = static_cast< index_t >( DomainType::UNDEFINED ) + 1 ;

            // count blocks and sidesets
            Vector< index_t > tCount( n, 0 );
            for ( Block * tBlock : mMesh->blocks() )
            {
                if ( mEnrichedBlocks.key_exists( tBlock->id() ) ) continue ;
                tCount( static_cast< index_t >( tBlock->domain_type() ) )++ ;
            }
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( mEnrichedSideSets.key_exists( tSideSet->id() ) ) continue ;
                tCount( static_cast< index_t >( tSideSet->domain_type() ) )++ ;
            }

            // allocate type maps
            for ( index_t i=0; i<n; ++i )
            {
                mTypeMap[ static_cast< DomainType >( i ) ] = new Vector< id_t >( tCount( i ) );
            }

            // collect types
            tCount.fill( 0 );
            for ( Block * tBlock : mMesh->blocks() )
            {
                if ( mEnrichedBlocks.key_exists( tBlock->id() ) ) continue ;

                // get the vector
                Vector< id_t > & tVector = *mTypeMap[ tBlock->domain_type() ] ;

                tVector( tCount( static_cast< index_t >( tBlock->domain_type() ) )++ ) = tBlock->id() ;
            }

            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( mEnrichedSideSets.key_exists( tSideSet->id() ) ) continue ;

                // get the vector
                Vector< id_t > & tVector = *mTypeMap[ tSideSet->domain_type() ] ;

                tVector( tCount( static_cast< index_t >( tSideSet->domain_type() ) )++ ) = tSideSet->id() ;
            }
        }

        void
        Topology::detect_sideset_types()
        {
            BELFEM_ASSERT( mCommRank == 0, "this function must be called on rank 0" );

            Cell< SideSet * > & tSideSets = mMesh->sidesets() ;

            for ( SideSet * tSideSet : tSideSets )
            {
                if ( tSideSet->domain_type() == DomainType::Default )
                {
                    if ( tSideSet->number_of_facets() == 0 )
                    {
                        tSideSet->set_domain_type( DomainType::Inactive );
                    }
                    else
                    {
                        // get the first facet
                        Facet * tFacet = tSideSet->facets()( 0 );

                        DomainType tMaster = mMesh->block( tFacet->master()->block_id() )->domain_type() ;

                        if ( tFacet->has_slave() )
                        {
                            DomainType tSlave = mMesh->block( tFacet->slave()->block_id() )->domain_type() ;
                            tSideSet->set_domain_type( this->sideset_type( tMaster, tSlave ) );
                        }
                        else
                        {
                            // One-sided buffer faces are phi boundaries:
                            // sideset_type(Buffer) classifies them as
                            // BufferAntiSymmetry.
                            tSideSet->set_domain_type( this->sideset_type( tMaster ) );
                        }
                    }
                }
            }
        }

        DomainType
        Topology::sideset_type( const DomainType aMaster )
        {
            switch ( aMaster )
            {
                case DomainType::Air :
                {
                    return DomainType::AirAntiSymmetry ;
                }
                case DomainType::Buffer :
                {
                    return DomainType::BufferAntiSymmetry ;
                }
                case DomainType::Conductor :
                {
                    return DomainType::ConductorAntiSymmetry ;
                }
                case DomainType::Ferro :
                {
                    return DomainType::FerroAntiSymmetry ;
                }
                default:
                {
                    return DomainType::Inactive ;
                }
             }
        }

        DomainType
        Topology::sideset_type( const DomainType aMaster, const DomainType aSlave )
        {
            // Buffer is a phi-block (Topology::select_blocks) and is treated
            // as Air for interface classification: a Buffer/Conductor face
            // is wired identically to an Air/Conductor face by the
            // hanging-edges machinery in MaxwellFactory::create_cuts (the
            // master is the conductor side, the slave is a phi-side block,
            // and that's the only thing the InterfaceCondAir path requires).
            const DomainType tMaster = ( aMaster == DomainType::Buffer ) ? DomainType::Air : aMaster ;
            const DomainType tSlave  = ( aSlave  == DomainType::Buffer ) ? DomainType::Air : aSlave  ;

            switch ( tMaster )
            {
                case DomainType::Air :
                {
                    switch ( tSlave )
                    {
                        case DomainType::Air :
                        {
                            return DomainType::Inactive ;
                        }
                        case DomainType::Conductor :
                        {
                            return DomainType::InterfaceCondAir ;
                        }
                        case DomainType::Ferro :
                        {
                            return DomainType::InterfaceFerroAir ;
                        }
                        case DomainType::Coil :
                        {
                            return DomainType::InterfaceAirCoil ;
                        }
                        default:
                        {
                            return DomainType::AirAntiSymmetry ;
                        }
                    }
                }
                case DomainType::Coil :
                {
                    switch ( tSlave )
                    {
                        case DomainType::Air :
                        {
                            return DomainType::InterfaceAirCoil ;
                        }
                        case DomainType::Conductor :
                        {
                            BELFEM_ERROR( false, "a conductor and a coil must not touch!" );
                            return DomainType::UNDEFINED ;
                        }
                        case DomainType::Ferro :
                        {
                            return DomainType::InterfaceFerroCoil ;
                        }
                        default:
                        {
                            return DomainType::Inactive ;
                        }
                    }
                }
                case DomainType::Conductor :
                {
                    switch ( tSlave )
                    {
                        case DomainType::Air :
                        {
                            return DomainType::InterfaceCondAir ;
                        }
                        case DomainType::Conductor :
                        {
                            return DomainType::Inactive ;
                        }
                        case DomainType::Ferro :
                        {
                            return DomainType::InterfaceCondFerro ;
                        }
                        case DomainType::Coil :
                        {
                            BELFEM_ERROR( false, "a conductor and a coil must not touch!" );
                            return DomainType::UNDEFINED ;
                        }
                        default:
                        {
                            return DomainType::ConductorAntiSymmetry ;
                        }
                    }
                }
                case DomainType::Ferro :
                {
                    switch ( tSlave )
                    {
                        case DomainType::Air :
                        {
                            return DomainType::InterfaceFerroAir ;
                        }
                        case DomainType::Conductor :
                        {
                            return DomainType::InterfaceCondFerro ;
                        }
                        case DomainType::Ferro :
                        {
                            return DomainType::Inactive ;
                        }
                        default:
                        {
                            return DomainType::FerroAntiSymmetry ;
                        }
                    }
                }
                default:
                {
                    switch ( tSlave )
                    {
                        case DomainType::Air :
                        {
                            return DomainType::AirAntiSymmetry ;
                        }
                        case DomainType::Conductor :
                        {
                            return DomainType::ConductorAntiSymmetry ;
                        }
                        case DomainType::Ferro :
                        {
                            return DomainType::FerroAntiSymmetry ;
                        }
                        default:
                        {
                            return DomainType::Inactive ;
                        }
                    }
                }
            }
        }

        void Topology::select_blocks()
        {
            BELFEM_ASSERT( mCommRank == 0, "this function must be called on rank 0" );

            Cell< Block * > tPhiBlocks ;
            Cell< Block * > tNonPhiBlocks ;

            for ( Block * tBlock : mMesh->blocks() )
            {
                // thin-shell layer and buffer blocks don't exist yet when the
                // fresh path builds this map; skip them on a reloaded mesh
                if ( mEnrichedBlocks.key_exists( tBlock->id() ) ) continue ;

                switch ( tBlock->domain_type() )
                {
                    case DomainType::Air :
                    case DomainType::Buffer :
                    case DomainType::Ferro :
                    {
                        tPhiBlocks.push( tBlock );
                        break ;
                    }
                    case DomainType::Conductor :
                    case DomainType::Coil :
                    {
                        tNonPhiBlocks.push( tBlock );
                        break ;
                    }
                    default:
                    {
                        // pass
                    }
                }
            }

            mPhiBlockIDs.set_size( tPhiBlocks.size() );
            mNonPhiBlockIDs.set_size( tNonPhiBlocks.size() );

            index_t tCount = 0 ;
            for ( Block * tBlock : tPhiBlocks )
            {
                mPhiBlockIDs( tCount++ ) = tBlock->id() ;
            }
            tCount = 0 ;
            for ( Block * tBlock : tNonPhiBlocks )
            {
                mNonPhiBlockIDs( tCount++ ) = tBlock->id() ;
            }
        }

        void
        Topology::select_sidesets()
        {
            BELFEM_ASSERT( mCommRank == 0, "this function must be called on rank 0" );

            Cell< SideSet * > tPhiInterfaces ;
            Cell< SideSet * > tPhiBoundaries ;
            Cell< SideSet * > tPhiPeriodic ;
            Cell< SideSet * > & tSideSets = mMesh->sidesets() ;

            for ( SideSet * tSideSet : tSideSets )
            {
                switch ( tSideSet->domain_type() )
                {
                    case DomainType::AirAntiSymmetry :
                    case DomainType::AirSymmetry :
                    case DomainType::BufferSymmetry :
                    case DomainType::BufferAntiSymmetry :
                    case DomainType::BackgroundField :
                    case DomainType::FerroAntiSymmetry :
                    case DomainType::FerroSymmetry :


                    {
                        tPhiBoundaries.push( tSideSet );
                        break ;
                    }
                    case DomainType::Periodic :
                    case DomainType::AirPeriodic :
                    case DomainType::BufferPeriodic :
                    case DomainType::FerroPeriodic :
                    {
                        tPhiPeriodic.push( tSideSet );
                        break;
                    }
                    case DomainType::InterfaceCondAir :
                    case DomainType::InterfaceCondFerro :
                    case DomainType::InterfaceAirCoil :
                    case DomainType::InterfaceFerroCoil :
                    case DomainType::ThinShell  :
                    {
                        tPhiInterfaces.push( tSideSet );
                        break ;
                    }
                    default:
                    {
                        // pass
                    }
                }
            }


            index_t tCount = 0 ;
            mPhiBoundaryIDs.set_size( tPhiBoundaries.size() );
            for ( SideSet * tSideSet : tPhiBoundaries )
            {
                mPhiBoundaryIDs( tCount++ ) = tSideSet->id() ;
            }


            tCount = 0 ;
            mPhiInterfaceIDs.set_size( tPhiInterfaces.size() );
            for ( SideSet * tSideSet : tPhiInterfaces )
            {
                mPhiInterfaceIDs( tCount++ ) = tSideSet->id() ;
            }

            tCount = 0 ;
            mPhiPeriodicIDs.set_size( tPhiPeriodic.size() );
            for ( SideSet * tSideSet : tPhiPeriodic )
            {
                mPhiPeriodicIDs( tCount++ ) = tSideSet->id() ;
            }
        }

        void
        Topology::update_block_map()
        {
            mBlockTypes.clear() ;

            for ( auto tGroup : mTypeMap )
            {
                switch ( tGroup.first )
                {
                    case DomainType::Air :
                    case DomainType::Buffer :
                    case DomainType::Ferro :
                    case DomainType::Conductor :
                    // case DomainType::Coil : no coil, because not used computation
                    {
                        const Vector< id_t > & tVector = *tGroup.second ;

                        for ( id_t tId : tVector )
                        {
                            mBlockTypes[ tId ] = tGroup.first ;
                        }
                        break ;
                    }
                    default:
                    {
                        continue;
                    }
                }
            }
        }


    }
}
