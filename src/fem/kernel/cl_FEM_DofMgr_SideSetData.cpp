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
#include "cl_FEM_DofMgr_SideSetData.hpp"
#include "commtools.hpp"
#include "cl_FEM_DofMgr_BlockData.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofMgr_DofData.hpp"
#include "cl_FEM_DofManager.hpp"

#include "fn_max.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//-----------------------------------------------------------------------------

            SideSetData::SideSetData(  DofManager * aParent ) :
                    mParent( aParent ),
                    mKernel( aParent->parent() ),
                    mMesh( aParent->parent()->mesh() ),
                    mCommRank( aParent->rank() )
            {

            }

//-----------------------------------------------------------------------------

            SideSetData::~SideSetData()
            {
                this->reset();
            }

//-----------------------------------------------------------------------------

            void
            SideSetData::reset()
            {
                // delete the maps
                mSideSetMap.clear() ;

                for( SideSet * tSideSet : mSideSets )
                {
                    delete tSideSet ;
                }
                mSideSets.clear() ;

                if( mEmptySideset != nullptr )
                {
                    delete mEmptySideset ;
                    mEmptySideset = nullptr ;
                }
            }

//-----------------------------------------------------------------------------

            void
            SideSetData::create_sidesets()
            {
                index_t tSideSetCount = 0;

                // count sidesets that exist
                for ( id_t  tSideSetID : mParent->iwg()->selected_sidesets() )
                {
                    if( mMesh->sideset_exists( tSideSetID ) )
                    {
                        // get sideset
                        mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                        // count owned facets
                        index_t tCount = 0;
                        for ( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            if ( tFacet->owner() == mCommRank )
                            {
                                ++tCount;
                            }
                        }
                        if( tCount > 0 )
                        {
                            ++tSideSetCount ;
                        }
                    }
                }

                mSideSets.set_size( tSideSetCount, nullptr );


                // reset the counter
                tSideSetCount = 0 ;

                for ( id_t  tSideSetID : mParent->iwg()->selected_sidesets() )
                {
                    if( mMesh->sideset_exists( tSideSetID ) )
                    {
                        // get sideset
                        mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                        // count owned facets
                        index_t tCount = 0;
                        for ( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            if ( tFacet->owner() == mCommRank )
                            {
                                ++tCount;
                            }
                        }

                        if ( tCount > 0 )
                        {
                            // allocate temporary container for facets
                            Cell< mesh::Facet * > tFacets( tCount, nullptr );

                            // reset counter
                            tCount = 0;

                            // collect facets
                            for ( mesh::Facet * tFacet : tSideSet->facets() )
                            {
                                if ( tFacet->owner() == mCommRank )
                                {
                                    tFacets( tCount++ ) = tFacet;
                                }
                            }

                            // get the type of the sideset
                            switch( mParent->iwg()->sideset_type( tSideSetID ) )
                            {
                                default :
                                {
                                    // create the sideset
                                    mSideSets( tSideSetCount ) = new SideSet(
                                            mParent,
                                            tSideSet->id(),
                                            tFacets );
                                    break ;
                                }
                            }

                            // query activation mode from IWG based on domain type
                            DomainType tDomainType = tSideSet->domain_type();
                            GroupActivationMode tActivationMode = mParent->iwg()->sideset_activation_mode( tDomainType );
                            mSideSets( tSideSetCount )->set_activation_mode( tActivationMode );

                            // add sideset to map
                            mSideSetMap[ tSideSet->id() ] = mSideSets( tSideSetCount++ );
                        }
                    }
                }

                // the empty sideset is a dummy that is exposed if
                // a set is accessed that doesn't exist on this proc
                Cell< mesh::Facet * > tEmpty;

                // create the empty sideset
                mEmptySideset = new SideSet(
                        mParent,
                        0,
                        tEmpty );
            }

//-----------------------------------------------------------------------------

            void
            SideSetData::collect_wetted_sidesets()
            {
                // get number of wetted sidesets on this proc
                if( mParent->iwg()->has_convection() )
                {
                    uint tNumDofs = mParent->iwg()->dof_entity_types().length() ;

                    index_t tCount = 0 ;
                    for( SideSet * tSideSet : mSideSets )
                    {
                        for ( uint k = 0; k < tNumDofs; ++k )
                        {
                            if ( tSideSet->bc_type( k ) == BoundaryConditionImposing::Neumann ||
                                 tSideSet->bc_type( k ) == BoundaryConditionImposing::Alpha )
                            {
                                ++tCount;
                                break; // <- break, because we need each sideset only once
                            }
                        }
                    }

                    if( tCount > 0 )
                    {

                        Vector< id_t > tWettedSidesets( tCount );
                        tCount = 0;
                        for ( SideSet * tSideSet : mSideSets )
                        {
                            for ( uint k = 0; k < tNumDofs; ++k )
                            {
                                if ( tSideSet->bc_type( k ) == BoundaryConditionImposing::Neumann ||
                                     tSideSet->bc_type( k ) == BoundaryConditionImposing::Alpha )
                                {
                                    tWettedSidesets( tCount++ ) = tSideSet->id();
                                    break; // <- break, because we need each sideset only once
                                }
                            }
                        }

                        mParent->iwg()->set_wetted_sidesets( tWettedSidesets );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SideSetData::count_wetted_nodes()
            {
                // unflag all nodes
                mMesh->unflag_all_nodes() ;

                // loop over all wetted sidesets
                for( id_t tID : mParent->iwg()->wetted_sidesets() )
                {
                    // select sideset
                    if( this->sideset_exists( tID ) )
                    {
                        mMesh->sideset( tID )->flag_all_nodes() ;
                    }
                }

                // reset node counter
                mNumberOfConvectionNodes = 0 ;

                // count nodes
                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        ++mNumberOfConvectionNodes ;
                    }
                }
            }

//------------------------------------------------------------------------

            void
            SideSetData::create_alpha_fields()
            {
                // check if an alpha BC (special BC for thermal convection) exists on any proc:

                uint tLocalFlag = 0 ;
                uint tGlobalFlag = 0 ;
                
                for( SideSet * tSideSet : mSideSets )
                {
                    uint tNumDofs = tSideSet->number_of_boundary_conditions();

                    for( uint k=0; k<tNumDofs; ++k )
                    {
                        if( tSideSet->bc_type( k ) == BoundaryConditionImposing::Alpha )
                        {
                            tLocalFlag = 1 ;
                            break ;
                        }
                    }
                    if( tLocalFlag > 0 )
                    {
                        break ;
                    }
                }

                if( mCommRank == 0 )
                {
                    Vector< uint > tFlags( comm_size(), 0 );
                    collect( tFlags ) ;

                    if( max( tFlags ) > 0 || tLocalFlag > 0 )
                    {
                        tGlobalFlag = 1 ;
                    }

                    tFlags.fill( tGlobalFlag );
                    distribute( tFlags );

                }
                else
                {
                    send( tLocalFlag );
                    receive( tGlobalFlag );
                }

                if( tGlobalFlag > 0 )
                {
                    // add fields to IWG
                    mParent->iwg()->add_fields( { "alpha", "Tinf" } );
                }
            }

//------------------------------------------------------------------------------

            void
            SideSetData::set_boundary_conditions()
            {
                comm_barrier() ;

                // loop over all sidesets and update data
                for ( mesh::SideSet * tSideset : mMesh->sidesets() )
                {
                    this->sideset( tSideset->id() )->set_boundary_conditions() ;
                }
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */




