//
// Created by christian on 7/14/21.
//

#include "cl_FEM_DofMgr_SideSetData.hpp"
#include "commtools.hpp"
#include "cl_FEM_DofMgr_BlockData.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofMgr_DofData.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Cut.hpp"
#include "cl_FEM_Shell.hpp"
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
                    mMyRank( aParent->rank() )
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
                            if ( tFacet->owner() == mMyRank )
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
                            if ( tFacet->owner() == mMyRank )
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
                                if ( tFacet->owner() == mMyRank )
                                {
                                    tFacets( tCount++ ) = tFacet;
                                }
                            }

                            // get the type of the sideset
                            switch( mParent->iwg()->sideset_type( tSideSetID ) )
                            {
                                case( DomainType::Cut ) :
                                {
                                    // create the sideset
                                    mSideSets( tSideSetCount ) = new Cut(
                                            mParent,
                                            tSideSet->id(),
                                            tFacets );
                                    break ;
                                }
                                case( DomainType::ThinShell ) :
                                {
                                    mSideSets( tSideSetCount ) = new Shell(
                                            mParent,
                                            tSideSet->id(),
                                            tFacets );
                                    break ;
                                }
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
            SideSetData::detect_wetted_sidesets()
            {
                // get number of wetted sidesets on this proc
                if( mParent->iwg()->wetted_sidesets().length() > 0 )
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

                if( mMyRank == mKernel->master() )
                {
                    Vector< uint > tFlags( mKernel->comm_table().length(), 0 );
                    receive( mKernel->comm_table(), tFlags ) ;

                    if( max( tFlags ) > 0 || tLocalFlag > 0 )
                    {
                        tGlobalFlag = 1 ;
                    }

                    tFlags.fill( tGlobalFlag );
                    send( mKernel->comm_table(), tFlags );

                }
                else
                {
                    send( mKernel->master(), tLocalFlag );
                    receive( mKernel->master(), tGlobalFlag );
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




