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

#include "cl_CutProcessorManual.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_IwgFactory.hpp"

namespace belfem
{
    namespace mesh
    {
        CutProcessorManual::CutProcessorManual( Mesh * aMesh,
                               const Vector< id_t > & aAirBlocks,
                               const Vector< id_t > & aCutSideSets,
                               const Vector< id_t > & aCutMasterVolumes ) :
            mMesh( aMesh ),
            mAirBlocks( aAirBlocks ),
            mCutSideSets( aCutSideSets ),
            mCutMasterVolumes( aCutMasterVolumes )
        {
            mSets.set_size( aCutSideSets.length(), nullptr );
            id_t tMaxNodeID = mMesh->max_node_id() ;
            mMesh->unflag_all_nodes();
            mMesh->unflag_all_elements() ;

            index_t tNcuts = aCutSideSets.length() ;

            mIncidentMatrix.set_size( tNcuts, tNcuts, 0.0 );

            index_t tCount = 0 ;

            mAbstractNodes.set_size( tNcuts, nullptr );
            for ( index_t i=0; i<tNcuts; ++i )
            {
                mSets( i ) = this->create_set( i, tMaxNodeID );
                mAbstractNodes( i ) = mSets( i )->mAbstractNode ;
                tCount += mSets( i )->mDuplicateNodes.size() ;
                mIncidentMatrix( i, i ) = 1.0 ;
            }

            // collect nodes
            mOriginalNodes.set_size( tCount, nullptr );
            mDuplicateNodes.set_size( tCount, nullptr );

            tCount = 0 ;
            for ( index_t i=0; i<tNcuts; ++i )
            {
                for ( Node * tNode : mSets( i )->mOriginalNodes )
                {
                    mOriginalNodes( tCount++ ) = tNode ;
                }
            }

            tCount = 0 ;
            for ( index_t i=0; i<tNcuts; ++i )
            {
                for ( Node * tNode : mSets( i )->mDuplicateNodes )
                {
                    mDuplicateNodes( tCount++ ) = tNode ;
                }
            }

            // add nodes to mesh
            append( mMesh->nodes(), mAbstractNodes );
            append( mMesh->nodes(), mDuplicateNodes );

            mMesh->unfinalize();
            mMesh->finalize();

        }

        CutProcessorManual::~CutProcessorManual()
        {
            for ( ManualSet * tSet : mSets )
            {
                delete tSet ;
            }
        }

        CutProcessorManual::ManualSet *
        CutProcessorManual::create_set( const index_t aIndex, id_t & aMaxNodeID )
        {

            SideSet * tSideSet = mMesh->sideset( mCutSideSets( aIndex ) );

            ManualSet * aSet = new ManualSet();

            Cell< Element * > & tElements = mMesh->block( mCutMasterVolumes( aIndex ) )->elements() ;

            index_t tCount = 0 ;

            for ( Facet * tFacet : tSideSet->facets() )
            {
                tFacet->flag_nodes() ;
            }

            tCount = 0 ;
            for ( Element * tElement : tElements )
            {
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

            aSet->mElements.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( Element * tElement : tElements )
            {
                if ( tElement->is_flagged() )
                {
                    aSet->mElements( tCount++ ) = tElement ;
                    tElement->unflag() ;
                }
            }

            // create abstract node
            aSet->mAbstractNode = new Node( ++aMaxNodeID, 0.0, 0.0, 0.0 );

            std::cout << "created abstract node " << aSet->mAbstractNode->id() << std::endl;

            // count original nodes
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }

            aSet->mOriginalNodes.set_size( tCount, nullptr );

            tCount = 0 ;

            // collect original nodes
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    aSet->mOriginalNodes( tCount++ ) = tNode ;
                }
            }

            // create duplicates
            aSet->mDuplicateNodes.set_size( tCount, nullptr );
            Cell< Node * > tSources( 2, nullptr );
            tSources( 0 ) = aSet->mAbstractNode ;
            Vector< real > tWeights( 2, 1.0 );

            for ( Node * tOrg : aSet->mOriginalNodes )
            {
                Node * tDup = new Node( ++aMaxNodeID, tOrg->x(), tOrg->y(), tOrg->z() );
                tSources( 1 ) = tOrg ;
                tDup->set_sources( tSources, tWeights );

                aSet->mDuplicateNodes( tOrg->index() ) = tDup ;
            }

            // relink nodes
            std::sort( tElements.begin(), tElements.end(), []( Element * a, Element * b ) { return a->id() < b->id(); } );

            for ( Element * tElement : tElements )
            {
                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    if ( tElement->node( k )->is_flagged() )
                    {
                        tElement->insert_node( aSet->mDuplicateNodes( tElement->node( k )->index() ), k );
                    }
                }
            }

            for ( Node * tNode : aSet->mOriginalNodes )
            {
                tNode->unflag() ;
            }
            return aSet ;
        }
    }
}
