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

#include "cl_SideSet.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        SideSet::SideSet( const id_t aID, const index_t aNumFacets ) :
                mID( aID )
        {

            mLabel = aID < 10 ?
                     sprint( "SideSet_0%u", ( unsigned int ) aID ) :
                     sprint( "SideSet_%u", ( unsigned int ) aID );

            if ( aNumFacets > 0 ) mFacets.set_size( aNumFacets, nullptr );
        }

//------------------------------------------------------------------------------

        SideSet::~SideSet()
        {
            for ( Facet * tFacet : mFacets )
            {
                delete tFacet ;
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::insert_facet( Facet * aFacet )
        {
            BELFEM_ASSERT( mFacetCounter < mFacets.size(),
                          "SideSet is full" );

            mFacets( mFacetCounter++ ) = aFacet;
        }

//------------------------------------------------------------------------------

        void
        SideSet::unflag_all_nodes()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->unflag_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::flag_all_nodes()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->flag_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::flag_corner_nodes()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->flag_corner_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::unflag_all_facets()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->unflag();
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::flag_all_facets()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->flag();
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::flag_edges()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->element()->flag_edges() ;
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::unflag_edges()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->element()->unflag_edges() ;
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::flag_faces()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->element()->flag_faces() ;
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::unflag_faces()
        {
            for( Facet * tFacet : mFacets )
            {
                tFacet->element()->unflag_faces() ;
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::collect_nodes()
        {
            this->unflag_all_nodes() ;

            // count nodes
            index_t tCount = 0 ;
            for( Facet * tFacet : mFacets )
            {
                uint tNumNodes = tFacet->number_of_nodes() ;
                for( uint k=0; k<tNumNodes; ++k )
                {
                    if( ! tFacet->node( k )->is_flagged() )
                    {
                        ++tCount ;
                        tFacet->node( k )->flag() ;
                    }
                }
            }

            // init container
            mNodes.set_size( tCount, nullptr );

            // reset counter and collect
            tCount = 0 ;
            this->unflag_all_nodes() ;
            for( Facet * tFacet : mFacets )
            {
                uint tNumNodes = tFacet->number_of_nodes() ;
                for( uint k=0; k<tNumNodes; ++k )
                {
                    if( ! tFacet->node( k )->is_flagged() )
                    {
                        mNodes( tCount++ ) = tFacet->node( k );
                        tFacet->node( k )->flag() ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::reset_node_container()
        {
            mNodes.clear();
        }

//------------------------------------------------------------------------------

        void
        SideSet::reset_facet_container()
        {
            mFacetCounter = 0 ;
            mFacets.clear();
        }

//------------------------------------------------------------------------------
    }
}