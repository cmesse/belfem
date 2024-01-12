//
// Created by Christian Messe on 2019-07-28.
//
#include "cl_SideSet.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        SideSet::SideSet( const id_t aID, const index_t aNumFacets ) :
                mID( aID ),
                mFacetCounter( 0 )
        {
            mFacets.set_size( aNumFacets, nullptr );
            mLabel = aID < 10 ?
                     sprint( "SideSet_0%u", ( unsigned int ) aID ) :
                     sprint( "SideSet_%u", ( unsigned int ) aID );
        }

//------------------------------------------------------------------------------

        SideSet::~SideSet()
        {
            // only clean map, facets are deleted bz mesh
            mFacetMap.clear() ;
        }

//------------------------------------------------------------------------------

        void
        SideSet::insert_facet( Facet * aFacet )
        {
            BELFEM_ASSERT( mFacetCounter < mFacets.size(),
                           "SideSet is full" );

            mFacets( mFacetCounter++ ) = aFacet;

            mFacetMap[ aFacet->id() ] = aFacet ;
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
            mFacetMap.clear();
            mFacetCounter = 0 ;
            mFacets.clear();
        }

//------------------------------------------------------------------------------
    }
}