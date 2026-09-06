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

#include "cl_Mesh_SourceExpander.hpp"
#include "fn_entity_type.hpp"
namespace belfem
{
    namespace mesh
    {
        SourceExpander::SourceExpander( Mesh *aMesh ) :
           mMesh( aMesh )
        {
            mNodeBitset    = new DynamicBitset( aMesh->number_of_nodes() ) ;
            mEdgeBitset    = new DynamicBitset( aMesh->number_of_edges() ) ;
            mFaceBitset    = new DynamicBitset( aMesh->number_of_faces() ) ;
            mElementBitset = new DynamicBitset( aMesh->number_of_elements() ) ;
            mFacetBitset   = new DynamicBitset( aMesh->number_of_facets() ) ;
        }

        SourceExpander::~SourceExpander()
        {
            if ( mNodeBitset != nullptr )
            {
                delete mNodeBitset ;
            }
            if ( mEdgeBitset != nullptr )
            {
                delete mEdgeBitset ;
            }
            if ( mFaceBitset != nullptr )
            {
                delete mFaceBitset ;
            }
            if ( mElementBitset != nullptr )
            {
                delete mElementBitset ;
            }
            if ( mFacetBitset != nullptr )
            {
                delete mFacetBitset ;
            }
        }

        void
        SourceExpander::reset_bitsets()
        {
            mNodeBitset->reset() ;
            mEdgeBitset->reset() ;
            mFaceBitset->reset() ;
            mElementBitset->reset() ;
            mFacetBitset->reset() ;
        }

        void
        SourceExpander::expand_sources( Basis * aBasis )
        {

            if ( ! aBasis->is_hanging() ) return ;

            // check if this source is already expanded
            bool tIsExpanded = true ;

            for ( uint k=0; k<aBasis->number_of_sources(); ++k )
            {
                if ( aBasis->source( k )->is_hanging() )
                {
                    tIsExpanded = false ;
                    break ;
                }
            }

            // nothing to do for this source
            if ( tIsExpanded ) return;

            this->reset_bitsets() ;
            this->flag_all_sources( aBasis ) ;

            // this will contain only independent sources
            mSources.clear();

            // collect nodes
            Cell< Node * > & tNodes = mMesh->nodes() ;
            mNodeBitset->where( mBasisIndices );
            for ( index_t k : mBasisIndices )
            {
                mSources.push( tNodes(k) );
            }

            // collect edges
            Cell< Edge * > & tEdges = mMesh->edges() ;
            mEdgeBitset->where( mBasisIndices );
            for ( index_t k : mBasisIndices )
            {
                mSources.push( tEdges(k) );
            }

            //  collect faces
            Cell< Face * > & tFaces = mMesh->faces() ;
            mFaceBitset->where( mBasisIndices );
            for ( index_t k : mBasisIndices )
            {
                mSources.push( tFaces(k) );
            }

            // collect elements
            Cell< Element * > & tElements = mMesh->elements() ;
            mElementBitset->where( mBasisIndices );
            for ( index_t k : mBasisIndices )
            {
                mSources.push( tElements(k) );
            }

            // collect facets
            Cell< Facet * > & tFacets = mMesh->facets() ;
            mFacetBitset->where( mBasisIndices );
            for ( index_t k : mBasisIndices )
            {
                mSources.push( tFacets(k) );
            }
            mBasisIndices.clear();

            // initiate weights
            mWeights.clear();
            for ( Basis * tBasis : mSources )
            {
                mWeights[ tBasis ] = 0.0 ;
            }

            // expand weights
            this->expand_weights( aBasis ) ;

            // convert map to vector
            mWork.set_size( mSources.size() );

            for ( index_t k = 0; k< mSources.size(); ++k )
            {
                mWork( k ) = mWeights( mSources( k ) );
            }
#if !defined( NDEBUG ) || defined( DEBUG )
            // sanity check
            for ( Basis * tSource : mSources )
            {
                BELFEM_ERROR( ! tSource->is_hanging(), "Source is hanging" );
            }
#endif

            // now we need to reallocate the basis
            aBasis->reset_source_container();
            aBasis->set_sources( mSources, mWork );
        }


    }
}