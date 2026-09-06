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

#include "cl_FEM_Dof.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID, const uint aType, mesh::Node * aNode ) :
            Vertex(),
            mTypeID( aType ),
            mMeshBasis( aNode ),
            mIndexOnEntity( 0 ),
            mIndexOnField( aNode->index() )

        {
            this->set_id( aID );
            this->set_owner( aNode->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( id_t aID,
                const uint aType,
                mesh::Edge * aEdge,
                const uint aIndexOnEdge,
                const index_t aDofIndexOnField  ) :
                Vertex(),
                mTypeID( aType ),
                mMeshBasis( aEdge ),
                mIndexOnEntity( aIndexOnEdge ),
                mIndexOnField( aDofIndexOnField )
        {
            BELFEM_ERROR( aDofIndexOnField != gNoIndex, "Dof index on field must not be gNoIndex" );

            this->set_id( aID );
            this->set_owner( aEdge->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID,
                const uint aType,
                mesh::Face * aFace,
                uint aIndexOnFace,
                const index_t aDofIndexOnField ) :
                Vertex(),
                mTypeID( aType ),
                mMeshBasis( aFace ),
                mIndexOnEntity( aIndexOnFace ),
                mIndexOnField( aDofIndexOnField )
        {
            this->set_id( aID );
            this->set_owner( aFace->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof(
            const id_t aID,
            const uint aType,
            mesh::Element * aElement,
            const uint aIndexOnElement,
            const index_t aDofIndexOnField ) :
                Vertex(),
                mTypeID( aType ),
                mMeshBasis( aElement ),
                mIndexOnEntity( aIndexOnElement ),
                mIndexOnField( aDofIndexOnField )
        {
            this->set_id( aID );
            this->set_owner( aElement->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID,
                    const uint aType,
                mesh::Facet * aFacet ,
                const uint aIndexOnFacet,
                const index_t aDofIndexOnField ) :
                Vertex(),
                mTypeID( aType ),
                mMeshBasis( aFacet ),
                mIndexOnEntity( aIndexOnFacet ),
                mIndexOnField( aDofIndexOnField )

        {
            this->set_id( aID );
            this->set_owner( aFacet->owner() );
        }

//------------------------------------------------------------------------------

        Dof::~Dof()
        {
            this->reset_sources();
        }

//------------------------------------------------------------------------------

        void
        Dof::set_sources( Cell< Dof * > & aSources, Vector< real > & aWeights )
        {
            BELFEM_ASSERT( mNumberOfSources == 0 , "Sources have already been allocated for dof %lu",
                           ( long unsigned int ) this->id() );

            BELFEM_ASSERT( aSources.size() == aWeights.length() ,
                           "Lengths of source and weight vectors for dof %lu do not match (%u vs %u).",
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) aSources.size(),
                           ( unsigned int ) aWeights.length() );

            mNumberOfSources = aSources.size() ;
            mSources = ( Dof ** ) malloc( mNumberOfSources * sizeof( Dof * ) );
            mCoefficients =  ( real * ) malloc( mNumberOfSources * sizeof( real ) );

            for( uint k=0; k<mNumberOfSources; ++k )
            {
                BELFEM_ASSERT( aSources( k ) != nullptr , "Source dof %u for dof %lu is nullptr",
                    ( unsigned int ) k, ( long unsigned int ) this->id() );
                BELFEM_ASSERT( ! aSources( k )->is_hanging() , "Source dofs must not hang!");

                mSources[ k ] = aSources( k );
            }
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mCoefficients[ k ] = aWeights( k );
            }
        }

        void
        Dof::set_sources( Cell< Dof * > & aSources, Cell< real > & aWeights )
        {
            BELFEM_ASSERT( mNumberOfSources == 0 , "Sources have already been allocated for dof %lu",
                           ( long unsigned int ) this->id() );

            BELFEM_ASSERT( aSources.size() == aWeights.size() ,
                           "Lengths of source and weight vectors for dof %lu do not match (%u vs %u).",
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) aSources.size(),
                           ( unsigned int ) aWeights.size() );

            mNumberOfSources = aSources.size() ;
            mSources = ( Dof ** ) malloc( mNumberOfSources * sizeof( Dof * ) );
            mCoefficients =  ( real * ) malloc( mNumberOfSources * sizeof( real ) );

            for( uint k=0; k<mNumberOfSources; ++k )
            {
                BELFEM_ASSERT( aSources( k ) != nullptr , "Source dof %u for dof %lu is nullptr",
                    ( unsigned int ) k, ( long unsigned int ) this->id() );
                BELFEM_ASSERT( ! aSources( k )->is_hanging() , "Source dofs must not hang!");

                mSources[ k ] = aSources( k );
            }
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mCoefficients[ k ] = aWeights( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Dof::set_source( Dof * aSource, const real aWeight)
        {
            BELFEM_ASSERT( mNumberOfSources == 0 , "Sources have already been allocated for dof %lu",
                           ( long unsigned int ) this->id() );

            mNumberOfSources = 1 ;
            mSources = ( Dof ** ) malloc( mNumberOfSources * sizeof( Dof * ) );
            mCoefficients =  ( real * ) malloc( mNumberOfSources * sizeof( real ) );

            mSources[ 0 ] = aSource;
            mCoefficients[ 0 ] = aWeight;
        }

        void
        Dof::reset_sources()
        {
            if ( mNumberOfSources == 0 ) return ;

            ::free( mSources );
            ::free( mCoefficients );
            mNumberOfSources = 0;
            mSources = nullptr;
            mCoefficients = nullptr;
        }

//------------------------------------------------------------------------------

        size_t
        Dof::memory() const
        {
            size_t aMem = sizeof( Dof );

            aMem += this->number_of_vertices() * sizeof( Vertex * );
            aMem += this->number_of_sources()  * sizeof( Dof * );
            aMem += this->number_of_sources()  * sizeof( real );

            return aMem;
        }

//------------------------------------------------------------------------------
    }
}
