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

#ifndef BELFEM_CL_MESH_SOURCEEXPANDER_HPP
#define BELFEM_CL_MESH_SOURCEEXPANDER_HPP

#include "assert.hpp"
#include "cl_Cell.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Mesh.hpp"
#include "cl_Vector.hpp"
namespace belfem
{
    namespace mesh
    {
        class SourceExpander
        {
            Mesh * mMesh ;

            DynamicBitset * mNodeBitset = nullptr ;
            DynamicBitset * mEdgeBitset = nullptr ;
            DynamicBitset * mFaceBitset = nullptr ;
            DynamicBitset * mElementBitset = nullptr ;
            DynamicBitset * mFacetBitset = nullptr ;

            Cell< index_t > mBasisIndices ;

            Cell< Basis * > mSources ;
            Map< const Basis *, real > mWeights ;

            Vector< real > mWork ;

        public:

            SourceExpander( Mesh * aMesh );

            ~SourceExpander();

            void
            run();

            void
            expand_sources( Basis * aBasis );

        private:

            DynamicBitset *
                bitset( const EntityType aEntityType );

            void
            reset_bitsets();

            void
            flag_all_sources( const Basis * aBasis );

            void
            expand_weights( const Basis * aBasis, const real aWeight=1.0 );


        };


        inline DynamicBitset *
        SourceExpander::bitset( const EntityType aEntityType )
        {
            switch( aEntityType )
            {
                case EntityType::NODE :
                {
                    return mNodeBitset ;
                }
                case EntityType::EDGE :
                {
                    return mEdgeBitset ;
                }
                case EntityType::FACE :
                {
                    return mFaceBitset ;
                }
                case EntityType::ELEMENT :
                {
                    return mElementBitset ;
                }
                case EntityType::FACET :
                {
                    return mFacetBitset ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown entity type" );
                    return nullptr ;
                }
            }
        }

        inline
        void SourceExpander::flag_all_sources( const Basis * aBasis )
        {
            for ( uint k=0; k<aBasis->number_of_sources() ; ++k )
            {
                const Basis * tSource = aBasis->source( k ) ;
                BELFEM_ASSERT( aBasis != tSource, "Basis and source are the same" );
                this->bitset( tSource->entity_type() )->set( tSource->index() );
                this->flag_all_sources( tSource );
            }

            if ( aBasis->is_hanging() )
            {
                this->bitset( aBasis->entity_type() )->reset( aBasis->index() );
            }
        }

        inline void
        SourceExpander::expand_weights( const Basis * aBasis, const real aWeight )
        {
            for ( uint k=0; k<aBasis->number_of_sources(); ++k )
            {
                const Basis * tSource = aBasis->source( k ) ;
                if ( tSource->is_hanging() )
                {
                    this->expand_weights( tSource, aWeight * aBasis->weight( k ) );
                }
                else
                {
                    mWeights( tSource ) +=  aWeight * aBasis->weight( k )  ;
                }
            }
        }

    }
}
#endif // BELFEM_CL_MESH_SOURCEEXPANDER_HPP