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
#include "cl_FEM_Block.hpp"


#include "meshtools.hpp"
#include "assert.hpp"

#include "fn_IF_initialize_integration_points.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        // creates an empty block
        Block::Block(  DofManagerBase * aParent, const ElementType aType ) :
                Group( aParent,
                       GroupType::BLOCK,
                       aType,
                       0, 0 )
        {

        }

//------------------------------------------------------------------------------

        Block::Block(
            DofManagerBase * aParent,
            mesh::Block    * aBlock,
            const Vector< index_t > & aOwnedElementIndices,
            const Vector< index_t > & aAuraElementIndices ):
            Group( aParent, GroupType::BLOCK, aBlock->element_type(),
                                              aBlock->id(),
                                              aOwnedElementIndices.length() ),
            mBlock( aBlock )
        {

            mCalc->initialize_integration( this->element_type(), aParent->iwg()->interpolation_type() );

            // set order to auto (default)
            this->set_integration_order( 0 );

            Cell< mesh::Element * > & tElements = mParent->mesh()->elements() ;

            DofManager * tParent = reinterpret_cast< DofManager * >( aParent );


            index_t tCount = 0;
            if ( aOwnedElementIndices.length() > 0 )
            {
                mElements.set_size( aOwnedElementIndices.length(), nullptr );
                for ( index_t e : aOwnedElementIndices )
                {
                    mElements( tCount++ ) = new Element(
                            this, tParent, tElements( e ) );
                }
            }

            if ( aAuraElementIndices.length() > 0 )
            {
                tCount = 0;
                mAuraElements.set_size( aAuraElementIndices.length(), nullptr );
                for ( index_t e : aAuraElementIndices )
                {
                    mAuraElements( tCount++ ) =
                        new Element( this,  tElements( e ) );
                }
            }

            this->create_element_map();

            this->set_activation_mode( tParent->iwg()->block_activation_mode( aBlock->domain_type() ) );
        }

//------------------------------------------------------------------------------

        Block::~Block()
        {
            this->delete_pointers();
        }

//------------------------------------------------------------------------------

        void
        Block::set_integration_order( const uint aOrder )
        {
            // set the number of nodes per element
            mNumberOfNodesPerElement = mesh::number_of_nodes( mElementType ) ;

            mIntegrationOrder = aOrder ;

            this->initialize_lookup_tables( aOrder );

            if( mCalc != nullptr )
            {
                mCalc->set_integration_order( aOrder );
            }
        }
 //------------------------------------------------------------------------------

        void
        Block::initialize_lookup_tables( const uint aIntegrationOrder )
        {
            IntegrationScheme tScheme = IntegrationScheme::GAUSSCLASSIC ;

            // using this switch, we make sure that only tri and tet elements are enriched
            // because quads, hexes and pentas, pyramids etc are not implemented
            bool tEnrich  = mesh::geometry_type( mElementType ) == GeometryType::TRI
                         || mesh::geometry_type( mElementType ) == GeometryType::TET ;

            if( mParent != nullptr )
            {
                tScheme = mParent->integration_scheme() ;
                tEnrich = mParent->iwg()->enrich_sidesets() && tEnrich ;
            }
            if( tEnrich )
            {
                uint tNumFacets = mesh::number_of_facets( mElementType );

                for( IntegrationData * tData : mEnrichmentData )
                {
                    delete tData ;
                }

                InterpolationFunctionFactory tFactory ;

                mEnrichmentData.set_size( tNumFacets, nullptr );

                for( uint f=0; f<tNumFacets; ++f )
                {
                    InterpolationFunction * tBubble = tFactory.create_bubble_function( mElementType, f );

                    mEnrichmentData( f ) = new IntegrationData( mElementType, tBubble, true );
                    mEnrichmentData( f )->populate( aIntegrationOrder, tScheme );
                }
            }
        }

//------------------------------------------------------------------------------
    }
}