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

#include "cl_Block.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Block::Block( const id_t aID, const index_t aNumElements ) :
            mID( aID ),
            mElementCounter( 0 )
        {
            mElements.set_size( aNumElements, nullptr );

            mLabel = aID < 10 ?
                     sprint( "Block_0%u", ( unsigned int ) aID ) :
                     sprint( "Block_%u", ( unsigned int ) aID );
        }

//------------------------------------------------------------------------------

        Block::~Block()
        {
            for ( Element * tElement : mElements )
            {
                delete tElement ;
            }
        }

//------------------------------------------------------------------------------

        void
        Block::insert_element( Element * aElement )
        {
            BELFEM_ASSERT( mElementCounter < mElements.size(),
                "Block %lu is full", ( long unsigned int ) mID );

            mElements( mElementCounter++ ) = aElement;
        }

//------------------------------------------------------------------------------

        void
        Block::flag_elements( const uint8_t aFlagIndex )
        {
            for( Element* tElement : mElements )
            {
                tElement->flag( aFlagIndex );
            }
        }

//------------------------------------------------------------------------------

        void
        Block::unflag_elements( const uint8_t aFlagIndex )
        {
            for( Element* tElement : mElements )
            {
                tElement->unflag( aFlagIndex );
            }
        }

//------------------------------------------------------------------------------

        void
        Block::flag_nodes()
        {
            for( Element* tElement : mElements )
            {
                tElement->flag_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        Block::unflag_nodes()
        {
            for( Element* tElement : mElements )
            {
                tElement->unflag_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        Block::flag_corner_nodes()
        {
            for( Element* tElement : mElements )
            {
                tElement->flag_corner_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        Block::flag_edges()
        {
            if ( mElements.size() == 0 ) return;
            if ( ! mElements( 0 )->has_edges() ) return;

            for( Element* tElement : mElements )
            {
                tElement->flag_edges();
            }
        }

//------------------------------------------------------------------------------

        void
        Block::flag_faces()
        {
            if ( mElements.size() == 0 ) return ;
            if ( ! mElements( 0 )->has_faces() ) return ;

            for( Element* tElement : mElements )
            {
                tElement->flag_faces();
            }
        }

//------------------------------------------------------------------------------

        void
        Block::set_edges_flag( const bool aFlag )
        {
            mHasEdges = aFlag;
        }

        void
        Block::set_faces_flag( const bool aFlag )
        {
            mHasFaces = aFlag;
        }

//------------------------------------------------------------------------------
    }
}
