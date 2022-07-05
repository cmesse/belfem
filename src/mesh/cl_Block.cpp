//
// Created by Christian Messe on 2019-07-28.
//

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

        void
        Block::insert_element( Element * aElement )
        {
            BELFEM_ASSERT( mElementCounter < mElements.size(),
                "Block %lu is full", ( long unsigned int ) mID );

            mElements( mElementCounter++ ) = aElement;
        }

//------------------------------------------------------------------------------

        void
        Block::flag_elements()
        {
            for( Element* tElement : mElements )
            {
                tElement->flag();
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
            for( Element* tElement : mElements )
            {
                tElement->flag_edges();
            }
        }

//------------------------------------------------------------------------------

        void
        Block::flag_faces()
        {
            for( Element* tElement : mElements )
            {
                tElement->flag_faces();
            }
        }

//------------------------------------------------------------------------------
    }
}