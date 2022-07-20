//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_BLOCK_HPP
#define BELFEM_CL_BLOCK_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
        class GmshReader;

        class Block
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            id_t mID;

            index_t mElementCounter = 0;

            Cell< Element * > mElements;

            string mLabel;

            friend GmshReader;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Block( const id_t aID, const index_t aNumElements );

//------------------------------------------------------------------------------

            ~Block() = default;

//------------------------------------------------------------------------------

            inline id_t
            id() const
            {
                return mID;
            }

//------------------------------------------------------------------------------

            void
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            inline Element *
            element( const index_t aIndex )
            {
                BELFEM_ASSERT( aIndex < mElements.size(),
                              "Invalid Element index %lu in block %lu",
                              ( long unsigned int ) aIndex,
                              ( long unsigned int ) mID );

                return mElements( aIndex );
            }

//------------------------------------------------------------------------------

            inline Cell< Element * > &
            elements()
            {
                return mElements;
            }

//------------------------------------------------------------------------------

            inline index_t
            number_of_elements()
            {
                return mElements.size();
            }

//------------------------------------------------------------------------------

            inline string &
            label()
            {
                return mLabel;
            }

//------------------------------------------------------------------------------

            inline ElementType
            element_type()
            {
                BELFEM_ASSERT( mElements.size() > 0, "Block %u seems to be empty.",
                              ( unsigned int ) mID );

                return mElements( 0 )->type();
            }

//------------------------------------------------------------------------------

            void
            flag_elements();

//------------------------------------------------------------------------------

            void
            unflag_elements();

//------------------------------------------------------------------------------

            void
            flag_nodes();

//------------------------------------------------------------------------------

            void
            flag_corner_nodes();

//------------------------------------------------------------------------------

            void
            flag_edges();

//------------------------------------------------------------------------------

            void
            flag_faces();

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_BLOCK_HPP
