//
// Created by Christian Messe on 2019-07-30.
//

#ifndef BELFEM_CL_STRINGLIST_HPP
#define BELFEM_CL_STRINGLIST_HPP

#include "typedefs.hpp"
#include "assert.hpp"

#ifdef BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-private-field"
#endif

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * a help class needed for the Exodus Writer
     */
    class StringList
    {
        const uint mMemory;

              uint mCount;

        char ** mData;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        StringList( const uint aNumberOfStrings );

//------------------------------------------------------------------------------

        ~StringList();

//------------------------------------------------------------------------------

        void
        push( const string & aString );

//------------------------------------------------------------------------------

        const char *
        item( const uint aIndex ) const;

//------------------------------------------------------------------------------

        /**
         * expose the data container
         */
        char **
        data();

//------------------------------------------------------------------------------
    };

    inline const char *
    StringList::item( const uint aIndex ) const
    {
        BELFEM_ASSERT( aIndex < mMemory,
                      "Error reading stringlist item: index %u out of bounds ( must be < %u )",
                      ( unsigned int ) aIndex,
                      ( unsigned int ) mMemory );

        return mData[ aIndex ];
    }

//-----------------------------------------------------------------------------

    inline char **
    StringList::data()
    {
        return mData;
    }

//------------------------------------------------------------------------------
}

#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#endif

#endif //BELFEM_CL_STRINGLIST_HPP
