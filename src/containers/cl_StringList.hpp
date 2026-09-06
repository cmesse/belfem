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
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    class StringList
    {
        const uint mMemory;

              uint mCount;

        char ** mData;

        // non-copyable, non-movable (owns raw heap memory)
        StringList( const StringList & ) = delete ;
        StringList & operator=( const StringList & ) = delete ;
        StringList( StringList && ) = delete ;
        StringList & operator=( StringList && ) = delete ;

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
        BELFEM_ASSERT( aIndex < mCount,
                      "Error reading stringlist item: index %u out of bounds ( must be < %u )",
                      ( unsigned int ) aIndex,
                      ( unsigned int ) mCount );

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
