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

#ifndef BELFEM_CL_XML_HPP
#define BELFEM_CL_XML_HPP

#ifdef BELFEM_XML
#include <tinyxml2.h>
#endif

#include "typedefs.hpp"
#include "filetools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    /**
     * @brief XML document interface.
     *
     * @ingroup grp_io
     * @see @ref io_io_usage_guide
     */
    class XML
    {
        const string          mPath;
#ifdef BELFEM_XML
        tinyxml2::XMLDocument mFile;

        tinyxml2::XMLElement* mActiveElement = nullptr;

        // current level in XML file
        uint mLevel = 0;
#endif

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        XML( const string & aPath,
             const FileMode aMode = FileMode::OPEN_RDONLY );

//------------------------------------------------------------------------------

        ~XML();

//------------------------------------------------------------------------------

        /**
         * return the path of this file
         */
         inline const string &
         path() const
         {
            return mPath ;
         }

//------------------------------------------------------------------------------

        /**
         * test if at least one child of that name exists
         */
         bool
         child_exists( const string & aLabel );

//------------------------------------------------------------------------------

        /**
         * cont the number of children in current tag
         */
        uint
        number_of_children();

//------------------------------------------------------------------------------

        /**
         * cont the number of children with that name
         */
        uint
        number_of_children( const string & aLabel );

//------------------------------------------------------------------------------

        /**
         * move to the next sibling of the same name; returns false
         * ( selection unchanged ) if there is none
         */
         bool
         next_sibling_of_same_name();

//------------------------------------------------------------------------------

        void
        select_first_child( const string & aLabel );

//------------------------------------------------------------------------------

        void
        select_parent();

//------------------------------------------------------------------------------

        /**
         * directly jump down to this element if you know that it exists
         * separate element names using a "/"
         */
        void
        select_subtree( const string & aTree );

//------------------------------------------------------------------------------

        string
        get_string( const string & aKey );

//------------------------------------------------------------------------------

        int
        get_int( const string & aKey );

//------------------------------------------------------------------------------

        real
        get_real( const string & aKey );

//------------------------------------------------------------------------------

        bool
        get_bool( const string & aKey );

//------------------------------------------------------------------------------

        bool
        key_exists( const string & aKey );

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_XML_HPP
