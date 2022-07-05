//
// Created by Christian Messe on 2018-12-27.
//

#ifndef BELFEM_CL_XML_HPP
#define BELFEM_CL_XML_HPP

#include <tinyxml2.h>

#include "typedefs.hpp"
#include "filetools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    class XML
    {
        const string          mPath;
        tinyxml2::XMLDocument mFile;

        tinyxml2::XMLElement* mActiveElement = nullptr;

        // current level in XML file
        uint mLevel = 0;

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
         * count number of siblings of same name
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
