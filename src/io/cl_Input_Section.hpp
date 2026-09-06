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

#ifndef BELFEM_CL_INPUT_SECTION_HPP
#define BELFEM_CL_INPUT_SECTION_HPP

#include "typedefs.hpp"
#include "stringtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace input
    {

//------------------------------------------------------------------------------
        /**
         * @brief One hierarchical section of a configuration file.
         *
         * @ingroup grp_io
         * @see @ref io_io_usage_guide
         */
        class Section
        {
            const Section * mParent ;
            const int mLevel ;
            const Cell< string > & mBuffer ;
            const string mType ;
            const string mLabel ;
            const string mKey ;
            const index_t mStartFlag ;
            const index_t mEndFlag ;
            Cell < Section * > mData ;


            Map< string, Section * > mSections ;
            Map< string, string > mKeys ;
            Map< string, real >   mRealKeys ;
            Map< string, value >  mValueKeys ;

            Cell< string > mKeyLabels ;

//------------------------------------------------------------------------------
                public:
//------------------------------------------------------------------------------

            Section(
                    const Section        * aParent,
                    const Cell< string > & aBuffer,
                    const string         & aKey,
                    const index_t          aStartFlag,
                    const index_t          aEndFlag
                );

//------------------------------------------------------------------------------

            ~Section() ;

//------------------------------------------------------------------------------

            /**
             * return the type of this Section
             */
            const string &
            type() const;

//------------------------------------------------------------------------------

            /**
             * return the name of this Section
             */
            const string &
            label() const;

//------------------------------------------------------------------------------

            const string &
            key() const;

//------------------------------------------------------------------------------

            /**
             * tell if a section exists
             */
             bool
             section_exists( const string & aType ) const;

//------------------------------------------------------------------------------

            /**
             * tell if a section exists
             */
            bool
            section_exists( const string & aType, const string & aLabel ) const;

//------------------------------------------------------------------------------

            /**
             * check if a key exists
             */
            bool
            key_exists( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * check if a key is real
             */
            bool
            key_is_real( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * get the value of a key
             */
            const string &
            get_string( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * get the bool value of a key
             */
            bool
            get_bool( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * get the real value of a key
             */
            real
            get_real( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * get the real value of a key with unit
             */
            value
            get_value( const string & aKey , const string & aUnit ) const;

//------------------------------------------------------------------------------

            /**
             * get the key/value pair for a certain index
             */
            void
            get_value( const index_t aIndex, string & aKey, value & aValue ) const ;

//------------------------------------------------------------------------------

            /**
             * get the units of a key
             */
            string
            get_units( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * get the int value of a key
             */
            int
            get_int( const string & aKey ) const;

//------------------------------------------------------------------------------

            void
            get_ids( const string & aKey, Vector< id_t > & aIDs ) const ;

//------------------------------------------------------------------------------

            void
            get_intersection_list( const string & aKey, Vector< id_t > & aIDs1, Vector< id_t > & aIDs2 ) const ;

//------------------------------------------------------------------------------

            std::pair< id_t, id_t >
            get_intersection_ids( const string & aKey ) const ;

//------------------------------------------------------------------------------

            void
            get_ids( const string & aKey, const string & aTape, Vector< id_t > & aIDs ) const ;

//------------------------------------------------------------------------------

            void
            get_id_groups( const string & aKey, Cell<Cell< id_t >> & aIDs ) const ;

//------------------------------------------------------------------------------

            void
            get_reals( const string & aKey, Vector< real > & aReals ) const ;

//------------------------------------------------------------------------------

            /**
             * return a subsection by type string (requires empty label)
             */
             const Section *
             section( const string & aType  ) const ;

//------------------------------------------------------------------------------

            /**
             * return a subsection by type and label
             */
            const Section *
            section( const string & aType, const string & aLabel  ) const ;

//------------------------------------------------------------------------------

             /**
              * access a section by index
              */
             const Section *
             section( const index_t aIndex ) const ;

//------------------------------------------------------------------------------

             /**
              * number of sections
              */
             index_t
             num_sections() const ;

//------------------------------------------------------------------------------

            /**
             * number of keys
             */
             index_t
             num_keys() const ;

//------------------------------------------------------------------------------

            /**
             * return the name of this key
             */
             const string &
             key( const index_t aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * the level of this section
             */
             int
             level() const ;

//------------------------------------------------------------------------------

            /**
             * the parent of this section
             */
             const Section *
             parent() const ;

//------------------------------------------------------------------------------

             /**
               * returns the tree of this section
               */
             string
             tree() const ;

//------------------------------------------------------------------------------

             /**
              * start flag in buffer
              */
             index_t
             start() const ;

//------------------------------------------------------------------------------

             /**
              * end flag in buffer
              */
             index_t
             end() const ;

//------------------------------------------------------------------------------

            void
            ids_from_string( const string & aString, Vector< id_t > & aIDs  ) const ;

//------------------------------------------------------------------------------

            const Cell< string > &
            buffer() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_children();

//------------------------------------------------------------------------------

            void
            create_key( const string & aLine );

//------------------------------------------------------------------------------

            string
            error_section_not_exisits( const string & aSection ) const ;

//------------------------------------------------------------------------------

            string
            error_key_not_exists( const string & aKey ) const ;

//------------------------------------------------------------------------------

            string
            error_key_not_real( const string & aKey ) const ;

//------------------------------------------------------------------------------

            string
            get_type( const string & aLine ) const ;

            string
            get_label( const string & aLine ) const ;

//------------------------------------------------------------------------------
        };



//------------------------------------------------------------------------------

        inline index_t
        Section::num_sections() const
        {
            return mData.size() ;
        }

//------------------------------------------------------------------------------

        inline index_t
        Section::num_keys() const
        {
            return mKeyLabels.size() ;
        }

//------------------------------------------------------------------------------

        inline const string &
        Section::key( const index_t aIndex ) const
        {
            return mKeyLabels( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Section *
        Section::section( const index_t aIndex ) const
        {
            return mData( aIndex );
        }

//------------------------------------------------------------------------------

        inline index_t
        Section::start() const
        {
            return mStartFlag ;
        }

//------------------------------------------------------------------------------

        inline index_t
        Section::end() const
        {
            return mEndFlag ;
        }

//------------------------------------------------------------------------------

        inline bool
        Section::section_exists( const string & aSection ) const
        {
            // section types are lowercased when the file is parsed
            // ( get_type ), exactly as key names are, so the lookup must
            // lowercase its argument too. Keys have always done this; the
            // section lookups did not, which made 'Defect' miss a 'defect'
            // that the parser had already folded to lower case
            return mSections.key_exists( string_to_lower( aSection ) );
        }

//------------------------------------------------------------------------------

        inline bool
        Section::section_exists( const string & aType, const string & aLabel ) const
        {
            // both halves are lowercased at parse ( get_type / get_label )
            string tSection = string_to_lower( aType ) + ":"
                            + string_to_lower( aLabel ) ;
            return mSections.key_exists( tSection );
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        Section::buffer() const
        {
            return mBuffer ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_INPUT_SECTION_HPP
