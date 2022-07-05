//
// Created by christian on 9/17/21.
//

#ifndef BELFEM_CL_INPUT_SECTION_HPP
#define BELFEM_CL_INPUT_SECTION_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"

namespace belfem
{
    namespace input
    {

//------------------------------------------------------------------------------
        class Section
        {
            const Section * mParent ;
            const int mLevel ;
            const Cell< string > & mBuffer ;
            const string mLabel ;
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
                    const string         & aLabel,
                    const index_t          aStartFlag,
                    const index_t          aEndFlag
                );

//------------------------------------------------------------------------------

            ~Section() ;

//------------------------------------------------------------------------------

            /**
             * return the name of this Section
             */
            const string &
            label() const;

//------------------------------------------------------------------------------

            /**
             * tell if a section exists
             */
             bool
             section_exists( const string & aSection ) const;

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
              * get the int value of a key
              */
             int
             get_int( const string & aKey ) const;

//------------------------------------------------------------------------------

            /**
             * return a subsection by label
             */
             const Section *
             section( const string & aSection  ) const ;

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
            return mSections.key_exists( aSection );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_INPUT_SECTION_HPP
