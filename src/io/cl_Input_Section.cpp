//
// Created by christian on 9/17/21.
//

#include "cl_Input_Section.hpp"
#include "stringtools.hpp"
#include "fn_check_unit.hpp"

namespace belfem
{
    namespace input
    {
        Section::Section( const Section *           aParent,
                 const Cell< string > & aBuffer,
                 const string         & aLabel,
                 const index_t          aStartFlag,
                 const index_t          aEndFlag
                 ) :
                 mParent( aParent ),
                 mLevel( aParent == nullptr ? 0 : aParent->level() + 1 ),
                 mBuffer( aBuffer ),
                 mLabel( aLabel ),
                 mStartFlag( aStartFlag ),
                 mEndFlag( aEndFlag )
        {
            this->create_children();
        }

//------------------------------------------------------------------------------

        Section::~Section()
        {
            for ( Section * tSection : mData )
            {
                delete tSection ;
            }
        }

//------------------------------------------------------------------------------

        void
        Section::create_children()
        {
            string tLabel ;
            index_t tStart = 0 ;
            index_t tSectionCount = 0 ;

            for( index_t k=mStartFlag; k<mEndFlag; ++k )
            {
                char tChar = mBuffer( k ).c_str()[ 0 ];
                if( tChar == 123 )
                {
                    if( tSectionCount == 0 )
                    {
                        tLabel = mBuffer( k - 1 );
                        tStart = k+1 ;
                    }
                    ++tSectionCount ;
                }
                else if ( tChar == 125 )
                {
                    --tSectionCount ;
                    if( tSectionCount == 0 )
                    {
                        // create a new section
                        Section * tSection = new Section( this, mBuffer, tLabel, tStart, k );

                        // add section to data container
                        mData.push( tSection );

                        // add section to map
                        mSections[ tLabel ] = tSection ;
                    }
                }
                else if ( mBuffer( k ).find( ":" ) < mBuffer( k ).length() && tSectionCount == 0 )
                {
                    this->create_key( mBuffer( k ) );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Section::create_key( const string & aLine )
        {
            size_t tPos = aLine.find(":") ;
            string tKey = string_to_lower( clean_string( aLine.substr( 0, tPos ) ) );
            string tString =  aLine.substr( tPos+1 ) ;
            tString = clean_string( tString.substr( 0, tString.find(";") ) );

            mKeys[ tKey ] = tString ;

            // add key to cell container
            mKeyLabels.push( tKey );


            // check if this could be a numerical value with unit
            Cell< string > tWords = string_to_words( tString );

            real tReal = to_real( tString );

            // check if value is numeric
            if( ! std::isnan( tReal ) )
            {

                // check if value has a unit
                if( tWords.size() == 1 )
                {
                    // assume SI unit
                    mRealKeys[ tKey ] = tReal ;
                    value tValue = unit_to_si("-");
                    tValue.first *= tReal ;
                    mValueKeys[ tKey ] = tValue ;

                }
                else if( tWords.size() == 2 )
                {
                    const string & tUnit = tWords( 1 );

                    // check if second word is not numeric
                    if( std::isnan( to_real( tUnit ) ) )
                    {

                        value tValue = unit_to_si("K");

                        // catch temperature
                        if( tUnit == "C" || tUnit == "°C" )
                        {
                            tReal += 273.15 ;
                        }
                        else if ( tUnit == "F" ||  tUnit == "°F" )
                        {
                            tReal -= 32.0 ;
                            tReal /= 1.8 ;
                            tReal += 273.15 ;
                        }
                        else if ( tUnit == "R" ||  tUnit == "°R" )
                        {
                            tReal /= 1.8 ;
                        }
                        else if ( tUnit != "K" &&  tUnit != "°K" )
                        {
                            tValue = unit_to_si( tUnit );
                        }

                        // convert units
                        tValue.first *= tReal ;
                        tReal = tValue.first ;

                        mRealKeys[ tKey ]  = tReal ;
                        mValueKeys[ tKey ] = tValue ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        const string &
        Section::label() const
        {
            return mLabel ;
        }

//------------------------------------------------------------------------------

        int
        Section::level() const
        {
            return mLevel ;
        }

//------------------------------------------------------------------------------

        bool
        Section::key_exists( const string & aKey ) const
        {
            return mKeys.key_exists( string_to_lower( aKey ) );
        }

//------------------------------------------------------------------------------

        bool
        Section::key_is_real( const string & aKey ) const
        {
            return mRealKeys.key_exists( string_to_lower( aKey ) );
        }

//------------------------------------------------------------------------------

        const string &
        Section::get_string( const string & aKey ) const
        {
            string tKey = string_to_lower( aKey );

            if( ! mKeys.key_exists( tKey ) )
            {
                BELFEM_ERROR( false,
                             this->error_key_not_exists( aKey ).c_str() );

            }

            return mKeys( tKey );
        }

//------------------------------------------------------------------------------

        bool
        Section::get_bool( const string & aKey ) const
        {
            return string_to_bool( this->get_string( aKey ) );
        }

//------------------------------------------------------------------------------

        real
        Section::get_real( const string & aKey ) const
        {
            string tKey = string_to_lower( aKey );

            if( !  mKeys.key_exists( tKey ) )
            {
                BELFEM_ERROR( false,
                             this->error_key_not_exists( aKey ).c_str() );

            }
            else if ( ! mRealKeys.key_exists( tKey ) )
            {
                BELFEM_ERROR( false,
                             this->error_key_not_real( aKey ).c_str() );
            }

            return mRealKeys( tKey );
        }

//------------------------------------------------------------------------------

        value
        Section::get_value( const string & aKey, const string & aUnit ) const
        {
            string tKey = string_to_lower( aKey );

            if( !  mKeys.key_exists( tKey ) )
            {
                BELFEM_ERROR( false,
                             this->error_key_not_exists( aKey ).c_str() );

            }
            else if ( ! mValueKeys.key_exists( tKey ) )
            {
                BELFEM_ERROR( false,
                             this->error_key_not_real( aKey ).c_str() );
            }

            // get key
            value aValue = mValueKeys( tKey );

            // check unit
            BELFEM_ERROR( check_unit( aValue, aUnit ),
                         "Invalid unit for value %s, expect %s or same dimension",
                         aKey.c_str(), aUnit.c_str() );

            return aValue ;
        }

//------------------------------------------------------------------------------

        int
        Section::get_int( const string & aKey ) const
        {
            return round( this->get_real( aKey ) );
        }

//------------------------------------------------------------------------------

        const Section *
        Section::section( const string & aSection ) const
        {
            if( !  mSections.key_exists( aSection ) )
            {
                BELFEM_ERROR( false,
                             this->error_section_not_exisits( aSection ).c_str() );

            }

            return mSections( aSection );
        }

//------------------------------------------------------------------------------


        const Section *
        Section::parent() const
        {
            return mParent ;
        }

//------------------------------------------------------------------------------

        string
        Section::tree() const
        {
            const Section * tSection = this ;
            string aTree = mLabel ;

            for( int k=mLevel; k>0; k-- )
            {
                tSection = tSection->parent();

                aTree = tSection->label() + "->" + aTree ;
            }
            return aTree ;
        }

//------------------------------------------------------------------------------

        string
        Section::error_section_not_exisits( const string & aSection ) const
        {
            string aMessage = mLabel + "->" + aSection ;

            const Section * tSection = this ;

            for( int k=mLevel; k>1; k-- )
            {
                tSection = tSection->parent();

                aMessage = tSection->label() + "->" + aMessage ;
            }

            return "Section '" + aMessage + " in file "
                    + tSection->parent()->label() + " does not exist " ;
        }

//------------------------------------------------------------------------------

        string
        Section::error_key_not_exists( const string & aKey ) const
        {
            string aMessage = mLabel ;

            const Section * tSection = this ;

            for( int k=mLevel; k>1; k-- )
            {
                tSection = tSection->parent();

                aMessage = tSection->label() + "->" + aMessage ;
            }

            return "Key '" + aKey + "' in section " + aMessage + " in file "
                    + tSection->parent()->label() + " does not exist " ;
        }

//------------------------------------------------------------------------------

        string
        Section::error_key_not_real( const string & aKey ) const
        {
            string aMessage = mLabel ;

            const Section * tSection = this ;

            for( int k=mLevel; k>1; k-- )
            {
                tSection = tSection->parent();

                aMessage = tSection->label() + "->" + aMessage ;
            }

            return "Key '" + aKey = "' in section " + aMessage + " in file "
                    + tSection->parent()->label() + " exists but is not a number : %s",
                    mKeys( aKey ).c_str() ;
        }

//------------------------------------------------------------------------------
    }
}