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

#include "cl_Input_Section.hpp"
#include "stringtools.hpp"
#include "fn_check_unit.hpp"
#include "petsctools.hpp"
#include "fn_unique.hpp"
namespace belfem
{
    namespace input
    {
        Section::Section( const Section *           aParent,
                 const Cell< string > & aBuffer,
                 const string         & aKey,
                 const index_t          aStartFlag,
                 const index_t          aEndFlag
                 ) :
                 mParent( aParent ),
                 mLevel( aParent == nullptr ? 0 : aParent->level() + 1 ),
                 mBuffer( aBuffer ),
                 mType( this->get_type( aKey ) ) ,
                 mLabel( this->get_label( aKey ) ),
                 mKey( mLabel == "" ? mType : mType + ":" + mLabel),
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
                        // The header is the line immediately BEFORE the
                        // brace, so it must lie inside this section's own scan
                        // range -- that is the invariant, and mStartFlag is
                        // where the range begins.
                        //
                        // At the root mStartFlag is 0, and k - 1 then wraps:
                        // index_t is unsigned, so mBuffer( k - 1 ) reads out of
                        // bounds -- a debug assert, and undefined behaviour in
                        // release, where Cell's bounds check is compiled out.
                        // Reachable from ordinary input, because
                        // remove_comments() and tidy_up() drop comment and
                        // blank lines: a deck whose first section header is
                        // lost to an edit begins its buffer with '{'. Observed
                        // when examples/corc_solder lost its "mesh" header.
                        //
                        // Nested, the same read is in bounds but reaches the
                        // PARENT's brace and names the child "{". Not a crash,
                        // so it survived unnoticed; one predicate covers both,
                        // and no shipped deck has a brace without a header.
                        BELFEM_ERROR( k > mStartFlag,
                            "In '%s': an opening brace has no section header on the line above it. "
                            "A section is named on the line before its '{', and a comment does not "
                            "count as that line.",
                            mKey.c_str() );

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
                        if ( tSection->label() == "" )
                        {
                            mSections[ tSection->type() ] = tSection ;
                        }
                        else
                        {
                            mSections[ tSection->type() + ":" + tSection->label() ] = tSection ;
                        }
                    }
                }
                else if ( mBuffer( k ).find( ";" ) < mBuffer( k ).length() && tSectionCount == 0 )
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

            string tKey ;
            string tString ;

            // catch special case if there is no value assigned to the key
            if( tPos > aLine.length() )
            {
                tPos = aLine.find(";") ;
                tKey = string_to_lower( clean_string( aLine.substr( 0, tPos ) ) );
                tString = "true" ;
            }
            else
            {
                tKey = string_to_lower( clean_string( aLine.substr( 0, tPos )));
                tString = aLine.substr( tPos + 1 );
                tString = clean_string( tString.substr( 0, tString.find( ";" )));
            }

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
                        else if ( tUnit == "°F" )
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
        Section::type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------

        const string &
        Section::label() const
        {
            return mLabel ;
        }

//------------------------------------------------------------------------------

        const string &
        Section::key() const
        {
            return mKey ;
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

        void Section::get_value( const index_t aIndex, string & aKey, value & aValue ) const
        {
            auto tPair = mValueKeys.get_entry( aIndex );

            aKey = tPair.first ;
            aValue = tPair.second ;
        }

//------------------------------------------------------------------------------

        string
        Section::get_units( const string & aKey ) const
        {
            string tString = this->get_string(aKey) ;
            string tUnit ;

            Cell< string > tWords = string_to_words( tString );
            if( tWords.size() == 2 )
            {
                tUnit = tWords( 1 );
            }
            else
            {
                tUnit = "" ;
            }

            return tUnit ;
        }

//------------------------------------------------------------------------------

        int
        Section::get_int( const string & aKey ) const
        {
            return round( this->get_real( aKey ) );
        }

//------------------------------------------------------------------------------

        const Section *
        Section::section( const string & aType ) const
        {
            // see section_exists: the parser lowercases section types, so the
            // lookup lowercases its argument to match, as the key lookups do
            const string tType = string_to_lower( aType ) ;

            BELFEM_ERROR( mSections.key_exists( tType ),
                          this->error_section_not_exisits( aType ).c_str() );

            return mSections( tType );
        }

//------------------------------------------------------------------------------

        const Section *
        Section::section( const string & aType, const string & aLabel ) const
        {
            string tKey = string_to_lower( aType ) + ":"
                        + string_to_lower( aLabel ) ;

            BELFEM_ERROR( mSections.key_exists( tKey ),
                         this->error_section_not_exisits(
                             aType + " : " + aLabel ).c_str() );

            return mSections( tKey );
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
            string aTree = mType ;

            for( int k=mLevel; k>0; k-- )
            {
                tSection = tSection->parent();

                aTree = tSection->key() + "->" + aTree ;
            }
            return aTree ;
        }

//------------------------------------------------------------------------------

        string
        Section::error_section_not_exisits( const string & aSection ) const
        {
            string aMessage = mType + "->" + aSection ;

            const Section * tSection = this ;

            for( int k=mLevel; k>1; k-- )
            {
                tSection = tSection->parent();

                aMessage = tSection->key() + "->" + aMessage ;
            }

            return "Section '" + aMessage + "' does not exist " ;
        }

//------------------------------------------------------------------------------

        string
        Section::error_key_not_exists( const string & aKey ) const
        {
            string aMessage = mType ;

            const Section * tSection = this ;

            for( int k=mLevel; k>1; k-- )
            {
                tSection = tSection->parent();

                aMessage = tSection->key() + "->" + aMessage ;
            }

            return "Key '" + aKey + "' in section " + aMessage + " in file "
                    + tSection->parent()->key() + " does not exist " ;
        }

//------------------------------------------------------------------------------

        string
        Section::error_key_not_real( const string & aKey ) const
        {
            string aMessage = mType ;

            const Section * tSection = this ;

            for( int k=mLevel; k>1; k-- )
            {
                tSection = tSection->parent();

                aMessage = tSection->key() + "->" + aMessage ;
            }

            return "Key '" + aKey = "' in section " + aMessage + " in file "
                    + tSection->parent()->key() + " exists but is not a number : %s",
                    mKeys( aKey ).c_str() ;
        }

        string
        Section::get_type( const string & aLine ) const
        {
            return string_to_lower( clean_string(  aLine.substr( 0, aLine.find( ":" ) ) ) );
        }

        string
        Section::get_label( const string & aLine ) const
        {
            return string_to_lower( clean_string( aLine.find( ":" ) < aLine.length() ?
                aLine.substr( aLine.find( ":" ) + 1, aLine.length() ) : "" ) );
        }


//------------------------------------------------------------------------------

        void
        Section::get_ids( const string & aKey, Vector< id_t > & aIDs ) const
        {
            this->ids_from_string( this->get_string( aKey ), aIDs );
        }

//------------------------------------------------------------------------------

        void
        Section::get_intersection_list( const string & aKey, Vector< id_t > & aIDs1, Vector< id_t > & aIDs2 ) const
        {
            std::string tFullString = this->get_string( aKey ) ;
            index_t tAt = tFullString.find( "@" ) ;

            std::string tLeftString = tFullString.substr( 0, tAt ) ;
            std::string tRightString = tFullString.substr( tAt + 1, tFullString.length() ) ;

            this->ids_from_string( tLeftString, aIDs1 );
            this->ids_from_string( tRightString, aIDs2 ) ;
        }

//------------------------------------------------------------------------------

        std::pair< id_t, id_t >
        Section::get_intersection_ids( const string & aKey ) const
        {
            std::pair< id_t, id_t > aPair ;
            std::string tFullString = this->get_string( aKey ) ;
            index_t tAt = tFullString.find( "@" ) ;

            aPair.first = std::stoi(tFullString.substr( 0, tAt )) ;
            aPair.second = std::stoi(tFullString.substr( tAt + 1, tFullString.length() )) ;

            return aPair ;
        }

//------------------------------------------------------------------------------

        void
        Section::get_ids( const string & aKey, const string & aTape, Vector< id_t > & aIDs )  const
        {
            Cell< string > tWords = string_to_words( string_to_lower( this->get_string( aKey ) ), ',' );

            Cell< id_t > tIDs ;

            for ( string & tWord : tWords )
            {
                if ( tWord.find( "@" ) < tWord.length() )
                {
                    // check for the name
                    if ( tWord.substr( tWord.find( "@" ) + 1, tWord.length() )  == aTape )
                    {
                        Vector< id_t > tIds ;
                        this->ids_from_string( tWord.substr( 0, tWord.find( "@" )), tIds );
                        for ( id_t tId : tIds )
                        {
                            tIDs.push( tId );
                        }
                    }
                }
            }


            // make sure that IDs are unique
            //unique( tIDs );

            // convert to Vector
            aIDs.set_size( tIDs.size() );
            for( id_t k = 0; k<tIDs.size(); ++k )
            {
                aIDs( k ) = tIDs( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Section::ids_from_string( const string & aString, Vector< id_t > & aIDs  ) const
        {
             Cell< string > tWords = string_to_words(
             search_and_replace(search_and_replace(search_and_replace( search_and_replace(
                         aString, " ","" ),","," "),"[", " "),"]", " "));

            Cell< id_t > tGroupIDs;

            for( string & tWord : tWords )
            {
                // let's check if this string contains a colon
                auto tColon = tWord.find( ":" );
                if( tColon != string::npos )
                {
                    string tWordA = tWord.substr( 0, tColon );
                    BELFEM_ERROR( is_integer( tWordA ),
                                  "Error reading id list in section %s : %s",
                                  mType.c_str(), tWordA.c_str() );

                    string tWordB = tWord.substr( tColon+1, tWord.size() );
                    BELFEM_ERROR( is_integer( tWordB ),
                                  "Error reading id list in section %s : %s",
                                  mType.c_str(), tWordB.c_str() );

                    // skip this word if we have an at-sign. These refer to thin shells and will be added later
                    if ( tWordA.find( "@" ) < tWordA.length() || tWordB.find( "@" ) < tWordB.length() )
                    {
                        continue;
                    }

                    id_t tA = std::stoi( tWordA );
                    id_t tB = std::stoi( tWordB );

                    if( tA < tB )
                    {
                        for( id_t tID = tA; tID <= tB; ++tID )
                        {
                            tGroupIDs.push( tID );
                        }
                    }
                    else
                    {
                        for( id_t tID = tB; tID <= tA; ++tID )
                        {
                            tGroupIDs.push( tID );
                        }
                    }
                }
                else
                {
                    string tWordA = tWord.substr( 0, tColon );
                    if ( tWordA.find( "@" ) < tWordA.length()  )
                    {
                        continue;
                    }

                    BELFEM_ERROR( is_integer( tWordA ),
                                  "Error reading id list in section %s : %s",
                                  mType.c_str(), tWordA.c_str() );
                    tGroupIDs.push( std::stoi( tWordA ) );
                }
            }

            //unique( tGroupIDs );

            aIDs.set_size( tGroupIDs.size() );
            for( id_t k = 0; k<tGroupIDs.size(); ++k )
            {
                aIDs( k ) = tGroupIDs( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Section::get_id_groups( const string & aKey, Cell<Cell< id_t >> & aIDs ) const
        {
            const string & tStringFull = this->get_string( aKey );

            //If it is a thin shell, read only what's before the @
            id_t tInd = tStringFull.find( "@" ) ;
            const string tString = tStringFull.substr( 0,tInd ) ;

            Cell< string > tWords = string_to_words(
                    search_and_replace(search_and_replace(search_and_replace( search_and_replace(
                            tString, " ","" ),","," "),"[", " [ "),"]", " ] "));

            bool tGroup = false ;
            Cell< id_t > tGroupIDs;
            for( string & tWord : tWords )
            {
                if (tWord == "[") // open the terminal group
                {
                    tGroup = true;
                }
                else if (tWord == "]") // close the terminal group and add the ids from the group
                {
                    //unique( tGroupIDs );
                    aIDs.push(tGroupIDs);
                    tGroupIDs.clear() ;

                    //Reset the group flag to false
                    tGroup = false;
                }
                else
                {
                    // let's check if this string contains a colon
                    auto tColon = tWord.find( ":" );
                    if( tColon != string::npos )
                    {
                        string tWordA = tWord.substr( 0, tColon );
                        BELFEM_ERROR( is_integer( tWordA ),
                                      "Error reading id list in section %s : %s",
                                      mLabel.c_str(), tWordA.c_str() );

                        string tWordB = tWord.substr( tColon+1, tWord.size() );
                        BELFEM_ERROR( is_integer( tWordB ),
                                      "Error reading id list in section %s : %s",
                                      mLabel.c_str(), tWordB.c_str() );


                        id_t tA = std::stoi( tWordA );
                        id_t tB = std::stoi( tWordB );

                        id_t tFirst = tA < tB ? tA : tB;
                        id_t tLast = tA > tB ? tA : tB;
                        for( id_t tID = tFirst; tID <= tLast; ++tID )
                        {
                            tGroupIDs.push(tID) ;
                            if(!tGroup)
                            {
                                //unique( tGroupIDs );
                                aIDs.push(tGroupIDs);
                                tGroupIDs.clear() ;
                            }
                        }
                    }
                    else
                    {
                        string tWordA = tWord.substr( 0, tColon );
                        BELFEM_ERROR( is_integer( tWordA ),
                                      "Error reading id list in section %s : %s",
                                      mLabel.c_str(), tWordA.c_str() );
                        tGroupIDs.push( std::stoi( tWordA ) );
                        if(!tGroup)
                        {
                            //unique( tGroupIDs );
                            aIDs.push(tGroupIDs);
                            tGroupIDs.clear() ;
                        }
                    }
                }
            }

        }


//------------------------------------------------------------------------------

        void
        Section::get_reals( const string & aKey, Vector< real > & aReals ) const
        {
            const string & tString = this->get_string( aKey );

            Cell< string > tWords = string_to_words(
                    search_and_replace( search_and_replace(
                            tString, " ","" ),","," ") );

            Cell< real > tGroupReals;

            for( string & tWord : tWords )
            {
                string tWordA = tWord.substr( 0, string::npos );
                tGroupReals.push( std::stoi( tWordA ) );
            }


            aReals.set_size( tGroupReals.size() );
            for( id_t k = 0; k<tGroupReals.size(); ++k )
            {
                aReals( k ) = tGroupReals( k );
            }
        }

//------------------------------------------------------------------------------
    }
}
