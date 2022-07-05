//
// Created by christian on 9/17/21.
//

#include "cl_InputFile.hpp"
#include "commtools.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    InputFile::InputFile( const string & aPath ) :
            Ascii( aPath, FileMode::OPEN_RDONLY )
    {
        this->remove_comments() ;
        this->tidy_up() ;

        mData = new input::Section( nullptr,
                                mBuffer,
                                basename( aPath ),
                                0,
                                mBuffer.size() );
    }

//------------------------------------------------------------------------------

    InputFile::~InputFile()
    {
        if( mData != nullptr )
        {
            delete mData ;
        }
    }

//------------------------------------------------------------------------------

    void
    InputFile::print()
    {
        uint tCount = 0 ;
        for ( string & tLine : mBuffer )
        {
            std::cout << tCount++ << ": " << tLine << std::endl ;
        }
    }

//------------------------------------------------------------------------------

    void
    InputFile::remove_comments()
    {
        for ( string & tLine : mBuffer )
        {
            // fix colons
            tLine = search_and_replace( tLine, ":", " : ");

            // fix semicolons
            tLine = search_and_replace( tLine, ";", " ; " );

            // remove comments and clean up
            tLine = clean_string( tLine.substr( 0, tLine.find( "//" , 0 )  ) );
        }
    }

//-----------------------------------------------------------------------------

    void
    InputFile::tidy_up()
    {
        Cell < string > tBuffer ;

        // make sure that all lines are split for "{" and "}"
        string tTag = "{";
        for( uint k=0; k<2; ++k )
        {
            tBuffer.clear() ;

            for ( string & tLine : mBuffer )
            {
                size_t tLength = tLine.length() ;
                size_t tPos = 0 ;
                size_t tFound = tLine.find(tTag, tPos );

                while( tFound < tLength )
                {
                    tBuffer.push( tLine.substr( tPos, tFound - tPos ) );
                    tBuffer.push(tTag);
                    tPos = tFound + 1 ;
                    tFound = tLine.find(tTag, tPos );
                }
                tBuffer.push( tLine.substr( tPos, tFound-tPos ) );
            }
            tTag = "}";
            mBuffer = tBuffer ;

        }

        // handle line breaks for semicolon
        tBuffer.clear() ;

        for ( string & tLine : mBuffer )
        {
            size_t tLength = tLine.length() ;
            size_t tPos = 0 ;
            size_t tFound = tLine.find(tTag, tPos );

            while( tFound < tLength )
            {
                tBuffer.push( tLine.substr( tPos, tFound - tPos+1 ) );
                tPos = tFound + 1 ;
                tFound = tLine.find(tTag, tPos );
            }
            tBuffer.push( tLine.substr( tPos, tFound-tPos ) );
        }
        tTag = "}";
        mBuffer = tBuffer ;

        // tidy up
        mBuffer.clear() ;

        for(  string & tLine : tBuffer )
        {
            string tNewLine = clean_string( tLine ) ;
            if ( tNewLine.size() > 0 )
            {
                char tChar = tNewLine.c_str()[ 0 ];
                if( tChar != 32 && tChar != 9 && tChar != 13 && tChar != 10 )
                {
                    mBuffer.push(  tNewLine );
                }
            }
        }
    }


//-----------------------------------------------------------------------------

    const input::Section *
    InputFile::section( const string & aSection  ) const
    {
        return mData->section( aSection );
    }

//-----------------------------------------------------------------------------

    const input::Section *
    InputFile::section( const index_t aIndex  ) const
    {
        return mData->section( aIndex );
    }

//------------------------------------------------------------------------------

    index_t
    InputFile::num_sections() const
    {
        return mData->num_sections() ;
    }

//------------------------------------------------------------------------------

    bool
    InputFile::section_exists( const string & aSection ) const
    {
        return mData->section_exists( aSection );
    }

//------------------------------------------------------------------------------
}