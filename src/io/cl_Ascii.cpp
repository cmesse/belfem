//
// Created by Christian Messe on 2018-12-26.
//

#include <fstream>
#include <iostream>

#include "cl_Ascii.hpp"
#include "assert.hpp"
#include "stringtools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Ascii::Ascii( const string & aPath, const FileMode & aMode ) :
        mMode( aMode )
    {

        // test if path is absolute
        if( aPath.substr( 0,1 ) == "/" )
        {
            mPath = aPath;
        }
        else
        {
            mPath = sprint( "%s/%s", std::getenv( "PWD" ), aPath.c_str() );
        }

        switch ( aMode )
        {
            case( FileMode::OPEN_RDONLY ) :
            {
                this->load_buffer();
                break;
            }
            case( FileMode::NEW ) :
            {
                mBuffer.clear();
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown file mode for ASCII file" );
            }
        }
    }

//------------------------------------------------------------------------------

    Ascii::~Ascii()
    {
        BELFEM_ERROR( ! mChangedSinceLastSave,
                   "File %s was changed but never saved.",
                    mPath.c_str() );

        mBuffer.clear();
    }
//------------------------------------------------------------------------------

    bool
    Ascii::save()
    {
        BELFEM_ERROR( mMode != FileMode::OPEN_RDONLY,
                "File %s can't be saved since it is opened in write protected mode.",
                mPath.c_str() );


        // open file
        std::ofstream tFile( mPath.c_str(),  std::ofstream::trunc );

        if( tFile )
        {
            // save buffer to file
            for( string & tLine : mBuffer )
            {
                tFile << tLine << std::endl;
            }

            tFile.close();
        }
        else
        {
            BELFEM_ERROR( false,
                       "Something went wrong while trying to save %s.",
                       mPath.c_str() );
        }

        mChangedSinceLastSave = false;

        return mChangedSinceLastSave;
    }

//------------------------------------------------------------------------------

    index_t
    Ascii::length() const
    {
        return mBuffer.size();
    }

//------------------------------------------------------------------------------

    string &
    Ascii::line( const index_t aLineNumber )
    {
        return mBuffer( aLineNumber );
    }

//------------------------------------------------------------------------------

    const string &
    Ascii::line( const index_t aLineNumber ) const
    {
        return mBuffer( aLineNumber );
    }

//------------------------------------------------------------------------------

    void
    Ascii::print( const std::string & aLine )
    {
        mBuffer.push( aLine );

        mChangedSinceLastSave = true;
    }

//------------------------------------------------------------------------------

    void
    Ascii::load_buffer()
    {
        // tidy up buffer
        mBuffer.clear();

        // make sure that file exists
        BELFEM_ERROR( file_exists( mPath ),
            "File %s does not exist.",
            mPath.c_str() );

        // open file
        std::ifstream tFile( mPath );

        // test if file can be opened
        if( tFile )
        {
            // temporary container for string
            string tLine;

            while ( std::getline( tFile, tLine ) )
            {
                mBuffer.push( tLine );
            }

            // close file
            tFile.close();
        }
        else
        {
            BELFEM_ERROR( false, "Someting went wrong while opening file\n %s",
                    mPath.c_str() );
        }
    }

//------------------------------------------------------------------------------

}