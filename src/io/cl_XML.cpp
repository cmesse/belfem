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

#include "cl_XML.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "stringtools.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    XML::XML( const string & aPath,
              const FileMode aMode ) : mPath( aPath )
    {
#ifdef BELFEM_XML
        // load XML file
        mFile.LoadFile( aPath.c_str() );

        BELFEM_ERROR( ! mFile.Error(),
                "Something went wrong while trying to load file %s : \n Error Core %d \n Message: %s",
                aPath.c_str(),
                mFile.ErrorID(),
                mFile.ErrorIDToName( mFile.ErrorID() ) );

        // set the active element to current file
        mActiveElement = nullptr;

        // set level to zero
        mLevel = 0;
#else
        BELFEM_ERROR( false, "XML support not compiled in" );
#endif
    }

//------------------------------------------------------------------------------

    XML::~XML()
    {
        // reset element

        //mActiveElement = nullptr;
    }

//------------------------------------------------------------------------------

    bool
    XML::child_exists( const string & aLabel )
    {
#ifdef BELFEM_XML
        return ( mActiveElement->FirstChildElement( aLabel.c_str() ) != NULL );
#else
        return false;
#endif
    }

//------------------------------------------------------------------------------

    uint
    XML::number_of_children()
    {
#ifdef BELFEM_XML
        tinyxml2::XMLNode * tNode = mActiveElement->FirstChild();

        if( tNode != NULL )
        {
            uint aCount = 0;

            while( tNode )
            {
                ++aCount;
                tNode = tNode->NextSibling();
            }

            return aCount;
        }
        else
        {
            return 0;
        }
#else
        return 0;
#endif
    }

//-----------------------------------------------------------------------------

    uint
    XML::number_of_children( const string & aLabel )
    {
#ifdef BELFEM_XML
        tinyxml2::XMLElement * tElement;

        if ( mLevel == 0 )
        {
            tElement = mFile.FirstChildElement( aLabel.c_str());
        }
        else
        {
            tElement = mActiveElement->FirstChildElement( aLabel.c_str());
        }

        if( tElement != NULL )
        {
            uint aCount = 0;
            while( tElement )
            {
                if( string(tElement->Name()) == aLabel )
                {
                    ++aCount;
                }
                tElement = tElement->NextSiblingElement( aLabel.c_str() );
            }

            return aCount;
        }
        else
        {
            return 0;
        }
#else
        return 0;
#endif
    }

//-----------------------------------------------------------------------------

    void
    XML::select_first_child( const string & aLabel )
    {
#ifdef BELFEM_XML

        BELFEM_ERROR( this->number_of_children( aLabel ) > 0,
                "Could not find any child named %s in current tree",
                 aLabel.c_str() );


        // select child element
        if( mLevel == 0 )
        {
            mActiveElement = mFile.FirstChildElement( aLabel.c_str() );
        }
        else
        {
            mActiveElement = mActiveElement->FirstChildElement( aLabel.c_str() );
        }

        // increment level
        ++mLevel;
#endif
    }

//-----------------------------------------------------------------------------

    bool
    XML::next_sibling_of_same_name()
    {
#ifdef BELFEM_XML
        tinyxml2::XMLElement * tElement =  mActiveElement->NextSiblingElement(
                mActiveElement->Name() );


        if( tElement != NULL )
        {
            mActiveElement = tElement;
            return true;
        }
        else
        {
            return false;
        }
#else
        return false;
#endif
    }

//-----------------------------------------------------------------------------

    void
    XML::select_parent()
    {
#ifdef BELFEM_XML
        if( mLevel > 0 )
        {
            mActiveElement = mActiveElement->Parent()->ToElement();
            --mLevel;

        }
#endif
    }

//-----------------------------------------------------------------------------

    void
    XML::select_subtree( const string & aTree )
    {
#ifdef BELFEM_XML
        // replace slashes with space
        Cell< string > tWords = string_to_words( search_and_replace( aTree, "/", " " ) );


        mActiveElement = mFile.FirstChildElement( tWords( 0 ).c_str() );

        BELFEM_ERROR( ! mFile.Error(),
                   "Something went wrong while trying to load key %s from file %s: \n Error Core %d \n Message: %s",
                   tWords( 0 ).c_str(),
                   mPath.c_str(),
                   mFile.ErrorID(),
                   mFile.ErrorIDToName( mFile.ErrorID() ) );


        for( size_t k=1; k<tWords.size(); ++k )
        {
            mActiveElement = mActiveElement->FirstChildElement( tWords( k ).c_str() );

            BELFEM_ERROR( ! mFile.Error(),
                       "Something went wrong while trying to load key %s from file %s: \n Error Core %d \n Message: %s",
                       tWords( k ).c_str(),
                       mPath.c_str(),
                       mFile.ErrorID(),
                       mFile.ErrorIDToName( mFile.ErrorID() ) );
        }

        mLevel = tWords.size();
#endif
    }

//-----------------------------------------------------------------------------

    bool
    XML::key_exists( const string & aKey )
    {
#ifdef BELFEM_XML
         return mActiveElement->FirstChildElement( aKey.c_str() ) != NULL;
#else
        return false ;
#endif
    }
//-----------------------------------------------------------------------------
    string
    XML::get_string( const string & aKey )
    {
#ifdef BELFEM_XML
        BELFEM_ERROR( this->key_exists( aKey ), "Key %s does not exist in element %s",
            aKey.c_str(), mActiveElement->Name() );


        string aValue = std::string(
                mActiveElement->FirstChildElement( aKey.c_str() )->GetText() );


        return clean_string( aValue );
#else
        return "";
#endif
    }

//-----------------------------------------------------------------------------

    int
    XML::get_int( const string & aKey )
    {
        return std::stoi( this->get_string( aKey ) );
    }

//-----------------------------------------------------------------------------

    real
    XML::get_real( const string & aKey )
    {
        return std::stod( this->get_string( aKey ) );
    }

//-----------------------------------------------------------------------------

    bool
    XML::get_bool( const string & aKey )
    {
        return string_to_bool( this->get_string( aKey ) );
    }
}