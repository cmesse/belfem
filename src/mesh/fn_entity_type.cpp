//
// Created by christian on 7/7/21.
//

#include "fn_entity_type.hpp"
#include "stringtools.hpp"

namespace belfem
{
    /**
     * guesses the type of the entity based on the passed field
     */
    EntityType
    entity_type( const string & aFieldLabel )
    {
        const string tFieldLabel = string_to_lower( aFieldLabel );

        // check if this could be an edge field
        if( aFieldLabel.length() > 4 )
        {
            if( tFieldLabel.substr( 0, 4 ) == "edge" )
            {
                return EntityType::EDGE ;
            }
            else if( tFieldLabel.substr( 0, 4 ) == "face" )
            {
                return EntityType::FACE ;
            }
            else if( tFieldLabel.substr( 0, 4 ) == "cell" )
            {
                return EntityType::CELL ;
            }
            else if ( aFieldLabel.length() > 5 )
            {
                if( tFieldLabel.substr( 0, 6 ) == "lambda" )
                {
                    return EntityType::FACET ;
                }
                else if ( aFieldLabel.length() > 6 )
                {
                    if( tFieldLabel.substr( 0, 7 ) == "element" )
                    {
                        return EntityType::ELEMENT ;
                    }
                    else
                    {
                        return EntityType::NODE ;
                    }
                }
                else
                {
                    return EntityType::NODE ;
                }
            }
            else
            {
                return EntityType::NODE ;
            }
        }
        else
        {
            return EntityType::NODE ;
        }
    }

    string
    to_string( const EntityType aEntityType )
    {
        switch( aEntityType )
        {
            case( EntityType::NODE ) :
            {
                return "node" ;
            }
            case( EntityType::FACET ) :
            {
                return "facet" ;
            }
            case( EntityType::ELEMENT ) :
            {
                return "element" ;
            }
            case( EntityType::EDGE ) :
            {
                return "edge" ;
            }
            case( EntityType::FACE ) :
            {
                return "face" ;
            }
            case( EntityType::CELL ) :
            {
                return "cell" ;
            }
            default:
            {
                return "undefined" ;
            }
        }
    }
}
