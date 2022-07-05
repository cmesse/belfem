//
// Created by Christian Messe on 14.12.20.
//

#include "en_IntegrationScheme.hpp"
#include "stringtools.hpp"
#include "assert.hpp"
namespace belfem
{
//------------------------------------------------------------------------------


    string
    to_string( const IntegrationScheme & aIntegrationScheme )
    {
        switch ( aIntegrationScheme )
        {
            case( IntegrationScheme::GAUSSCLASSIC ) :
            {
                return "GaussClassic" ;
            }
            case( IntegrationScheme::GAUSS ) :
            {
                return "GaussModern" ;
            }
            case( IntegrationScheme::LOBATTO ) :
            {
                return "Lobatto" ;
            }
            default:
            {
                return "undefined" ;
            }
        }
    }

//------------------------------------------------------------------------------

    IntegrationScheme
    string_to_integration_scheme( const string & aString )
    {
        string tString = string_to_lower( search_and_replace( search_and_replace(
                search_and_replace( aString,"ß","ss" ), " ", "" ), "_", "" ) );


        if ( tString == "gaussclassic")
        {
            return IntegrationScheme::GAUSSCLASSIC ;
        }
        else if ( tString == "gauss" || tString == "gaussmodern" )
        {
            return IntegrationScheme::GAUSS ;
        }
        else if ( tString == "lobatto" )
        {
            return IntegrationScheme::LOBATTO ;
        }
        else
        {
            BELFEM_ERROR( false, "Unknown Integration Scheme: %s", aString.c_str() );
            return IntegrationScheme::UNDEFINED ;
        }
    }

//------------------------------------------------------------------------------
}
