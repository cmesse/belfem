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
