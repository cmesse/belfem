//
// Created by Christian Messe on 02.09.19.
//
#include "fn_GT_fix_label.hpp"
#include "stringtools.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        string
        fix_label( const string & aLabel )
        {
            if( aLabel == "C2H3,vinyl" )
            {
                return "C2H3";
            }
            else
            {
                return search_and_replace(
                        search_and_replace( search_and_replace( aLabel,
                        "AR", "Ar" ),
                        "AL", "Al" ),
                        "CL", "Cl" );
            }
        }

//------------------------------------------------------------------------------
    }
}