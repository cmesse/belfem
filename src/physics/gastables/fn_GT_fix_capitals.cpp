//
// Created by Christian Messe on 02.09.19.
//
#include "fn_GT_fix_capitals.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        string
        fix_capitals( const string & aLabel )
        {
            uint tN = aLabel.length();
            if ( tN == 0 )
            {
                return aLabel;
            }
            else if ( tN == 1 )
            {
                return string_to_upper( aLabel );
            }
            else
            {
                string aResult = string_to_upper( aLabel.substr( 0, 1 ) );
                aResult += string_to_lower( aLabel.substr( 1, tN - 1 ) );
                return aResult;
            }
        }

//------------------------------------------------------------------------------
    }
}