//
// Created by christian on 2/19/25.
//

#include "assert.hpp"
#include "cl_Input_Settings.hpp"
namespace belfem
{
    namespace input
    {
        Settings::Settings( const Section * aSection ) :
            mSection( aSection )
        {

        }

        Settings::~Settings()
        {
            // pass
        }


    }
}