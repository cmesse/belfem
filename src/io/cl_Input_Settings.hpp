//
// Created by christian on 2/19/25.
//

#ifndef CL_INPUT_SETTINGS_HPP
#define CL_INPUT_SETTINGS_HPP

#include "cl_Input_Section.hpp"

namespace belfem
{
    namespace input
    {
        class Settings
        {

        protected:

            const Section * mSection = nullptr ;

//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            Settings( const Section * aSection );

//------------------------------------------------------------------------------

            virtual ~Settings() ;

//------------------------------------------------------------------------------
        };
    }
}
#endif //CL_INPUT_SETTINGS_HPP
