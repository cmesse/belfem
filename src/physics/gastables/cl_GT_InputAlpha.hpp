//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_CL_GT_INPUTALPHA_HPP
#define BELFEM_CL_GT_INPUTALPHA_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_GT_GasData.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        // reads the alpha coefficients for a cubic gas
        class InputAlpha : public Ascii
        {
            // map to line in buffer
            Map <string, uint> mMap;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            InputAlpha( const string & aPath );

//----------------------------------------------------------------------------

            ~InputAlpha() = default;

//----------------------------------------------------------------------------

            bool
            entry_exists( GasData * aData );

//----------------------------------------------------------------------------

            void
            read_data( GasData * aData );

//----------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_CL_GT_INPUTALPHA_HPP
