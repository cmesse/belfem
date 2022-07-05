//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_CL_GT_INPUTDATA_HPP
#define BELFEM_CL_GT_INPUTDATA_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_GT_RefGas.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        class InputData : public Ascii
        {
            Map <string, uint> mMap;
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            InputData( const string & aPath );

//----------------------------------------------------------------------------

            ~InputData() = default;

//----------------------------------------------------------------------------

            void
            read_data( GasData * aData );

//----------------------------------------------------------------------------

            bool
            entry_exists( GasData * aData );

//----------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_CL_GT_INPUTDATA_HPP
