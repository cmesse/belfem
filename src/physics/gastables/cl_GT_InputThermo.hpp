//
// Created by Christian Messe on 2019-08-18.
//

#ifndef BELFEM_CL_GT_INPUTTHERMO_HPP
#define BELFEM_CL_GT_INPUTTHERMO_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_GT_RefGas.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

    /**
     * a class that can read the thermo.imp
     */
    class InputThermo : public Ascii
    {
        // map connecting gas labes to line in buffer
        Map< string, uint > mMap;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
        * Default Constructor. Loads an input file from a path
        * @param aPath
           */
        InputThermo( const string & aPath );

        ~InputThermo() = default;

//------------------------------------------------------------------------------

        bool
        entry_exists( RefGas * aRefgas );

//------------------------------------------------------------------------------

        void
        read_data( RefGas * aRefgas );

    };

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_INPUTTHERMO_HPP
