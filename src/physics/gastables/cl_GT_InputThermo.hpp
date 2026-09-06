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
     * a class that can read the thermo.inp
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
