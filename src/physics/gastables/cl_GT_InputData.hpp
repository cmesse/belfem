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
        /**
         * Reader for gasdata.inp, the table of critical point and molecular
         * data.
         *
         * The file is fixed column, not whitespace delimited, so the field
         * offsets in the implementation are part of the format contract and
         * must move together with the file. Records are indexed by species
         * label into mMap; the map keeps the last record of a repeated label,
         * which the shipped table avoids by carrying exactly one row per label.
         *
         * Reads into GasData and converts to SI on the way: molar mass g/mol to
         * kg/mol, critical pressure bar to Pa. The dipole moment stays in
         * debye, which is what the correlations that use it expect.
         */
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
