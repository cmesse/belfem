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

#ifndef BELFEM_CL_GT_INPUTTRANSPORT_HPP
#define BELFEM_CL_GT_INPUTTRANSPORT_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_TransportPoly.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        /**
         * Reader for trans.inp, the CEA transport property table.
         *
         * One species may carry several viscosity and conductivity intervals,
         * and the two are interleaved in the file, so the reader hands RefGas
         * both sets separately. The file also holds interaction pairs - records
         * keyed on two species names - used for the binary viscosity of
         * mixtures.
         *
         * Like the other tables this one is fixed column; the temperature range
         * fields in particular are narrow enough that the parse window has to
         * be exact, or a fractional bound such as 63.7 K silently reads as 63.
         */
        class InputTransport: public Ascii
        {
            Map <string, uint> mMap;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * Default Constructor. Loads an input file from a path
             * @param aPath
             */
            InputTransport( const string & aPath );

            ~InputTransport() = default;

//------------------------------------------------------------------------------

            void
            read_data( RefGas * aRefgas );

//------------------------------------------------------------------------------

            uint
            entry_exists( RefGas * aRefgas  );

//------------------------------------------------------------------------------
            /*
             * test if interaction parameter exists
             */
            bool
            interaction_parameter_exists(
                    const string & aLabelA,
                    const string & aLabelB );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            TransportPoly *
            read_polynomial( uint & aLineCount );

//------------------------------------------------------------------------------

            real
            word_to_real( const string & aWord );

        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_INPUTTRANSPORT_HPP
