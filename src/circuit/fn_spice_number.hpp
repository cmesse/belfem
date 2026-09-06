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

#ifndef BELFEM_FN_SPICE_NUMBER_HPP
#define BELFEM_FN_SPICE_NUMBER_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        /**
         * Convert one SPICE number token into its SI value.
         *
         * Implements the ngspice number grammar (ngspice manual, "Ngspice
         * scale factors" table and the "Letters following a number" rule):
         * an integer (12, -44), a floating point number (3.14159), either
         * followed by an integer exponent (1e-14, 2.65e3), optionally
         * followed by a case-insensitive scale factor
         *
         *   T = 1e12   G = 1e9    Meg = 1e6   K = 1e3    mil = 25.4e-6
         *   m = 1e-3   u = 1e-6   n = 1e-9    p = 1e-12  f = 1e-15  a = 1e-18
         *
         * with longest match first (Meg and mil win over m). Letters
         * immediately after the number or the scale factor are ignored, so
         * 10, 10V, 10Volts and 10Hz are the same number and 1k, 1kOhm and
         * 1kHz the same value. SPICE values are dimensionless SI: F is the
         * femto scale factor, never Farad -- a 100 F capacitor is written
         * "100", while "100F" is 100 fF. The same trap exists for currents:
         * A is the atto scale factor, never Ampere -- a 1 A source is
         * written "1", while "1A" is 1e-18.
         *
         * Errors (BELFEM_ERROR): empty token, no leading number ("abc",
         * "inf", "0x10"), a magnitude outside [DBL_MIN, DBL_MAX] other than
         * exact zero, before or after the scale factor is applied ("1e999",
         * "1e-999", "1e300T", and subnormals like "1e-320"), or any
         * non-letter character after the number/scale factor -- including
         * the old-style embedded suffix "2k5" (write "2.5k") and embedded
         * whitespace. Conversion is locale-independent (std::from_chars).
         *
         * @param aToken   one whitespace-free token from a netlist card
         * @return         the value in SI units
         */
        real
        spice_number_to_si( const string & aToken );

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_SPICE_NUMBER_HPP
