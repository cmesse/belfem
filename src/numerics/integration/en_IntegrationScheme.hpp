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

#ifndef BELFEM_EN_INTEGRATIONSCHEME_HPP
#define BELFEM_EN_INTEGRATIONSCHEME_HPP

#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    enum class IntegrationScheme
    {
        GAUSS,          // best points as found in literature
        GAUSSCLASSIC,   // the classic way of interpolating HEX
        LOBATTO,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    string
    to_string( const IntegrationScheme & aIntegrationScheme );

//------------------------------------------------------------------------------

    IntegrationScheme
    string_to_integration_scheme( const string & aString );

//------------------------------------------------------------------------------
}
#endif //BELFEM_EN_INTEGRATIONSCHEME_HPP
