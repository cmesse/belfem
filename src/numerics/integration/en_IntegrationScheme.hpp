//
// Created by Christian Messe on 14.12.20.
//

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
