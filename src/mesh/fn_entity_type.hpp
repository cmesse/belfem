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

#ifndef BELFEM_FN_ENTITY_TYPE_HPP
#define BELFEM_FN_ENTITY_TYPE_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
    /**
     * guesses the type of the entity based on the passed field
     */
    EntityType
    entity_type( const string & aFieldLabel );

    string
    to_string( const EntityType aEntityType );

}

#endif //BELFEM_FN_ENTITY_TYPE_HPP
