//
// Created by christian on 7/7/21.
//

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
