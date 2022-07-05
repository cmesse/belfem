//
// Created by christian on 10/25/21.
//

#include "cl_FEM_Group.hpp"
#include "cl_FEM_Cut.hpp"
#include "en_FEM_DomainType.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Cut::Cut(
                DofManager * aParent,
                const id_t aID,
                Cell< mesh::Facet * > & aFacets ) :
                SideSet(  aParent, aID, aFacets, GroupType::CUT )
        {
            this->set_domain_type( DomainType::Cut );
        }

//------------------------------------------------------------------------------

        Cut::~Cut()
        {
        }

//------------------------------------------------------------------------------
    }
}