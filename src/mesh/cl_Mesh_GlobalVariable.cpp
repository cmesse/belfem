//
// Created by Christian Messe on 09.10.19.
//

#include "cl_Mesh_GlobalVariable.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        GlobalVariable::GlobalVariable(
                const string & aLabel,
                const id_t   & aID,
                const real     aValue ) :
                mLabel( aLabel ),
                mID( aID ),
                mValue( aValue )
        {
        }

//------------------------------------------------------------------------------
    }
}