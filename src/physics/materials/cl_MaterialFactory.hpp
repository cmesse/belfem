//
// Created by Christian Messe on 14.11.19.
//

#ifndef BELFEM_CL_MATERIALFACTORY_HPP
#define BELFEM_CL_MATERIALFACTORY_HPP

#include "en_Materials.hpp"
#include "cl_Material.hpp"

namespace belfem
{
//----------------------------------------------------------------------------
    class MaterialFactory
    {
//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

         MaterialFactory() = default;

//----------------------------------------------------------------------------

        ~MaterialFactory() = default;

//----------------------------------------------------------------------------

        Material *
        create_material( const MaterialType aMaterial );

//----------------------------------------------------------------------------

        Material *
        create_material( const string & aLabel );

//----------------------------------------------------------------------------
    private:
//----------------------------------------------------------------------------

        Material *
        create_isotropic_material( const MaterialType aMaterial );

//----------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_MATERIALFACTORY_HPP
