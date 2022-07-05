//
// Created by Christian Messe on 28.07.20.
//

#ifndef BELFEM_EN_FEM_BOUNDARYCONDITIONIMPOSING
#define BELFEM_EN_FEM_BOUNDARYCONDITIONIMPOSING

namespace belfem
{
    namespace fem
    {
        enum class BoundaryConditionImposing
        {
            Free      = 0,
            Dirichlet = 1,     // nodal value such as displacement or temperature
            Neumann   = 2,     // nodal flux, such as heat load or force
            Alpha     = 3,     // special type, for convective heat flux, no RADIATION !
            Lambda    = 4,     // special type for maxwell
            Weak      = 5,     // special type for maxwell
            UNDEFINED = 6
        };
    }
}
#endif // BELFEM_EN_FEM_BOUNDARYCONDITIONIMPOSING
