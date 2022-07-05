//
// Created by Christian Messe on 22.02.22.
//

#ifndef BELFEM_FN_FEM_COMPUTE_BRANDT_HPP
#define BELFEM_FN_FEM_COMPUTE_BRANDT_HPP

#include "cl_FEM_DofManager.hpp"

/**
 * This function computes the Brandt solution
 * see https://doi.org/10.1103/PhysRevB.54.4246
 *
 * It requires a special 2d mesh that looks as follows
 *
 *
 *    8-------------------7
 *    |                   |
 *    |    4---------3    |
 *    |    |  SC( 1) |    |
 *    |    1---------2    |
 *    | AIR ( block 2 )   |
 *    5-------------------6
 */

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * @param    aField  The dof manager that we use
         * @param    aHaHp   the ratio of the applied field and the penentration field
         * @return   aBa     the value of the applied B-Field
         */
        real
        compute_brandt( DofManager * aField, const real aHaHp );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_FEM_COMPUTE_BRANDT_HPP
