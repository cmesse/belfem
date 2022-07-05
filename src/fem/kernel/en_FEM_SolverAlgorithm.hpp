//
// Created by christian on 8/24/21.
//

#ifndef BELFEM_EN_FEM_SOLVERALGORITHM_HPP
#define BELFEM_EN_FEM_SOLVERALGORITHM_HPP

namespace belfem
{
    namespace fem
    {
        enum class SolverAlgorithm
        {
            Direct        = 0,
            NewtonRaphson = 1,
            Picard        = 2,
            UNDEFINED     = 3
        };
    }
}

#endif //BELFEM_EN_FEM_SOLVERALGORITHM_HPP
