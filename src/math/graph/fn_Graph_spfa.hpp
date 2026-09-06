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

#ifndef BELFEM_FN_GRAPH_SPFA_HPP
#define BELFEM_FN_GRAPH_SPFA_HPP

#include <cstdint>

#include "typedefs.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
         * Feasibility solver for a system of difference constraints
         *
         *     theta( head( a ) ) - theta( tail( a ) ) <= weight( a )
         *
         * for every arc a, using the queue-based Bellman-Ford-Moore
         * algorithm (SPFA), warm-started from a BFS spanning forest
         * ( theta( head ) = theta( tail ) + w along tree arcs ) so that
         * only off-tree violations need relaxing. Infeasibility is
         * detected by amortized periodic scans of the predecessor graph;
         * only a strictly negative predecessor cycle is accepted as proof.
         *
         * Vertices are addressed by index in [ 0, aNumVertices );
         * vertices that appear in no arc are left untouched at theta = 0.
         * The returned theta is A feasible solution ( difference
         * constraints are shift-invariant ), not distances from a source.
         *
         * @param aNumVertices    size of the vertex index space
         * @param aArcTails       tail vertex of each arc
         * @param aArcHeads       head vertex of each arc
         * @param aArcWeights     signed weight of each arc
         * @param aTheta          output: feasible potentials, sized to
         *                        aNumVertices ( valid only if return is true )
         * @param aNegativeCycle  output: if infeasible, the arc indices of
         *                        one negative cycle, in cycle order
         *
         * @return true if the system is feasible, false if a negative
         *         cycle exists ( certificate in aNegativeCycle )
         */
        bool
        spfa_difference_constraints(
                const index_t           aNumVertices,
                const Cell< index_t > & aArcTails,
                const Cell< index_t > & aArcHeads,
                const Cell< int64_t > & aArcWeights,
                Cell< int64_t >       & aTheta,
                Cell< index_t >       & aNegativeCycle );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GRAPH_SPFA_HPP
