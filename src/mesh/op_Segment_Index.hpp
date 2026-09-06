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

#ifndef OP_SEGMENT_INDEX_HPP
#define OP_SEGMENT_INDEX_HPP

#include "cl_Segment.hpp"

namespace belfem
{
    namespace mesh
    {
        // comparison object for sorting
        inline struct OpSegmentIndex
        {
            inline bool
            operator()( const Segment * aA, const Segment * aB )
            {
                return aA->index() < aB->index();
            }
        } opSegmentIndex;
    }
}


#endif //OP_SEGMENT_INDEX_HPP
