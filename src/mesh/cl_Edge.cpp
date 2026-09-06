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

#include "cl_Edge.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Edge::Edge() :
                Vertex()
        {
        }

//------------------------------------------------------------------------------

        Edge::~Edge()
        {
            this->delete_containers();
        }

//-----------------------------------------------------------------------------
    }
}