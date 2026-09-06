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

#ifndef FN_to_master_orientation_HPP
#define FN_to_master_orientation_HPP
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Face.hpp"
namespace belfem
{
    namespace mesh
    {
        // orders the node from the slave in the orientation of the master
        // this is necessary for duplicate nodes
        void
        to_master_orientation( Facet * aFacet,
            Cell< Node * > & aSlaveOrientation,
            Cell< Node * > & aMasterOrientation );

        void
        to_master_orientation( Face * aFace,
            Cell< Node * > & aSlaveOrientation,
            Cell< Node * > & aMasterOrientation );

        void
        to_master_orientation( Facet * aFacet,
           Cell< Edge * > & aSlaveOrientation,
           Cell< Edge * > & aMasterOrientation );

        void
        to_master_orientation( Face * aFace ,
           Cell< Edge * > & aSlaveOrientation,
           Cell< Edge * > & aMasterOrientation );

    }
}
#endif // FN_to_master_orientation_HPP
