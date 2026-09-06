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

#ifndef BELFEM_CL_MESH_TMATRIXREADER_HPP
#define BELFEM_CL_MESH_TMATRIXREADER_HPP

#include "cl_HDF5.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
        void
        read_node_tmatrix_from_hdf5( Mesh * aMesh, const string & aFile );
    }
}

#endif //BELFEM_CL_MESH_TMATRIXREADER_HPP
