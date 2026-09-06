/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_CREATE_DATABASE_MESH_HPP
#define BELFEM_FN_CREATE_DATABASE_MESH_HPP

#include "constants.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    inline Mesh * create_database_mesh( const uint aOrder=2 )
    {
        return new Mesh( aOrder, { 95, 35, 37 },
                    { 4.0, 0.1, 5 * constant::deg }, { 0.0, -2, 0.0 } );
    }
}
#endif //BELFEM_FN_CREATE_DATABASE_MESH_HPP