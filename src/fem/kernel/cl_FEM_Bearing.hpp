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

#ifndef BELFEM_CL_FEM_BEARING_HPP
#define BELFEM_CL_FEM_BEARING_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_FEM_Dof.hpp"

namespace belfem
{
    namespace fem
    {
        class DofManagerBase;
        class Field ;
        class DofManager ;

//------------------------------------------------------------------------------

        class Bearing
        {
            DofManagerBase * mParent    = nullptr;

            const id_t mID;

            mesh::Node * mNode = nullptr;

            uint mNumDofs = 0 ;

            Dof ** mDOFs = nullptr ;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // constructor for empty bearing
            Bearing( DofManagerBase * aParent );

//------------------------------------------------------------------------------

            Bearing( DofManager * aParent, const id_t aID, mesh::Node * aNode );

//------------------------------------------------------------------------------

            ~Bearing();

//------------------------------------------------------------------------------

            id_t
            id() const;

//------------------------------------------------------------------------------

            mesh::Node *
            node();

//------------------------------------------------------------------------------

            void
            impose_dirichlet( const real aValue, const uint aDofType=0 );

//------------------------------------------------------------------------------

            void
            allocate_dof_container( const index_t aNumDofs );

//------------------------------------------------------------------------------

            void
            insert_dof( Dof * aDof, const index_t aIndex );

//------------------------------------------------------------------------------

            void
            free();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline id_t
        Bearing::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline mesh::Node *
        Bearing::node()
        {
            return mNode;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_BEARING_HPP
