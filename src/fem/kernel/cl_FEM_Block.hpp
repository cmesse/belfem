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
#ifndef BELFEM_CL_FEM_BLOCK_HPP
#define BELFEM_CL_FEM_BLOCK_HPP

#include "cl_Mesh.hpp"
#include "cl_Block.hpp"
#include "cl_Cell.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_Cell.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"


#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace fem
    {
        class Block : public Group
        {
            // block on mesh
            mesh::Block *  mBlock = nullptr ;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // creates an empty block
            Block( DofManagerBase * aParent, const ElementType aType = ElementType::EMPTY );

//------------------------------------------------------------------------------

            Block(  DofManagerBase * aParent,
                    mesh::Block    * aBlock,
                    const Vector< index_t > & aOwnedElementIndices,
                    const Vector< index_t > & aAuraElementIndices );

//------------------------------------------------------------------------------

            ~Block() override;

//------------------------------------------------------------------------------

            /**
             * expose block object on mesh
             */
            mesh::Block *
            block();

//------------------------------------------------------------------------------

            void
            set_integration_order( const uint aOrder );

//------------------------------------------------------------------------------

            void
            initialize_lookup_tables( const uint aIntegrationOrder ) override;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline mesh::Block *
        Block::block()
        {
            return mBlock;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEM_BLOCK_HPP
