//
// Created by Christian Messe on 26.10.19.
//

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
        class Field;

        class Block : public Group
        {
            // block on mesh
            mesh::Block *  mBlock = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // creates an empty block
            Block( DofManagerBase * aParent );

//------------------------------------------------------------------------------

            Block(  DofManagerBase * aParent,
                    mesh::Block    * aBlock,
                    const index_t  aNumberOfElements );

//------------------------------------------------------------------------------

            ~Block();

//------------------------------------------------------------------------------

            /**
             * expose block object on mesh
             */
            mesh::Block *
            block();

//------------------------------------------------------------------------------

            ElementType
            element_type() const;

//------------------------------------------------------------------------------

            void
            set_integration_order( const uint aOrder );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            initialize_elements();

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
