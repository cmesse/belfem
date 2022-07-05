//
// Created by christian on 7/9/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_BLOCKDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_BLOCKDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_Block.hpp"
namespace belfem
{
    class Mesh ;

    namespace fem
    {
        class Kernel ;
        class DofManager ;

        namespace dofmgr
        {
            class BlockData
            {
                //! the parent object
                DofManager * mParent ;

                //! the kernel
                Kernel * mKernel ;

                //! the mesh this problem runs on
                Mesh * mMesh ;

                // my rank
                const proc_t mMyRank ;

                // map for fem blocks
                Map< id_t, Block * > mBlockMap ;

                // container for blocks
                Cell< Block * > mBlocks ;

                // the empty block is a dummy that is exposed if
                // a set is accessed that doesn't exist on this proc
                Block * mEmptyBlock = nullptr ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                BlockData( DofManager * aParent ) ;

//------------------------------------------------------------------------------

                ~BlockData();
//------------------------------------------------------------------------------

                void
                create_blocks();

//------------------------------------------------------------------------------

                /**
                 * Access a block by its ID. If it doesn't exist on the proc,\
                 * an empty block is returned */
                 Block *
                 block( const id_t aID );

//------------------------------------------------------------------------------

                /**
                 * get the block container
                 **/
                Cell< Block * > &
                blocks();

//------------------------------------------------------------------------------

                /**
                 * tells if the block exists on current mesh
                 */
                bool
                block_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------
                void
                reset();

//------------------------------------------------------------------------------
            };

//------------------------------------------------------------------------------

            inline Block *
            BlockData::block( const id_t aID )
            {
                if ( mBlockMap.key_exists( aID ) )
                {
                    return mBlockMap( aID );
                }
                else
                {
                    mEmptyBlock->set_mesh_id( aID );
                    return mEmptyBlock ;
                }
            }

//------------------------------------------------------------------------------

            inline Cell< Block * > &
            BlockData::blocks()
            {
                return mBlocks ;
            }

//------------------------------------------------------------------------------

            inline bool
            BlockData::block_exists( const id_t aID ) const
            {
                return mBlockMap.key_exists( aID ) ;
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_BLOCKDATA_HPP
