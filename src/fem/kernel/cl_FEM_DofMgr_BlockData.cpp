//
// Created by christian on 7/9/21.
//

#include "commtools.hpp"
#include "cl_FEM_DofMgr_BlockData.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//------------------------------------------------------------------------------

            BlockData::BlockData(  DofManager * aParent ) :
                mParent( aParent ),
                mKernel( aParent->parent() ),
                mMesh( aParent->parent()->mesh() ),
                mMyRank( aParent->rank()  )
            {

            }

//------------------------------------------------------------------------------

            BlockData::~BlockData()
            {
               this->reset() ;
            }

//------------------------------------------------------------------------------

            void
            BlockData::reset()
            {
                if( mEmptyBlock != nullptr )
                {
                    delete mEmptyBlock;
                    mEmptyBlock = nullptr ;
                }


                for( Block * tBlock : mBlocks )
                {
                    delete tBlock ;
                }

                mBlocks.clear() ;
                mBlockMap.clear() ;
            }

//------------------------------------------------------------------------------

            void
            BlockData::
            create_blocks()
            {
                // restore factory settings
                this->reset() ;

                // get the block IDs from the IWG
                const Vector< id_t > & tBlockIDs = mParent->iwg()->selected_blocks() ;

                // get number of blocks
                uint tNumberOfBlocks = tBlockIDs.length();

                // counter for number of elements per block
                Vector< index_t > tElementsPerBlock( tNumberOfBlocks, 0 );

                // loop over all blocks and count owned elements
                for ( uint b = 0; b < tNumberOfBlocks; ++b )
                {
                    if( mMesh->block_exists( tBlockIDs( b ) ) )
                    {
                        Cell< mesh::Element * > & tElements
                                = mMesh->block( tBlockIDs( b ))->elements();

                        for ( mesh::Element * tElement : tElements )
                        {
                            if ( tElement->owner() == mMyRank )
                            {
                                ++tElementsPerBlock( b );
                            }
                        }
                    }
                }

                // count number of blocks that contain at least one element
                index_t tCount = 0;
                for ( uint b = 0; b < tNumberOfBlocks; ++b )
                {
                    if ( tElementsPerBlock( b ) > 0 )
                    {
                        ++tCount;
                    }
                }

                // reset block map
                mBlockMap.clear() ;

                // allocate block container
                mBlocks.set_size( tCount, nullptr );

                // reset counter
                tCount = 0;

                comm_barrier() ;
                // create blocks that contain at least one used element
                for ( uint b = 0; b < tNumberOfBlocks; ++b )
                {
                    if ( tElementsPerBlock( b ) > 0 )
                    {

                        // get mesh block
                        mesh::Block * tBlock = mMesh->block( tBlockIDs( b ) );

                        // add block to container
                        mBlocks( tCount  ) =
                                new fem::Block( mParent, tBlock, tElementsPerBlock( b ) );

                        // add entry to map
                        mBlockMap[ tBlock->id() ] = mBlocks( tCount++ );

                    }
                }

                comm_barrier() ;


                // create the empty block
                mEmptyBlock = new fem::Block( mParent );

            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */