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

#include "cl_FEM_KernelParameters.hpp"
#include "commtools.hpp"
#include "cl_Map.hpp"
#include "assert.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        KernelParameters::KernelParameters( Mesh & aMesh  ) :
                mCommRank( comm_rank() ),
                mCommSize( comm_size() ),
                mKernel( nullptr ),
                mMesh( & aMesh )
        {

            this->init_defaults();
        }


//------------------------------------------------------------------------------

        KernelParameters::KernelParameters( Mesh * aMesh  ) :
                mCommRank( comm_rank() ),
                mCommSize( comm_size() ),
                mKernel( nullptr ),
                mMesh( aMesh )
        {
            this->init_defaults();
        }

        KernelParameters::KernelParameters( Kernel * aKernel ) :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mKernel( aKernel ),
            mMesh( aKernel->mesh() )
        {
            this->init_defaults();
        }

//------------------------------------------------------------------------------
        
        void
        KernelParameters::init_defaults()
        {
            if( mMesh->master() == mCommRank )
            {
                // select all blocks by default
                uint tNumBlocks = mMesh->number_of_blocks();

                mBlockIndices.set_size( tNumBlocks );

                for ( uint b = 0; b < tNumBlocks; ++b )
                {
                    mBlockIndices( b ) = b;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        KernelParameters::select_blocks( const Vector< id_t > & aBlockIDs )
        {
            if( mMesh->master() == mCommRank )
            {
                // remember block IDs
                mBlockIDs = aBlockIDs ;

                // reserve memory for block indices
                mBlockIndices.set_size( aBlockIDs.length() );

                // populate block indices
                uint tCount = 0 ;
                for ( id_t b : aBlockIDs )
                {
                    mBlockIndices( tCount++ ) = mMesh->block( b )->index();
                }
            }
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::select_blocks( const Cell< id_t > & aBlockIDs )
        {
            if( mMesh->master() == mCommRank )
            {
                // remember block IDs
                mBlockIDs.set_size( aBlockIDs.size() );
                uint tCount = 0 ;
                for ( id_t b : aBlockIDs )
                {
                    mBlockIDs( tCount++ ) = b;
                }

                // reserve memory for block indices
                mBlockIndices.set_size( aBlockIDs.size() );

                // populate block indices
                tCount = 0 ;
                for ( id_t b : aBlockIDs )
                {
                    mBlockIndices( tCount++ ) = mMesh->block( b )->index();
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        KernelParameters::select_sidesets( const Vector< id_t > & aSidesetIDs )
        {
            if( mMesh->master() == mCommRank )
            {
                // remember block IDs
                mSideSetIDs = aSidesetIDs ;

                // reserve memory for block indices
                mSidesetIndices.set_size( aSidesetIDs.length() );

                // populate block indices
                uint tCount = 0 ;
                for ( id_t s : aSidesetIDs )
                {
                    mSidesetIndices( tCount++ ) = mMesh->sideset( s )->index();
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        KernelParameters::select_sidesets( const Cell< id_t > & aSidesetIDs )
        {
            if( mMesh->master() == mCommRank )
            {
                // remember block IDs
                uint tCount = 0 ;

                mSideSetIDs.set_size( aSidesetIDs.size() );
                for ( id_t s : aSidesetIDs )
                {
                    mSideSetIDs( tCount++ ) = s;
                }

                // reserve memory for block indices
                mSidesetIndices.set_size( aSidesetIDs.size() );

                // populate block indices
                tCount = 0 ;
                for ( id_t s : aSidesetIDs )
                {
                    mSidesetIndices( tCount++ ) = mMesh->sideset( s )->index();
                }
            }
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_integration_scheme( const IntegrationScheme & aIntegrationScheme )
        {
            mIntegrationScheme = aIntegrationScheme ;
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_reordering_method( const ReorderingMethod aMethod )
        {
            mReorderingMethod = aMethod ;
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_auto_partition( const bool aFlag )
        {
            mAutoPartition = aFlag ;
        }

//------------------------------------------------------------------------------

        Kernel *
        KernelParameters::kernel()
        {
            return mKernel ;
        }

//------------------------------------------------------------------------------

        Mesh *
        KernelParameters::mesh()
        {
            return mKernel == nullptr ? mMesh :  mKernel->mesh() ;
        }

//------------------------------------------------------------------------------
    }
}
