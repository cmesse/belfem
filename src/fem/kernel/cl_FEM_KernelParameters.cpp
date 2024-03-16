//
// Created by Christian Messe on 29.10.19.
//

#include "cl_FEM_KernelParameters.hpp"
#include "commtools.hpp"
#include "cl_Map.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        KernelParameters::KernelParameters( Mesh & aMesh  ) :
                mMyRank( comm_rank() ),
                mMesh( & aMesh )
        {

            this->init_defaults();
        }


//------------------------------------------------------------------------------

        KernelParameters::KernelParameters( Mesh * aMesh  ) :
            mMyRank( comm_rank() ),
            mMesh( aMesh )
        {
            this->init_defaults();
        }
        
//------------------------------------------------------------------------------
        
        void
        KernelParameters::init_defaults()
        {
            mNumberOfProcs = comm_size();

            if( mMesh->master() == mMyRank )
            {
                // allocate the id container
                mSelectedProcs.set_size( mNumberOfProcs );

                // initialize counter
                uint tCount = 0;

                for ( proc_t p = 0; p < mNumberOfProcs; ++p )
                {
                    mSelectedProcs( tCount++ ) = p;
                }

                // select all blocks by default
                uint tNumBlocks = mMesh->number_of_blocks();

                mBlockIndices.set_size( tNumBlocks );

                for ( uint b = 0; b < tNumBlocks; ++b )
                {
                    mBlockIndices( b ) = b;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::select_procs( const Vector < proc_t > & aProcList )
        {
            if( mMesh->master() == mMyRank )
            {
                mSelectedProcs = aProcList;

                mNumberOfProcs = aProcList.length();
            }

            broadcast( mMesh->master(), mNumberOfProcs );
        }

//-----------------------------------------------------------------------------

        void
        KernelParameters::select_blocks( const Vector< id_t > & aBlockIDs )
        {
            if( mMesh->master() == mMyRank )
            {
                // remember block IDs
                mBlockIDs = aBlockIDs ;

                // the map that links the block ids with indices
                Map <id_t, index_t> tMap;

                // get number of blocks from mesh
                index_t tNumberOfBlocks = mMesh->number_of_blocks();

                Cell< mesh::Block * > & tBlocks = mMesh->blocks();

                // loop over all blocks
                for ( index_t b = 0; b < tNumberOfBlocks; ++b )
                {
                    tMap[ tBlocks( b )->id() ] = b;
                }

                tNumberOfBlocks = aBlockIDs.length();

                // reserve memory for block indices
                mBlockIndices.set_size( tNumberOfBlocks );

                // populate block indices
                for ( index_t b = 0; b < tNumberOfBlocks; ++b )
                {
                    mBlockIndices( b ) = tMap( aBlockIDs( b ) );
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        KernelParameters::set_field_dimensions( const uint aNumDofsPerNode,
                                                const uint aNumDofsPerEdge )
        {
            mNumDofsPerNode.set_size( 1, aNumDofsPerNode );
            mNumDofsPerEdge.set_size( 1, aNumDofsPerEdge );

            this->set_field_dimensions( mNumDofsPerNode, mNumDofsPerEdge );
        }

//-----------------------------------------------------------------------------

        void
        KernelParameters::set_field_dimensions( const Vector< uint > & aNumDofsPerNode )
        {
            mNumDofsPerEdge.set_size( aNumDofsPerNode.length(), 0 );

            this->set_field_dimensions( aNumDofsPerNode, mNumDofsPerEdge );
        }


//-----------------------------------------------------------------------------

        void
        KernelParameters::set_field_dimensions(
                const Vector< uint > & aNumDofsPerNode,
                const Vector< uint > & aNumDofsPerEdge )
        {

            mNumDofsPerNode = aNumDofsPerNode ;
            mNumDofsPerEdge = aNumDofsPerEdge ;

            // check input
            BELFEM_ERROR( aNumDofsPerNode.length() == aNumDofsPerEdge.length(),
                "Length of input vectors does not match" );

            uint tN = mNumDofsPerNode.length();

            // Adjust integration order list
            uint tM = mSideSetIntegrationOrders.length();
            if( tM < tN )
            {
                // backup data
                Vector< uint > tOrders = mSideSetIntegrationOrders;

                // initialize data and set integration order to auto
                mSideSetIntegrationOrders.set_size( tN, 0 );

                // copy data into new vector
                for( uint k=0; k<tM; ++k )
                {
                    mSideSetIntegrationOrders( k ) = tOrders( k );
                }
            }

            // Adjust integration order list
            tM = mBlockIntegrationOrders.length();
            if( tM < tN )
            {
                // backup data
                Vector< uint > tOrders = mBlockIntegrationOrders;

                // initialize data and set integration order to auto
                mBlockIntegrationOrders.set_size( tN, 0 );

                // copy data into new vector
                for( uint k=0; k<tM; ++k )
                {
                    mBlockIntegrationOrders( k ) = tOrders( k );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_sideset_integration_orders
            ( const Vector< uint > & aSideSetIntegrationOrders )
        {
            mSideSetIntegrationOrders = aSideSetIntegrationOrders;
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_block_integration_orders( const uint & aBlockIntegrationOrder )
        {
            mBlockIntegrationOrders.set_size(
                    mNumDofsPerNode.length(),
                    aBlockIntegrationOrder );
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_sideset_integration_orders( const uint & aSideSetIntegrationOrder )
        {
            mSideSetIntegrationOrders.set_size(
                    mNumDofsPerNode.length(),
                    aSideSetIntegrationOrder );
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_integration_scheme( const IntegrationScheme & aIntegrationScheme )
        {
            mIntegrationScheme = aIntegrationScheme ;
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_block_integration_orders
            ( const Vector< uint > & aBlockIntegrationOrders )
        {
            /*BELFEM_ERROR( aBlockIntegrationOrders.length() == mNumDofsPerNode.length(),
                         "length of element integraition orders must be equal field dimension list( %u and %u )",
                         ( unsigned int ) aBlockIntegrationOrders.length(),
                         ( unsigned int ) mNumDofsPerNode.length() ); */

            mBlockIntegrationOrders = aBlockIntegrationOrders;
        }

//------------------------------------------------------------------------------

        void
        KernelParameters::set_auto_partition( const bool aFlag )
        {
            mAutoPartition = aFlag ;
        }

//------------------------------------------------------------------------------
    }
}
