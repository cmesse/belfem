//
// Created by christian on 7/6/21.
//
#include "commtools.hpp"
#include "cl_IwgFactory.hpp"

#include "cl_IWG_StaticHeatConduction.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IwgFactory::IwgFactory( Mesh & aMesh ) :
            mMesh( & aMesh ),
            mNumberOfDimensions( aMesh.number_of_dimensions() )
        {
            this->populate_block_ids() ;
        }

//------------------------------------------------------------------------------

        IwgFactory::IwgFactory( Mesh * aMesh ) :
                mMesh( aMesh ),
                mNumberOfDimensions( aMesh->number_of_dimensions() )
        {
            this->populate_block_ids() ;
        }

//------------------------------------------------------------------------------

        IWG *
        IwgFactory::create_iwg( const IwgType aType, const ModelDimensionality aDimensionality ) const
        {
           BELFEM_ERROR( ! is_maxwell( aType ) , "Can't create Maxwell IWG with IwgFactory. Use MaxwellFactory instead");

            switch ( aType )
            {
                case( IwgType::StaticHeatConduction ) :
                {
                    return new IWG_StaticHeatConduction( aDimensionality );
                }
                default:
                {
                    BELFEM_ERROR( false, "invalid type");
                    return nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IwgFactory::populate_block_ids()
        {
            if( comm_rank() == 0 )
            {
                // grab blocks from mesh
                Cell< mesh::Block * > & tBlocks = mMesh->blocks() ;
                mAllBlockIDs.set_size( tBlocks.size() );

                index_t tCount = 0 ;
                for( mesh::Block * tBlock : tBlocks )
                {
                    mAllBlockIDs( tCount++ ) = tBlock->id() ;
                }

                // check if we run in parallel mode
                proc_t tNumProcs = comm_size() ;
                if( tNumProcs > 1 )
                {
                    // reset counter
                    tCount = 0 ;
                    Vector< proc_t > tCommlist( tNumProcs - 1 );

                    for( proc_t p=1; p<tNumProcs; ++p )
                    {
                        tCommlist( tCount++ ) = p ;
                    }

                    // communicate block IDs
                    send_same( tCommlist, mAllBlockIDs );
                }
            }
            else
            {
                receive( 0, mAllBlockIDs );
            }

            comm_barrier() ;
        }

//------------------------------------------------------------------------------
    }
}