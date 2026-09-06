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

#include "commtools.hpp"
#include "cl_IwgFactory.hpp"

#include "cl_IWG_Poisson.hpp"
#include "cl_IWG_StaticHeatConduction.hpp"
#include "cl_IWG_TransientHeatConduction.hpp"

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
                case( IwgType::Poisson ) :
                {
                    return new IWG_Poisson( aDimensionality );
                }

                case( IwgType::StaticHeatConduction ) :
                {
                    return new IWG_StaticHeatConduction( aDimensionality );
                }

                case( IwgType::TransientHeatConduction ) :
                {
                    return new IWG_TransientHeatConduction( aDimensionality );
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

            }

            broadcast( mAllBlockIDs );

            comm_barrier() ;
        }

//------------------------------------------------------------------------------
    }
}
