//
// Created by christian on 7/14/21.
//
#include "commtools.hpp"
#include "cl_FEM_DofMgr_BearingData.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Bearing.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//------------------------------------------------------------------------------

            BearingData::BearingData( DofManager * aParent ) :
                    mParent( aParent ),
                    mKernel( aParent->parent() ),
                    mMesh( aParent->parent()->mesh() ),
                    mMyRank( aParent->rank() )
            {

            }

//------------------------------------------------------------------------------

            BearingData::~BearingData()
            {
                this->reset();
            }

//------------------------------------------------------------------------------

            void
            BearingData::reset()
            {
                // delete map
                mBearingMap.clear() ;

                // delete pointers
                for( Bearing * tBearing : mBearings )
                {
                    delete tBearing ;
                }

                mBearings.clear();

                if( mEmptyBearing != nullptr )
                {
                    delete mEmptyBearing ;
                    mEmptyBearing = nullptr ;
                }
            }

//------------------------------------------------------------------------------

            void
            BearingData::create_bearings()
            {
                // restore factory settings
                this->reset() ;

                // initialize counter
                index_t tCount = 0 ;

                // allocate memory for bearings
                mBearings.set_size( mMesh->vertices().size(), nullptr );

                // loop over all vertices on mesh
                for( mesh::Element * tVertex : mMesh->vertices() )
                {
                    mBearings( tCount ) = new Bearing(
                            mParent,
                            tVertex->id(),
                            tVertex->node( 0 ) );

                    // add entry to map
                    mBearingMap[ tVertex->id() ] = mBearings( tCount++ );
                }

                // create the empty bearing
                mEmptyBearing = new Bearing( mParent );
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */