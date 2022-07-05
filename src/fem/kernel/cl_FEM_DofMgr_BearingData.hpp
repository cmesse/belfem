//
// Created by christian on 7/14/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_BEARINGDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_BEARINGDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_IWG.hpp"

namespace belfem
{
    class Mesh;

    namespace fem
    {
        class Bearing;

        class Kernel;

        class DofManager;

        namespace dofmgr
        {
            class BearingData
            {
                //! the parent object
                DofManager * mParent;

                //! the kernel
                Kernel * mKernel;

                //! the mesh this problem runs on
                Mesh * mMesh;

                // my rank
                const proc_t mMyRank;

                Cell< Bearing * > mBearings;

                Bearing * mEmptyBearing = nullptr;

                Map< id_t, Bearing * > mBearingMap;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                BearingData( DofManager * aParent );

//------------------------------------------------------------------------------

                ~BearingData();

//------------------------------------------------------------------------------

                //! called by dof manager
                void
                create_bearings();

//------------------------------------------------------------------------------

                /**
                 * return one specific bearing
                 */
                 Bearing *
                 bearing( const id_t aID );

//------------------------------------------------------------------------------

                void
                reset();

//------------------------------------------------------------------------------
            };

//------------------------------------------------------------------------------

            inline Bearing *
            BearingData::bearing( const id_t aID )
            {
                if ( mBearingMap.key_exists( aID ) )
                {
                    return mBearingMap( aID );
                }
                else
                {
                    return mEmptyBearing ;
                }
            }

//-----------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_BEARINGDATA_HPP
