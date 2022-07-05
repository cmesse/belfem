//
// Created by christian on 7/14/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_SIDESETDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_SIDESETDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "commtools.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_SideSet.hpp"

namespace belfem
{
    class Mesh;

    namespace fem
    {

        class Kernel;

        class DofManager;

        namespace dofmgr
        {
            class DofData ;

            class SideSetData
            {
                //! the parent object
                DofManager * mParent;

                //! the kernel
                Kernel * mKernel;

                //! the mesh this problem runs on
                Mesh * mMesh;

                // my rank
                const proc_t mMyRank;

                Cell< SideSet * > mSideSets;
                Map< id_t, SideSet * > mSideSetMap;

                // an empty sidset that is called if this set does not exist
                // on current proc
                SideSet * mEmptySideset = nullptr;

                // counts how many nodes on all sidesets experience convection
                index_t mNumberOfConvectionNodes = 0 ;

//-----------------------------------------------------------------------------
            public:
//-----------------------------------------------------------------------------

                SideSetData( DofManager * aParent ) ;

//-----------------------------------------------------------------------------

                ~SideSetData();

//-----------------------------------------------------------------------------

                void
                create_sidesets();

//-----------------------------------------------------------------------------

                /**
                 * return a specific sideset
                 */
                SideSet *
                sideset( const id_t aID ) ;

//-----------------------------------------------------------------------------

                /**
                 * return the sideset container
                 */
                Cell< SideSet * > &
                sidesets() ;

//-----------------------------------------------------------------------------

                /**
                 * check if a specific sideset exists on this proc
                 */
                 bool
                 sideset_exists( const id_t aID ) const;

//-----------------------------------------------------------------------------

                /**
                 *  check if user has set the wetted sidesets
                 */
                 void
                 detect_wetted_sidesets();

//-----------------------------------------------------------------------------

                 void
                 count_wetted_nodes();

//-----------------------------------------------------------------------------

                /**
                 * this function creates a special field for the alpha-bc
                 * if it exists. Called by DofManager::init_dofs();
                 */
                void
                create_alpha_fields();


//------------------------------------------------------------------------------

                /**
                 * called by DofManager::init_dofs();
                 */
                void
                set_boundary_conditions() ;

//-----------------------------------------------------------------------------

                //! restore factory settings
                void
                reset();

//----------------------------------------------------------------------------
            };

//----------------------------------------------------------------------------

            inline SideSet *
            SideSetData::sideset( const id_t aID )
            {
                if ( mSideSetMap.key_exists( aID ) )
                {
                    return mSideSetMap( aID );
                }
                else
                {
                    mEmptySideset->set_mesh_id( aID );
                    return mEmptySideset ;
                }
            }

//----------------------------------------------------------------------------

            inline Cell< SideSet * > &
            SideSetData::sidesets()
            {
                return mSideSets ;
            }

//------------------------------------------------------------------------------

            inline bool
            SideSetData::sideset_exists( const id_t aID ) const
            {
                return mSideSetMap.key_exists( aID );
            }

//------------------------------------------------------------------------------

        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_SIDESETDATA_HPP
