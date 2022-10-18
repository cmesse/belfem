//
// Created by christian on 7/7/21.
//

#ifndef BELFEM_CL_FEM_DOFMANAGER_HPP
#define BELFEM_CL_FEM_DOFMANAGER_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Map.hpp"
#include "cl_Mesh.hpp"
#include "cl_IWG.hpp"

#include "cl_FEM_Dof.hpp"
#include "cl_FEM_Bearing.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_SideSet.hpp"

#include "en_IntegrationScheme.hpp"

#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_DofMgr_Parameters.hpp"
#include "cl_FEM_DofMgr_DofData.hpp"
#include "cl_FEM_DofMgr_BlockData.hpp"
#include "cl_FEM_DofMgr_SideSetData.hpp"
#include "cl_FEM_DofMgr_BearingData.hpp"
#include "cl_FEM_DofMgr_FieldData.hpp"
#include "cl_FEM_DofMgr_SolverData.hpp"

#include "cl_Bitset.hpp"

namespace belfem
{
    namespace fem
    {
        class Kernel;

        /**
         * this class creates the DOFs based on the passed equation object.
         */
        class DofManager : public DofManagerBase
        {
            //! index of this manager on kernel
            const index_t mIndex ;

            //! the parameter object
            dofmgr::Parameters  * mParams ;

            //! data object for DOF handling
            dofmgr::DofData     * mDofData ;

            //! data object for block handling
            dofmgr::BlockData   * mBlockData ;

            //! data object for sideset handling
            dofmgr::SideSetData * mSideSetData ;

            //! data object for bearings
            dofmgr::BearingData * mBearingData ;

            //! data object for fields
            dofmgr::FieldData   * mFieldData ;

            //! data object for solver
            dofmgr::SolverData  * mSolverData ;

            //! flag telling if system has been initialized
            bool mInitializedFlag = false ;

            //! a DOF manager can contain other dof managers
            //! these are used for L2 projection
            //! postprocessors are owned and destroyed by the kernel
            Cell< DofManager * > mPostprocessors ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *
             * @param aParent
             * @param aMesh
             * @param aIndex :        corresponding index in kernel
             */
            DofManager(
                          Kernel         * aParent,
                    const index_t          aIndex );

//-----------------------------------------------------------------------------

            ~DofManager();

//-----------------------------------------------------------------------------

            /**
             * select the equation that is to be solved
             */
             void
             set_equation( IWG * aIWG );

//-----------------------------------------------------------------------------

            /**
             * a test routine for development
             */
             void
             print( const proc_t aRank=0 ) ;

//------------------------------------------------------------------------------

             uint
             sideset_integration_order() const ;

//------------------------------------------------------------------------------

             uint
             block_integration_order() const ;

//------------------------------------------------------------------------------

             IntegrationScheme
             integration_scheme() const ;

//------------------------------------------------------------------------------

            /**
             * return a specific dof
             */
            Dof *
            dof( const id_t aID );

//------------------------------------------------------------------------------

            bool
            dof_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------

            /**
             * return all dofs
             */
            Cell< Dof * > &
            dofs();

//-----------------------------------------------------------------------------

            /**
             * return a specific bearing
             */
            Bearing *
            bearing( const id_t aID );

//-----------------------------------------------------------------------------

            /**
             * return a specific sideset
             */
            SideSet *
            sideset( const id_t aID );

//------------------------------------------------------------------------------

            /**
             * return a specific block
             */
            Block *
            block( const id_t aID );

//------------------------------------------------------------------------------

            bool
            enforce_linear_interpolation() const ;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Node * aNode , const uint aDofType )  const;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Edge * aEdge , const uint aDofType )  const;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Face * aFace , const uint aDofType )  const;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Facet * aFacet , const uint aDofType )  const;

//-----------------------------------------------------------------------------

            void
            initialize();

//------------------------------------------------------------------------------

            void
            zero();

//-----------------------------------------------------------------------------

            void
            compute_jacobian( const bool aReset=true );

//-----------------------------------------------------------------------------

            void
            compute_rhs( const bool aReset=true );

//-----------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs( const bool aReset=true );

//------------------------------------------------------------------------------

             /**
             * set the solver type for this field
             */
            void
            set_solver( const SolverType aSolver );

//------------------------------------------------------------------------------

            /**
             * expose the solver
             */
            Solver *
            solver();

//-----------------------------------------------------------------------------

            void
            solve();

//-----------------------------------------------------------------------------

            real
            residual( const uint aIteration );

//-----------------------------------------------------------------------------

            /**
             * collects node data from the others and send it to master
             */
            void
            collect_fields( const Cell< string > & aFieldLabels ) ;

//-----------------------------------------------------------------------------

            /**
             * collects node data from the others and send it to master
             */
            void
            collect_field( const string & aFieldLabel ) ;

//-----------------------------------------------------------------------------

            /**
             * sends field data from master to the others
             */
            void
            distribute_fields( const Cell< string > & aFieldLabels ) ;

//-----------------------------------------------------------------------------

            /**
             * perform a collect first, then a distribute
             */
            void
            synchronize_fields( const Cell< string > & aFieldLabels ) ;

//-----------------------------------------------------------------------------

            void
            compute_volume_loads( const Vector< id_t > & aBlockIDs );

//-----------------------------------------------------------------------------

            SpMatrix *
            jacobian() ;

//-----------------------------------------------------------------------------

            Vector< real >  &
            rhs_vector();

//------------------------------------------------------------------------------

            void
            load_system( const string & aPath );

//------------------------------------------------------------------------------

            void
            save_system( const string & aPath );

//------------------------------------------------------------------------------

            void
            init_dof_values();

            void
            print_worst_dof() ;

            void
            write_residuals_to_mesh();

//------------------------------------------------------------------------------
#ifdef BELFEM_HDF5

            void
            load_system( HDF5 & aFile );

            void
            save_system( HDF5 & aFile );

#endif
//-----------------------------------------------------------------------------

            bool
            sideset_exists( const id_t & aID ) const ;

//-----------------------------------------------------------------------------

            bool
            block_exists( const id_t & aID ) const ;

//-----------------------------------------------------------------------------

            /**
             * check that the fields demanded by the IWG exist on the mesh
             * and create them otherwise
             */
            void
            create_fields( IWG * aIWG );

//-----------------------------------------------------------------------------

            /**
             * expose the projector arrays
             */
             Cell< DofManager * > &
             postprocessors();

//-----------------------------------------------------------------------------

            /**
             * compute the matrices for the projections
             */
             void
             initialize_postprocessors();

//-----------------------------------------------------------------------------

            /**
             * perform the L2 projections for the secondary fields
             */
            void
            postprocess();

//------------------------------------------------------------------------------

            void
            fix_node_dofs_on_ghost_sidesets() ;

//------------------------------------------------------------------------------

            /**
             * return the blocks on this dof manager
             * @return
             */
            Cell< Block * > &
            blocks();

//------------------------------------------------------------------------------

            /**
             * return the sidesets on this dof manager
             * @return
             */
            Cell< SideSet * > &
            sidesets();

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            reset();

//------------------------------------------------------------------------------

            void
            init_dofs();

//------------------------------------------------------------------------------

            void
            init_matrices();

//-----------------------------------------------------------------------------

            void
            compute_rhs_vector();

//-----------------------------------------------------------------------------

            void
            compute_rhs_matrix();

//-----------------------------------------------------------------------------

            void
            auto_set_materials();

//-----------------------------------------------------------------------------

        };

//-----------------------------------------------------------------------------

        inline uint
        DofManager::sideset_integration_order() const
        {
            return mParams->sideset_integration_order() ;
        }

//------------------------------------------------------------------------------

        inline IntegrationScheme
        DofManager::integration_scheme() const
        {
            return mParams->integration_scheme() ;
        }

//------------------------------------------------------------------------------

        inline uint
        DofManager::block_integration_order() const
        {
            return mParams->block_integration_order() ;
        }


//------------------------------------------------------------------------------

        inline Dof *
        DofManager::dof( const id_t aID )
        {
            return mDofData->dof( aID );
        }

//------------------------------------------------------------------------------

        inline Cell< Dof * > &
        DofManager::dofs()
        {
            return mDofData->dofs() ;
        }

//------------------------------------------------------------------------------

        inline bool
        DofManager::enforce_linear_interpolation() const
        {
            return mParams->enforce_linear();
        }

//------------------------------------------------------------------------------

        inline id_t
        DofManager::calculate_dof_id( const mesh::Node * aNode , const uint aDofType )  const
        {
            return mDofData->node_dof_id( aNode->id(), aDofType );
        }

//------------------------------------------------------------------------------

        inline id_t
        DofManager::calculate_dof_id( const mesh::Facet * aFacet , const uint aDofType )  const
        {
            return mDofData->lambda_dof_id( aFacet->id(), aDofType );
        }

//------------------------------------------------------------------------------

        inline id_t
        DofManager::calculate_dof_id( const mesh::Edge * aEdge , const uint aDofType ) const
        {
            return mDofData->edge_dof_id( aEdge->id(), aDofType );
        }

//------------------------------------------------------------------------------

        inline id_t
        DofManager::calculate_dof_id( const mesh::Face * aFace , const uint aDofType ) const
        {
            return mDofData->face_dof_id( aFace->id(), aDofType );
        }

//------------------------------------------------------------------------------

        inline Bearing *
        DofManager::bearing( const id_t aID )
        {
            return mBearingData->bearing( aID );
        }

//------------------------------------------------------------------------------

        inline SideSet *
        DofManager::sideset( const id_t aID )
        {
            return mSideSetData->sideset( aID );
        }

//------------------------------------------------------------------------------

        inline Block *
        DofManager::block( const id_t aID )
        {
            return mBlockData->block( aID );
        }

//------------------------------------------------------------------------------

        inline Solver *
        DofManager::solver()
        {
            return mSolverData->solver() ;
        }

//------------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::jacobian()
        {
            return mSolverData->jacobian() ;
        }

//------------------------------------------------------------------------------

        inline Vector< real >  &
        DofManager::rhs_vector()
        {
            return mSolverData->rhs_vector() ;
        }

//-----------------------------------------------------------------------------

        inline bool
        DofManager::sideset_exists( const id_t & aID ) const
        {
            return mSideSetData->sideset_exists( aID );
        }

//-----------------------------------------------------------------------------

        inline bool
        DofManager::block_exists( const id_t & aID ) const
        {
            return mBlockData->block_exists( aID );
        }

//------------------------------------------------------------------------------

        inline bool
        DofManager::dof_exists( const id_t aID ) const
        {
            return mDofData->dof_exists( aID );
        }

//-----------------------------------------------------------------------------

        /**
         * expose the projector arrays
         */
        inline Cell< DofManager * > &
        DofManager::postprocessors()
        {
            return mPostprocessors ;
        }

//------------------------------------------------------------------------------

        inline Cell< Block * > &
        DofManager::blocks()
        {
            return mBlockData->blocks() ;
        }

//------------------------------------------------------------------------------

        inline  Cell< SideSet * > &
        DofManager::sidesets()
        {
            return mSideSetData->sidesets() ;
        }

//-----------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_FEM_DOFMANAGER_HPP
