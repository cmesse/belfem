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
#include "cl_FEM_DofMgr_EigenValues.hpp"

namespace belfem
{
    namespace fem
    {
        class Kernel;
        class Postprocessor ;

        /**
         * this class creates the DOFs based on the passed equation object.
         */
        /**
         * @brief Main workhorse: DOF management, assembly coordination and the solver interface.
         *
         * @ingroup grp_fem_kernel
         * @see @ref fem_kernel_dof_manager_usage_guide
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

            //! data object for eigenvalues
            dofmgr::EigenValues  * mEigenValues ;

            //! flag telling if system has been initialized
            bool mInitializedFlag = false ;

            //! postprocessors ( L2 projections onto secondary fields ).
            //! Pushed in by the factory, owned and deleted by this dof manager
            Cell< Postprocessor * > mPostprocessors ;

            //! flag telling if the matrix has been computed.
            //! this flag will be reset once solve() is called
            bool mJacobianIsUpToDate = false ;

            //! this map is needed for the post processors
            Map< id_t, Element * > mElementMap ;

            Cell< string > mPostprocessorSourceFields ;
            Cell< string > mPostprocessorTargetFields ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *
             * @param aParent
             * @param aIndex :        corresponding index in kernel
             */
            DofManager(
                          Kernel         * aParent,
                    const index_t          aIndex );

//-----------------------------------------------------------------------------

            ~DofManager() override;

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
             sideset_integration_order() const override ;

//------------------------------------------------------------------------------

             uint
             block_integration_order() const override ;

//------------------------------------------------------------------------------

             IntegrationScheme
             integration_scheme() const override ;

//------------------------------------------------------------------------------

             bool
             is_initialized() const override
             {
                 return mInitializedFlag ;
             }

//------------------------------------------------------------------------------

            /**
             * return a specific dof
             */
            Dof *
            dof( const id_t aID ) override;

//------------------------------------------------------------------------------

            bool
            dof_exists( const id_t aID ) const override ;

//------------------------------------------------------------------------------

            /**
             * return the non-hanging dofs ( free and fixed; hanging dofs
             * live in DofData::hanging_dofs() )
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
            sideset( const id_t aID ) override;

//------------------------------------------------------------------------------

            /**
             * return a specific block
             */
            Block *
            block( const id_t aID ) override;

//------------------------------------------------------------------------------

            void
            extract_abstract_dofs_from_mesh() ;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Node * aNode , const uint aDofType )  const override;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Edge * aEdge , const uint aDofType )  const override;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Face * aFace , const uint aDofType )  const;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Facet * aFacet , const uint aDofType )  const override;

//-----------------------------------------------------------------------------

            void
            initialize() override;

//-----------------------------------------------------------------------------

            //! seed-mode overload: aSeedFreeDofsOnly = true skips the
            //! fixed-dof dof -> field write during init_dof_values(), so
            //! fields restored before this call ( load_memdump ) survive at
            //! Dirichlet nodes. The zero-arg override above must stay: it
            //! keeps base-pointer callers on full-mode virtual dispatch
            void
            initialize( const bool aSeedFreeDofsOnly );

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

//-----------------------------------------------------------------------------

            void
            compute_full_matrices();

//------------------------------------------------------------------------------

             /**
             * set the solver type for this field
             */
            void
            set_solver( const SolverParameters & aParams );

//------------------------------------------------------------------------------

            dofmgr::SolverData *
            solver_data() ;

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

            //! phase 1 of the two-phase iterative solve: residual
            //! of the committed state under the current assembly, no update
            void
            compute_residual();

//-----------------------------------------------------------------------------

            //! phase 2: linear solve + update consuming compute_residual()'s
            //! result; must be entered on every rank
            void
            solve_from_residual();

//-----------------------------------------------------------------------------

            real
            residual( const uint aIteration );

//-----------------------------------------------------------------------------

            /**
             * absolute residual norm of the last residual() call
             */
            real
            absolute_residual() const ;

//-----------------------------------------------------------------------------

            /**
             * pre-update residual of the entry state under the current
             * assembly ( Newton iterates only, QUIET_NAN otherwise ) — the
             * honest backtracking reference for the controller's line search
             */
            real
            pre_update_residual() const ;

//-----------------------------------------------------------------------------

            /**
             * fixed-point residual ||G(x)-x|| / ||x|| of the last Anderson
             * update ( REAL_MAX while Anderson is off )
             */
            real
            fixed_point_residual() const ;

//-----------------------------------------------------------------------------

            /**
             * true if the last solve failed softly ( soft-fail contract,
             * see solver::Wrapper ); rank-uniform
             */
            bool
            solve_failed() ;

//-----------------------------------------------------------------------------

            /**
             * collects node data from the others and send it to master
             */
            void
            collect_fields( const Cell< string > & aFieldLabels ) override ;

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
            distribute_fields( const Cell< string > & aFieldLabels ) override ;

//-----------------------------------------------------------------------------

            /**
             * perform a collect first, then a distribute
             */
            void
            synchronize_fields( const Cell< string > & aFieldLabels ) override ;

//-----------------------------------------------------------------------------

            void
            compute_volume_loads( const Vector< id_t > & aBlockIDs );

//-----------------------------------------------------------------------------
            /**
             * The JEDI Block Matrix System - May the Force be with you!
             *
             * Block form: [ J  D ] [ x ]   [ f ]
             *             [ E  I ] [ y ] = [ g ]
             *
             * where x = free DOFs, y = fixed DOFs (imposed values)
             *
             * Reduced system: J*x = f - D*y  (solves for unknowns)
             * Reaction forces: E*x + I*y = g (computes constraint forces!)
             *
             * In structural mechanics, the JEDI system literally computes
             * "the Force" at supports and constraints.
             */

            SpMatrix *
            system_matrix() ; // A: Free → Free coupling ( the Picard operator )

//-----------------------------------------------------------------------------

            SpMatrix *
            dirichlet() ;     // D matrix: Free → Fixed (Dirichlet BCs)

//-----------------------------------------------------------------------------

            SpMatrix *
            enforcement() ;   // E matrix: Fixed → Free (enforces constraints, computes forces)

//-----------------------------------------------------------------------------

            SpMatrix *
            imposition() ;    // I matrix: Fixed → Fixed (self-coupling of imposed values)

//-----------------------------------------------------------------------------

            SpMatrix *
            full_mass();

//-----------------------------------------------------------------------------

            SpMatrix *
            full_stiffness();

//-----------------------------------------------------------------------------

            void
            full_lhs( Vector< real > & aLHS );

//-----------------------------------------------------------------------------

            Vector< real >  &
            rhs_vector();

//------------------------------------------------------------------------------

            real
            rhs_norm();

//------------------------------------------------------------------------------

            void
            load_system( const string & aPath );

//------------------------------------------------------------------------------

            void
            save_system( const string & aPath );

//------------------------------------------------------------------------------

            void
            init_dof_values();

//------------------------------------------------------------------------------

            //! copy field values into the FREE dofs only ( fixed dofs keep
            //! their imposed values and do not write back into the field ).
            //! Call after the fields have been populated or restored — the
            //! initialize() pass seeds the dofs from whatever the fields
            //! held at that moment, and does not run again on this path
            //! ( only reset() via set_equation() clears the flag )
            void
            seed_dof_values();

//------------------------------------------------------------------------------

            void
            reset_convection();

//------------------------------------------------------------------------------
#ifdef BELFEM_HDF5

            void
            load_system( HDF5 & aFile );

            void
            save_system( HDF5 & aFile );

#endif
//-----------------------------------------------------------------------------

            bool
            sideset_exists( const id_t & aID ) const override ;

//-----------------------------------------------------------------------------

            bool
            block_exists( const id_t & aID ) const override ;

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
             Cell< Postprocessor * > &
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

//------------------------------------------------------------------------------
            /**
             * return the number of wetted nodes ( nodes on sidesets with a Neumann or Alpha BC )
             * @return
             */
             index_t
             number_of_wetted_nodes() const ;

//------------------------------------------------------------------------------

             void
             disconnect_dofs_from_mesh();

//------------------------------------------------------------------------------

            /**
             * expose the eigenvalue tool
             */
            dofmgr::EigenValues *
            eigen_values() ;

//------------------------------------------------------------------------------

            /**
             * returns direct access to one element. Needed for postprocessing
             */
             Element *
             element( const id_t aID );

             bool
             element_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------

            Cell< Dof * > &
            abstract_dofs();

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            reset();

//------------------------------------------------------------------------------

            void
            init_dofs( const bool aSeedFreeDofsOnly );

//------------------------------------------------------------------------------

            void
            init_work();

//------------------------------------------------------------------------------

            void
            init_matrices();

//-----------------------------------------------------------------------------

            void
            compute_rhs_vector();

//-----------------------------------------------------------------------------

            void
            compute_rhs_matrix();

//---------------------------------------------v--------------------------------

            void
            auto_set_materials();

//------------------------------------------------------------------------------

            void
            consolidate_dofs();

//------------------------------------------------------------------------------

            void
            consolidate_dofs( Element * aElement );

//-----------------------------------------------------------------------------

            void
            link_bearings_with_dofs();

//-----------------------------------------------------------------------------

            void
            create_element_map();

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

        inline dofmgr::SolverData *
        DofManager::solver_data()
        {
            return mSolverData ;
        }

//------------------------------------------------------------------------------

        inline Solver *
        DofManager::solver()
        {
            return mSolverData->solver() ;
        }

//------------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::system_matrix()
        {
            return mSolverData->system_matrix() ;
        }

//------------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::dirichlet()
        {
            return mSolverData->dirichlet() ;
        }

 //------------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::enforcement()
        {
            return mSolverData->enforcement() ;
        }

//------------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::imposition()
        {
            return mSolverData->imposition() ;
        }

//-----------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::full_mass()
        {
            return mSolverData->matrix( dofmgr::FullMass );
        }

//-----------------------------------------------------------------------------

        inline SpMatrix *
        DofManager::full_stiffness()
        {
            return mSolverData->matrix( dofmgr::FullStiffness );
        }

//------------------------------------------------------------------------------

        inline Vector< real >  &
        DofManager::rhs_vector()
        {
            return mSolverData->rhs_vector() ;
        }

//-----------------------------------------------------------------------------

        inline real
        DofManager::rhs_norm()
        {
            return mSolverData->rhs_norm() ;
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
        inline Cell< Postprocessor * > &
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

//------------------------------------------------------------------------------

        inline index_t
        DofManager::number_of_wetted_nodes() const
        {
            return mSideSetData->number_of_wetted_nodes() ;
        }

//------------------------------------------------------------------------------

        inline dofmgr::EigenValues *
        DofManager::eigen_values()
        {
            return mEigenValues ;
        }

        inline Element *
        DofManager::element( const id_t aID )
        {
            return mElementMap( aID );
        }

        inline bool
        DofManager::element_exists( const id_t aID ) const
        {
            return mElementMap.key_exists( aID );
        }
        
//------------------------------------------------------------------------------

        inline Cell< Dof * > &
        DofManager::abstract_dofs()
        {
            return mDofData->abstract_dofs() ;
        }

//-----------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_FEM_DOFMANAGER_HPP
