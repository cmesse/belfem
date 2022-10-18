//
// Created by christian on 7/15/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_SOLVERDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_SOLVERDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"

#include "cl_Solver.hpp"

#include "cl_IWG.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_HDF5.hpp"

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
            class BlockData ;
            class SideSetData ;

            class SolverData
            {
                //! the parent object
                DofManager * mParent;

                //! the kernel
                Kernel * mKernel;

                DofData     * mDofData ;
                BlockData   * mBlockData ;
                SideSetData * mSideSetData ;

                // my rank
                const proc_t mMyRank;

                Cell< Dof * > & mDOFs ;

                const index_t & mNumberOfFreeDofs ;
                const index_t & mNumberOfFixedDofs ;

                // master only
                Vector< index_t > mNumberOfFreeDofsPerProc;
                Vector< index_t > mNumberOfFixedDofsPerProc;

                // values for this proc, identical to mMyNumberOfFreeDofs for Master
                index_t mMyNumberOfFreeDofs ;
                index_t mMyNumberOfFixedDofs ;

                SpMatrix * mJacobian = nullptr ;
                SpMatrix * mDirichletMatrix = nullptr ;

                // indices on master matrix
                Cell< Vector< index_t > > mJacobianTable;
                Cell< Vector< index_t > > mDirichletTable;

                // left hand side ( field values or deltas )
                Vector< real > mLhsVector;
                Matrix< real > mLhsMatrix;

                // right hand side
                Vector< real > mRhsVector;
                Matrix< real > mRhsMatrix;

                // contains the norm of the rhs vector
                real mRhsNorm = BELFEM_QUIET_NAN ;

                // special field needed for timestepping
                Vector< real > mFieldValues ;

                // data container for convective data
                // contains e.g. the heatflux or external pressure
                Vector< real > mConvection ;

                // vector for volume loads, eg. imposed current
                Vector< real > mVolumeLoads ;

                //! the solver interface
                Solver * mSolver = nullptr ;

                //! contains values for initialization

                bool mUseResetValues = false ;
                Vector< real > mRhsVector0 ;
                Matrix< real > mRhsMatrix0 ;
                Vector< real > mJacobianValues0 ;
                Vector< real > mDirichletValues0 ;

                // needed for worst dof function
                Cell< Dof * > mFreeDofs ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                SolverData( DofManager * aParent,
                            DofData * aDofData,
                            BlockData * aBlockData,
                            SideSetData * aSideSetData );

//------------------------------------------------------------------------------

                ~SolverData() ;

//------------------------------------------------------------------------------

                void
                allocate_matrices();

//------------------------------------------------------------------------------

                void
                create_assembly_tables();

//------------------------------------------------------------------------------

                void
                reset();

//------------------------------------------------------------------------------

                void
                reset_matrices();

//------------------------------------------------------------------------------

                void
                reset_rhs_vector();

//------------------------------------------------------------------------------

                void
                reset_rhs_matrix();

//------------------------------------------------------------------------------

                void
                assemble_jacobian( Element * aElement,
                                   const Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

                void
                assemble_jacobian_and_rhs( Element * aElement,
                                   const Matrix< real > & aJacobian,
                                   const Vector< real > & aResidual );

//------------------------------------------------------------------------------

                void
                asseble_rhs( Element * aElement,
                             const Vector< real > & aRHS );
//------------------------------------------------------------------------------

                void
                assemble_volume_loads( Element * aElement,
                             const Vector< real > & aRHS );

//------------------------------------------------------------------------------

                void
                asseble_surface_loads( Element * aElement,
                                      const Vector< real > & aRHS );

//------------------------------------------------------------------------------

                void
                asseble_rhs( Element * aElement,
                             const Matrix< real > & aRHS );

//------------------------------------------------------------------------------

                void
                collect_jacobian();

//------------------------------------------------------------------------------

                void
                collect_rhs_vector();

//------------------------------------------------------------------------------

                void
                collect_vector( Vector< real > & aVector );

//------------------------------------------------------------------------------

                void
                collect_rhs_matrix();

//------------------------------------------------------------------------------

                void
                update_field_values();

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

//------------------------------------------------------------------------------

                /**
                 * solve the system
                 */
                void
                solve();

//------------------------------------------------------------------------------

                /**
                 * return the residual
                 */
                real
                residual( const uint aIteration ) ;

//------------------------------------------------------------------------------

                SpMatrix *
                jacobian();

//------------------------------------------------------------------------------

                SpMatrix *
                dirichlet();

//------------------------------------------------------------------------------

                Vector< real > &
                rhs_vector() ;

//------------------------------------------------------------------------------

                Vector< real > &
                volume_loads() ;

//------------------------------------------------------------------------------

                Vector< real > &
                surface_loads() ;

//------------------------------------------------------------------------------

                index_t
                my_number_of_free_dofs() const ;

//------------------------------------------------------------------------------

                index_t
                my_number_of_fixed_dofs() const ;

//------------------------------------------------------------------------------

                void
                save_system( const string & aPath );

//------------------------------------------------------------------------------

                void
                load_system( const string & aPath );

//------------------------------------------------------------------------------

                void
                remember_initialization_values( const bool aSaveRHS );

//------------------------------------------------------------------------------

                void
                use_reset_values( const bool aFlag );


//-----------------------------------------------------------------------------

                void
                print_worst_dof();

//------------------------------------------------------------------------------
#ifdef BELFEM_HDF5
                void
                save_system( HDF5 & aFile );

                void
                load_system( HDF5 & aFile );
#endif

//------------------------------------------------------------------------------
            private:
//------------------------------------------------------------------------------

                void
                compute_element_dof_connectivity( Vector< id_t > & aData );

//------------------------------------------------------------------------------

                void
                compute_dof_element_connectivity( Vector< id_t > & aData );

//----------------------------------------------------------------------------

                void
                compute_dof_dof_connectivity(
                        const Vector< id_t > & aDofWiseData,
                        const Vector< id_t > & aElementWiseData,
                              Vector< id_t > & aConnectivity );

//----------------------------------------------------------------------------

                void
                unite_dofs(
                        const Cell< Vector<id_t > > & aConnectivities,
                               Vector< id_t > & aConnectivity );

//----------------------------------------------------------------------------

                void
                populate_graph( const Vector< id_t > & aData,
                                const bool aFixedFlag,
                                Cell< graph::Vertex * > & aGraph );

//-----------------------------------------------------------------------------

                void
                collect_fields( Cell< mesh::Field * > & aFields );

//------------------------------------------------------------------------------
            };
//------------------------------------------------------------------------------

            inline Solver *
            SolverData::solver()
            {
                return mSolver ;
            }

//------------------------------------------------------------------------------

            inline SpMatrix *
            SolverData::jacobian()
            {
                return mJacobian ;
            }

//------------------------------------------------------------------------------

            inline SpMatrix *
            SolverData::dirichlet()
            {
                return mDirichletMatrix ;
            }

//------------------------------------------------------------------------------

            inline Vector< real > &
            SolverData::rhs_vector()
            {
                return mRhsVector ;
            }

//------------------------------------------------------------------------------

            inline Vector< real > &
            SolverData::volume_loads()
            {
                return mVolumeLoads ;
            }

//------------------------------------------------------------------------------

            inline Vector< real > &
            SolverData::surface_loads()
            {
                return mConvection ;
            }

//------------------------------------------------------------------------------

            inline index_t
            SolverData::my_number_of_free_dofs() const
            {
                return mMyNumberOfFreeDofs ;
            }

//------------------------------------------------------------------------------

            inline index_t
            SolverData::my_number_of_fixed_dofs() const
            {
                return mMyNumberOfFixedDofs ;
            }

//------------------------------------------------------------------------------

            inline void
            SolverData::use_reset_values( const bool aFlag )
            {
                mUseResetValues = aFlag ;
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_SOLVERDATA_HPP
