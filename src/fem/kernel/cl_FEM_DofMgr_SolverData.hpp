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

#ifndef BELFEM_CL_FEM_DOFMGR_SOLVERDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_SOLVERDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_ShiftRegister.hpp"
#include "cl_SpMatrix.hpp"

#include "cl_Solver.hpp"

#include "cl_IWG.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_FEM_DofMgr_DofData.hpp"
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
            /**
             * Matrix type enumeration in JEDI order
             *
             * The JEDI system computes "the Force" - literally!
             * In structural mechanics, E*x + I*y = reaction forces at constraints
             */
            enum MatrixType
            {
                System,       // A: Free → Free ( Picard operator; the matrix SolverData::system_matrix() returns )
                Jacobian,     // A + dA/dx term: Free → Free ( Newton tangent, shares the pattern of A;
                              //   assemble_newton() adds the derivative term )
                Enforcement,  // E: Fixed → Free (enforces constraints, computes forces)
                Dirichlet,    // D: Free → Fixed (Dirichlet boundary conditions)
                Imposition,   // I: Fixed → Fixed (self-coupling of imposed values)
                FullMass,         // full mass matrix, including fixed dofs
                FullStiffness     // full stiffness matrix, including fixed dofs
            };

            string
            to_string( const MatrixType aType );

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
                const proc_t mCommRank;
                const proc_t mCommSize ;
                
                Cell< Dof * > & mDOFs ;

                const index_t & mNumberOfFreeDofs ;
                const index_t & mNumberOfFixedDofs ;
                //const index_t & mNumberOfHangingDofs ;

                const index_t & mMyNumberOfFreeDofs ;
                const index_t & mMyNumberOfFixedDofs ;
                const index_t & mMyNumberOfHangingDofs ;

                // The JEDI Block Matrix System - Computing the Force Vector
                // Block system: [ J  D ] [ x ]   [ f ]
                //               [ E  I ] [ y ] = [ g ]
                //
                // Physical interpretation in structural mechanics:
                //   First equation:  J*x = f - D*y  (solve for free DOFs)
                //   Second equation: E*x + I*y = g  (computes reaction forces!)
                //
                // The JEDI system literally computes "the Force" at constraints!
                SpMatrix * mSystemMatrix   = nullptr ;  // A: Free → Free coupling
                SpMatrix * mEnforcementMatrix = nullptr ; // E: Fixed → Free (enforces constraints, computes forces)
                SpMatrix * mDirichletMatrix  = nullptr ;  // D: Free → Fixed (Dirichlet boundary conditions)
                SpMatrix * mImpositionMatrix = nullptr ;  // I: Fixed → Fixed (self-coupling of imposed values)

                SpMatrix * mFullMassMatrix = nullptr ;
                SpMatrix * mFullStiffnessMatrix = nullptr ;
                SpMatrix * mJacobianMatrix = nullptr ;

                // indices on master matrix
                Cell< Vector< index_t > > mSystemTable;
                Cell< Vector< index_t > > mJacobianTable;
                Cell< Vector< index_t > > mEnforcementTable;
                Cell< Vector< index_t > > mDirichletTable;
                Cell< Vector< index_t > > mImpositionTable;

                Cell< Vector< index_t > > mFullMatrixTable;

                // left hand side ( field values or deltas )
                Vector< real > mLhsVector;
                Matrix< real > mLhsMatrix;

                // right hand side
                Vector< real > mRhsVector;
                Matrix< real > mRhsMatrix;

                //! Picard-branch scratch: preserves the pre-update residual r
                //! across the increment solve so residual() can read it back.
                //! Newton does NOT use it — b lives in mRhsOriginal ( below )
                Vector< real > mRhsBackup;

                //! durable copy of the load-adjusted b, taken by
                //! compute_residual() on the Newton branch: the post-update
                //! recompute in solve_from_residual() restores b from HERE,
                //! never from mRhsBackup — the two backups carry incompatible
                //! meanings across the algorithms
                Vector< real > mRhsOriginal;

                //! retained scratch for the fixed-dof value scrape in
                //! compute_residual() — a per-call Vector on the iterate
                //! path violates the no-hidden-allocation rule
                Vector< real > mFixedValuesScratch;

                //! retained field-pointer collection for the two-phase
                //! solve; refilled by compute_residual(), consumed by
                //! solve_from_residual()
                Cell< mesh::Field * > mCollectedFields;

                //! two-phase handshake: set ( on every rank ) by
                //! compute_residual(), consumed by solve_from_residual(),
                //! invalidated by a fresh assembly. Guards the one-shot
                //! contract — a second compute_residual() on one assembly
                //! would double-add the convection / volume / Dirichlet
                //! load terms into the RHS
                bool mResidualReady = false ;

                // contains the norm of the rhs vector
                real mRhsNorm = BELFEM_QUIET_NAN ;

                //! absolute residual norm ||A x - b|| of the last residual()
                //! call ( the relative measure divides by mRhsNorm; the
                //! absolute-tolerance escape needs the raw value )
                real mAbsoluteResidual = BELFEM_REAL_MAX ;

                //! fixed-point residual ||G(x)-x|| / ||x|| of the last
                //! anderson_update() call ( REAL_MAX while Anderson is off;
                //! rank-uniform after the broadcast in residual() )
                real mFixedPointResidual = BELFEM_REAL_MAX ;

                //! PRE-update residual ||A x - b|| / ||b|| of the ENTRY state
                //! under the CURRENT assembly, captured in the Newton branch
                //! before the solve. This is the omega -> 0 limit of any
                //! line-search trial and therefore the only honest
                //! backtracking reference: the controller's previous-iterate
                //! epsilon was measured under the PREVIOUS assembly, and when
                //! the two assemblies disagree by more than the acceptance
                //! band, that stale reference rejects every trial flat in
                //! omega ( greg5 out2 trace ). QUIET_NAN on Picard iterates;
                //! rank-uniform after the broadcast in residual()
                real mPreUpdateResidual = BELFEM_QUIET_NAN ;

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
                Vector< real > mSystemValues0 ;
                Vector< real > mJacobianValues0 ;
                Vector< real > mDirichletValues0 ;
                Vector< real > mEnforcementValues0 ;
                Vector< real > mImpositionValues0 ;

                //! to be opened to user
                bool mUseJediForce  = false ;
                bool mUseFullForce  = false ;
                Cell< index_t > mWorkDofIndices ;

                //! Anderson mixing for the Picard branch ( opt-in, master
                //! side only; see todo/closed/anderson_picard_acceleration_plan.md ).
                //! 0 = off: the legacy relaxed update runs untouched
                uint mAndersonDepth = 0 ;

                //! committed history, ( 0 ) = newest; allocated only while
                //! depth > 0 ( ShiftRegister rejects zero capacity )
                ShiftRegister< Vector< real > > * mAndersonX = nullptr ;
                ShiftRegister< Vector< real > > * mAndersonR = nullptr ;

                //! staged candidate pair: written by anderson_update(), moved
                //! into the history by anderson_commit() once the controller
                //! accepts the trial, dropped by anderson_discard()
                Vector< real > mAndersonXStage ;
                Vector< real > mAndersonRStage ;
                bool mAndersonStaged = false ;

                //! scratch for the mixing step, sized on first use
                Matrix< real > mAndersonDeltaR ;
                Vector< real > mAndersonRhs ;
                Vector< real > mAndersonWork ;
                Vector< real > mAndersonGamma ;
                Vector< real > mAndersonColNorm ;
                Vector< real > mAndersonXNew ;

                //! Anderson-mixed Picard update, replaces the legacy relaxed
                //! update loop while depth > 0 ( master only )
                void
                anderson_update( IWG * aIWG, Cell< mesh::Field * > & aFields );

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
                extract_graph_from_mesh( Vector< id_t > & aGraphData );

//------------------------------------------------------------------------------

                void
                allocate_matrices(
                    const Vector< id_t > & aGraphData,
                    Graph & aFreeDofs,
                    Graph & aFixedDofs );

//------------------------------------------------------------------------------

                void
                compute_memory();

//------------------------------------------------------------------------------

                void
                create_assembly_tables();

//------------------------------------------------------------------------------

                void
                reset();

//------------------------------------------------------------------------------

                void
                reset_convection();

//------------------------------------------------------------------------------

                void
                reset_matrices( const bool aFullForce = false );

//------------------------------------------------------------------------------

                void
                reset_rhs_vector();

//------------------------------------------------------------------------------

                void
                reset_rhs_matrix();

//------------------------------------------------------------------------------

                void
                assemble_jacobian( Element * aElement,
                                   const Matrix< real > & aJacobian);

//------------------------------------------------------------------------------

                void
                assemble_full_matrices( Element * aElement,
                                   const Matrix< real > & aMass,
                                   const Matrix< real > & aStiffness );

//------------------------------------------------------------------------------

                void
                assemble_newton( Element * aElement,
                                   const Matrix< real > & adJdx );

//------------------------------------------------------------------------------

                void
                assemble_rhs( Element * aElement,
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
                assemble_rhs( Element * aElement,
                             const Matrix< real > & aRHS );

//------------------------------------------------------------------------------

                void
                collect_matrices( const bool aFullForce = false );

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
                set_solver( const SolverParameters & aParams );

//------------------------------------------------------------------------------

                /**
                 * expose the solver
                 */
                Solver *
                solver();

//------------------------------------------------------------------------------

                /**
                 * solve the system. For the iterative vector-RHS path this is
                 * the fused composition compute_residual() +
                 * solve_from_residual(); Direct mode and matrix-RHS callers
                 * keep the single-call contract unchanged.
                 */
                void
                solve();

//------------------------------------------------------------------------------

                /**
                 * Phase 1 of the two-phase iterative solve:
                 * prepare the load-adjusted RHS b, capture ||b||, and
                 * overwrite mRhsVector with the residual r = A x - b of the
                 * COMMITTED state under the current assembly — without
                 * solving. residual() then reports the head certificate.
                 * One-shot per assembly ( load terms would double-add );
                 * requires a fresh compute_jacobian_and_rhs.
                 */
                void
                compute_residual();

//------------------------------------------------------------------------------

                /**
                 * Phase 2: consume the residual left by compute_residual(),
                 * run the linear solve and the algorithm's update ( relaxed /
                 * Anderson Picard, or the Newton step with its post-update
                 * recompute ). Must be entered on EVERY rank — the collective
                 * barrier pairing lives here.
                 */
                void
                solve_from_residual();

//------------------------------------------------------------------------------

                /**
                 * invalidate the two-phase handshake; called by a fresh
                 * assembly
                 */
                void
                invalidate_residual()
                {
                    mResidualReady = false ;
                }

//------------------------------------------------------------------------------

                /**
                 * return the residual
                 */
                real
                residual( const uint aIteration ) ;

//------------------------------------------------------------------------------

                /**
                 * absolute residual norm of the last residual() call
                 * ( rank-uniform, broadcast alongside the relative value )
                 */
                real
                absolute_residual() const
                {
                    return mAbsoluteResidual ;
                }

//------------------------------------------------------------------------------

                /**
                 * pre-update residual of the entry state under the current
                 * assembly ( Newton iterates only, QUIET_NAN otherwise;
                 * rank-uniform after residual() ). The honest backtracking
                 * reference — see the member comment
                 */
                real
                pre_update_residual() const
                {
                    return mPreUpdateResidual ;
                }

//------------------------------------------------------------------------------

                /**
                 * fixed-point residual ||G(x)-x|| / ||x|| of the last
                 * Anderson update ( rank-uniform after residual() )
                 */
                real
                fixed_point_residual() const
                {
                    return mFixedPointResidual ;
                }

//------------------------------------------------------------------------------

                /**
                 * activate ( depth > 0 ) or deactivate ( depth = 0 ) Anderson
                 * mixing of the Picard branch; allocates the history registers
                 */
                void
                set_anderson_depth( const uint aDepth );

                /**
                 * move the staged ( x, r ) pair into the history; called by
                 * the controller once the trial iterate is accepted. No-op if
                 * mixing is off or nothing is staged.
                 */
                void
                anderson_commit();

                /**
                 * drop the staged pair ( trial was rejected )
                 */
                void
                anderson_discard();

                /**
                 * invalidate the history ( new attempt, algorithm switch,
                 * restore, ... ); also drops any staged pair
                 */
                void
                anderson_clear();

//------------------------------------------------------------------------------

                SpMatrix *
                system_matrix();

//------------------------------------------------------------------------------

                SpMatrix *
                enforcement();

//------------------------------------------------------------------------------

                SpMatrix *
                dirichlet();

//------------------------------------------------------------------------------

                SpMatrix *
                imposition();

//------------------------------------------------------------------------------

                SpMatrix *
                matrix( const MatrixType aType );

//------------------------------------------------------------------------------

                Cell< Vector< index_t > > &
                tables( const MatrixType aType );

//------------------------------------------------------------------------------

                Vector< real > &
                rhs_vector() ;

//------------------------------------------------------------------------------

                real
                rhs_norm() ;

//------------------------------------------------------------------------------

                Vector< real > &
                volume_loads() ;

//------------------------------------------------------------------------------

                Vector< real > &
                surface_loads() ;

//------------------------------------------------------------------------------

                index_t
                number_of_free_dofs() const ;

//------------------------------------------------------------------------------

                index_t
                number_of_fixed_dofs() const ;

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

//------------------------------------------------------------------------------

                void
                use_jedi_force( const bool aFlag ) ;

//------------------------------------------------------------------------------

                void
                use_full_force( const bool aFlag );

#ifdef BELFEM_HDF5
                void
                save_system( HDF5 & aFile );

                void
                load_system( HDF5 & aFile );
#endif

//----------------------------------------------------------------------------

                void
                populate_graph( const Vector< id_t > &    aData,
                                const MatrixType          aMatrixType,
                                Graph                & aFreeDofs,
                                Graph                & aFixedDofs,
                                bool aLinkToSelf = true );

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


//-----------------------------------------------------------------------------

                void
                collect_fields( Cell< mesh::Field * > & aFields );

//------------------------------------------------------------------------------

                void
                compute_hanging_dofs( Cell< mesh::Field * > & aFields );

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
            SolverData::system_matrix()
            {
                return mSystemMatrix ;
            }

//------------------------------------------------------------------------------

            inline SpMatrix *
            SolverData::enforcement()
            {
                return mEnforcementMatrix ;
            }

//------------------------------------------------------------------------------

            inline SpMatrix *
            SolverData::dirichlet()
            {
                return mDirichletMatrix ;
            }

//------------------------------------------------------------------------------

            inline SpMatrix *
            SolverData::imposition()
            {
                return mImpositionMatrix ;
            }

//------------------------------------------------------------------------------

            inline Vector< real > &
            SolverData::rhs_vector()
            {
                return mRhsVector ;
            }

//------------------------------------------------------------------------------

            inline real
            SolverData::rhs_norm()
            {
                return mRhsNorm ;
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
            SolverData::number_of_free_dofs() const
            {
                return mNumberOfFreeDofs ;
            }

//------------------------------------------------------------------------------

            inline index_t
            SolverData::number_of_fixed_dofs() const
            {
                return mNumberOfFixedDofs ;
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

            inline SpMatrix *
            SolverData::matrix( const MatrixType aType )
            {
                switch( aType )
                {
                    case System:      return mSystemMatrix ;
                    case Jacobian:      return mJacobianMatrix ;
                    case Enforcement:   return mEnforcementMatrix ;
                    case Dirichlet:     return mDirichletMatrix ;
                    case Imposition:    return mImpositionMatrix ;
                    case FullMass:      return mFullMassMatrix ;
                    case FullStiffness: return mFullStiffnessMatrix ;
                    default: return nullptr ;
                }
            }

//------------------------------------------------------------------------------

            inline Cell< Vector< index_t > > &
            SolverData::tables( const MatrixType aType )
            {
                switch ( aType )
                {
                    case System: return mSystemTable ;
                    case Jacobian: return mJacobianTable ;
                    case Enforcement: return mEnforcementTable ;
                    case Dirichlet: return mDirichletTable ;
                    case Imposition: return mImpositionTable ;
                    case FullMass:      return mFullMatrixTable ;
                    case FullStiffness: return mFullMatrixTable ;
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid matrix type");
                        return mSystemTable ;
                    }
                }
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_SOLVERDATA_HPP
