//
// Created by Christian Messe on 02.11.19.
//

#ifndef BELFEM_CL_FEM_FIELD_HPP
#define BELFEM_CL_FEM_FIELD_HPP

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

#include "cl_Solver.hpp"
#include "en_IntegrationScheme.hpp"
#include "cl_FEM_DofManagerBase.hpp"

namespace belfem
{
    class SpMatrix;

    namespace fem
    {
        class Kernel;

//------------------------------------------------------------------------------

        // todo: this class is deprecated and should not be used anymore
        class Field : public DofManagerBase
        {
            // Kernel * mParent;

            // Mesh * mMesh;

            // IWG * mIWG = nullptr;

            // const proc_t mMyRank;

            const Vector< id_t > mBlockIDs;
            const uint mNumberOfDOFsPerNode;
            const uint mNumberOfDOFsPerEdge;

            // special purpose flag to enforce a linear interpolation
            // regardless of the mesh
            const bool mEnforceLinear = false ;

            // flag on parent, where field is in mesh
            // const uint mFieldFlag;

            const uint mBlockIntegrationOrder = 0;
            const uint mSideSetIntegrationOrder = 0;
            const IntegrationScheme mIntegrationScheme = IntegrationScheme::UNDEFINED ;

            Solver * mSolver = nullptr ;


            // this map is linked to blocks on mesh
            Map< id_t, index_t > mBlockMap;
            Map< id_t, index_t > mSideSetMap;

            Map< id_t, index_t > mDofMap;
            Map< id_t, index_t > mBearingMap;

            Cell< fem::Block * >   mBlocks;

            // this map is linked to fem Blocks on field
            Map< id_t, index_t > mFemBlockMap;

            Cell< fem::SideSet * > mSideSets;
            Cell< fem::Bearing * > mBearings;

            fem::SideSet * mEmptySideset = nullptr;
            fem::Block   * mEmptyBlock = nullptr;
            fem::Bearing * mEmptyBearing = nullptr;

            // DOFs used by this proc
            Cell< Dof * > mDOFs;

            // DOFs used by each proc ( master only )
            Cell< Vector< index_t > > mDofIndices;
            Cell< Vector< id_t > > mDofIDs;

            // indices on master matrix
            Cell< Vector< index_t > > mJacobianTable;
            Cell< Vector< index_t > > mDirichletTable;

            // for convective heatflux ( currently not used )
            Cell< Vector< index_t > > mConvectionTable;

            // maps nodes on convection surfaces with indices
            Map< id_t, index_t >      mConvectionMap;

            // data container for convective data
            // contains e.g. the heatflux or external pressure
            // length: ( num nodes ) * ( num dofs per node )
            Vector< real > mConvection ;

            // vector for volume loads, eg. imposed current
            Vector< real > mVolumeLoads ;

            // number of convection nodes on this proc
            index_t mNumberOfConvectionNodes = 0;

            index_t mNumberOfFreeDofs;
            index_t mNumberOfFixedDofs;

            // values for this proc, identical to mMyNumberOfFreeDofs for Master
            index_t mMyNumberOfFreeDofs;
            index_t mMyNumberOfFixedDofs;

            // master only
            Vector< index_t > mNumberOfFreeDofsPerProc;
            Vector< index_t > mNumberOfFixedDofsPerProc;

            SpMatrix * mJacobian        = nullptr;
            SpMatrix * mDirichletMatrix = nullptr;

            //SpMatrix * mBearingMatrix = nullptr;

            // needed for newton raphson
            Vector< real > mFieldValues;

            // left hand side ( field values or deltas )
            Vector< real > mLhsVector;
            Matrix< real > mLhsMatrix;

            // right hand side
            Vector< real > mRhsVector;
            Matrix< real > mRhsMatrix;

            /**
             * this list contains the nodes as owned by each proc
             * exists only on master
             */
             Cell< Vector< index_t > > mNodeOwnerList ;

             // only non-master
             index_t mMyNumberOfOwnedNodes ;

             // offset for computing edge dof IDs
             id_t mEdgeIdOffset = gNoID ;

             // special purpose data for linear to higher projection
             Cell< Vector< index_t > > mAllCornerNodeIndices ;
             Vector< index_t > mMyCornerNodeIndices ;


             Cell< Vector< index_t > > mAllNonCornerNodeIndices ;
             Cell< mesh::Node * > mMyNonCornerNodes ;

             Vector< index_t > mMyNonCornerNodeIndices ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Field(
                    Kernel * aParent,
                    const Vector< id_t > & aBlockIDs,
                    const uint             aNumberOfDOFsPerNode,
                    const uint             aNumberOfDOFsPerEdge=0,
                    const bool             aEnforceLinear = false );

//------------------------------------------------------------------------------

            ~Field();

//------------------------------------------------------------------------------


             void
             init_dofs();

//------------------------------------------------------------------------------

            /**
             * initialize the dof values with the values from
             * the field as defined in iwg
             */
             void
             init_dof_values() ;

//------------------------------------------------------------------------------

            /**
             * update dofs sideset wise
             */
            void
            set_boundary_conditions();

//------------------------------------------------------------------------------

            /**
             * set the solver type for this field
             */
            void
            set_solver( const SolverType & aSolver );

//------------------------------------------------------------------------------

            /**
             * expose the solver
             */
            Solver *
            solver();

//------------------------------------------------------------------------------

            /**
             * return a specific dof
             */
            Dof *
            dof( const id_t aID );

//------------------------------------------------------------------------------

             /**
              * tell if that dof exists on current proc
              * @param aID
              * @return
              */
            bool
            dof_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------

            /**
             * tell if that block exists on current proc
             * @param aID
             * @return
             */
            bool
            block_exists( const id_t & aID ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_node() const;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_edge() const;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Node * aNode, const uint aDofType ) const;

//------------------------------------------------------------------------------

            id_t
            calculate_dof_id( const mesh::Edge * aEdge, const uint aDofType ) const;

//------------------------------------------------------------------------------

            /**
             * check of a FEM sideset exists on this proc
             */
            bool
            sideset_exists( const id_t & aID ) const ;

//------------------------------------------------------------------------------
             /**
              * expose a sideset
              */
            SideSet *
            sideset( const id_t aID );

//------------------------------------------------------------------------------

            /**
             * expose a block
             */
            Block *
            block( const id_t aID );

//------------------------------------------------------------------------------

            /**
             * expose block container
             */
            Cell< Block * > &
            blocks();

//------------------------------------------------------------------------------

            /**
             * expose a bearing
             */
            Bearing *
            bearing( const id_t & aID );

//------------------------------------------------------------------------------

            void
            set_integrated_weak_governing_equation( IWG * aIWG );

//------------------------------------------------------------------------------

            void
            initialize_jacobian();

//------------------------------------------------------------------------------

            void
            compute_jacobian( bool aReset=true );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs();

//------------------------------------------------------------------------------

            /**
             * expose the jacobian
             */
             inline SpMatrix &
             jacobian()
             {
                return * mJacobian ;
             }

//------------------------------------------------------------------------------

            /**
             * expose the dirichlet matrix
             */
            inline SpMatrix &
            dirichlet()
            {
                return * mDirichletMatrix ;
            }

//------------------------------------------------------------------------------

             /**
              * expose the rhs vector
               */
             inline Vector< real > &
             rhs_vector()
             {
                 return mRhsVector ;
             }

//------------------------------------------------------------------------------

            /**
             * expose the rhs matrix
             */
            inline Matrix< real > &
            rhs_matrix()
            {
                return mRhsMatrix ;
            }

//------------------------------------------------------------------------------

            /**
             * expose the lhs vector
              */
            inline Vector< real > &
            lhs_vector()
            {
                return mLhsVector ;
            }

//------------------------------------------------------------------------------

            /**
             * expose the lhs matrix
             */
            inline Matrix< real > &
            lhs_matrix()
            {
                return mLhsMatrix ;
            }

//-----------------------------------------------------------------------------

            inline bool
            enforce_linear_interpolation() const
            {
                return mEnforceLinear && mMesh->max_element_order() > 1 ;
            }

//------------------------------------------------------------------------------

            void
            solve();

//------------------------------------------------------------------------------

            real
            residual( const uint aIteration=1 );

//------------------------------------------------------------------------------

            void
            compute_rhs();

//------------------------------------------------------------------------------

            void
            free_all_dofs();

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

            // special function for tatcad ( old )
            // void
            // compute_convection_old();

//------------------------------------------------------------------------------

            // special function to compute convective term
            void
            compute_surface_loads();

//------------------------------------------------------------------------------

             /**
              * special table needed for convective heatflux
              * ( currently not in use, planned for parallelization of dot q
              *  computation )
              */
            void
            create_convection_table();

//------------------------------------------------------------------------------

            /**
             * initialize the value of all dofs to a specific uniform number
             */
             void
             set_all_dofs( const real aValue );

//------------------------------------------------------------------------------

            void
            distribute_fields( const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

            /**
             * resets the jacobian and the dirichlet matrix
             */
             void
             reset_jacobian();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_fields();

//------------------------------------------------------------------------------

            void
            create_blockmap();

//------------------------------------------------------------------------------

            void
            create_dofs();

//------------------------------------------------------------------------------

            void
            create_dof_table();

//------------------------------------------------------------------------------

            void
            collect_dof_ids_on_master( Vector< id_t > & aIDs );

//------------------------------------------------------------------------------

            void
            create_assembly_table();

//------------------------------------------------------------------------------

            void
            create_blocks();

//------------------------------------------------------------------------------

            void
            create_sidesets();

//------------------------------------------------------------------------------

            void
            create_bearings();

//------------------------------------------------------------------------------

            /**
             * compute the indices for fixed and free dofs in the sparse matrix
             */
            void
            compute_dof_indices();

//------------------------------------------------------------------------------

            void
            flag_all_dofs();

//------------------------------------------------------------------------------

            void
            flag_used_dofs();

//------------------------------------------------------------------------------

            void
            unflag_all_dofs();

//------------------------------------------------------------------------------

            void
            synchronize_jacobian();

//------------------------------------------------------------------------------

            void
            collect_rhs_vector();

//------------------------------------------------------------------------------

            void
            collect_vector( Vector< real > & aVector );

//------------------------------------------------------------------------------

            void
            synchronize_rhs_matrix();

//------------------------------------------------------------------------------

            /**
             * synchronize boundary conditions for dofs
             */
            void
            synchronize_dirichlet_bcs();

//------------------------------------------------------------------------------

            /**
             * take the BCs for alpha and Tinf from master
             * and send all to the others
             */
             void
             distribute_alpha_bcs();

//------------------------------------------------------------------------------

            void
            create_alpha_fields();

//------------------------------------------------------------------------------

            void
            compute_rhs_vector();

//------------------------------------------------------------------------------

            void
            compute_rhs_matrix();

//------------------------------------------------------------------------------

            /**
             * currently works only for scalar field
             */
            void
            detect_wetted_sidesets();

//------------------------------------------------------------------------------

            /**
             * currently only makes sense for scalar field
             */
            void
            count_wetted_nodes();

//------------------------------------------------------------------------------

            void
            add_to_system_matrices(
                          Element         * aElement,
                    const Matrix< real >  & aJel,
                    const Vector< real >  & aBel,
                    const uint            & aN,
                          SpMatrix        & aJ,
                          SpMatrix        & aD,
                          Vector < real > & aB ) ;

//------------------------------------------------------------------------------

            void
            collect_field( const string & aLabel );

//------------------------------------------------------------------------------

            void
            collect_fields( const Cell< string > & aLabels );

//------------------------------------------------------------------------------

            void
            collect_node_owners();

//------------------------------------------------------------------------------

            void
            compute_node_based_adjacency(
                    Cell< graph::Vertex *  > & aDOFs,
                    Cell< Dof * > & aAdjacency,
                    const bool aFlag );

//------------------------------------------------------------------------------

            void
            compute_edge_based_adjacency( Cell< graph::Vertex * > & aDOFs,
                                          Cell< Dof * > & aAdjacency,
                                          const bool aFlag );
//------------------------------------------------------------------------------

            void
            compute_node_and_edge_based_adjacency( Cell< graph::Vertex *  > & aDOFs,
                                                   Cell< Dof * > & aAdjacency,
                                                   const bool aFlag );

//------------------------------------------------------------------------------

            void
            flag_nodes();

//------------------------------------------------------------------------------

            void
            flag_edges();

//------------------------------------------------------------------------------

            void
            get_field_indices( Vector< index_t > & aFieldIndices, const EntityType aEntityType );

//------------------------------------------------------------------------------

            void
            update_field_indices();

//------------------------------------------------------------------------------

            /**
             * special function to project linear solution onto hiher order mesh
             */
             void
             initialize_linear_projection_lists();

//------------------------------------------------------------------------------

            /**
             * special function to send data from corner nodes to other procs
             */
             void
             communicate_corner_node_data( const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

            /**
             * special function to receive data from non corner nodes from other procs
             */
            void
            communicate_noncorner_node_data(
                    const Cell< string > & aFieldLabels,
                          Matrix< real > & aData );

//------------------------------------------------------------------------------
            void
            project_linear_field_to_higher_mesh(
                    const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

            void
            reset_convection_vector();

//------------------------------------------------------------------------------

            /**
             * adds computed matrix to system matrix
             * @param aElement
             * @param aJacobian
             */
            void
            assemble_jacobian( Element * aElement, Matrix< real > & aJacobian );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline Solver *
        Field::solver()
        {
            return mSolver ;
        }

//------------------------------------------------------------------------------

        inline Cell< Block * > &
        Field::blocks()
        {
            return mBlocks;
        }

//------------------------------------------------------------------------------

        inline Dof *
        Field::dof( const id_t aID )
        {
            return mDOFs( mDofMap( aID ) );
        }

//------------------------------------------------------------------------------

        inline bool
        Field::dof_exists( const id_t aID ) const
        {
            return mDofMap.key_exists( aID );
        }

//------------------------------------------------------------------------------

        inline bool
        Field::block_exists( const id_t & aID ) const
        {
            return mBlockMap.key_exists( aID );
        }

//------------------------------------------------------------------------------

        inline uint
        Field::number_of_dofs_per_node() const
        {
            return mNumberOfDOFsPerNode;
        }

//------------------------------------------------------------------------------

        inline uint
        Field::number_of_dofs_per_edge() const
        {
            return mNumberOfDOFsPerEdge;
        }

//------------------------------------------------------------------------------

        inline id_t
        Field::calculate_dof_id( const mesh::Node * aNode , const uint aDofType ) const
        {
            return mNumberOfDOFsPerNode * ( aNode->id() -1 ) + aDofType + 1;
        }

//------------------------------------------------------------------------------

        inline id_t
        Field::calculate_dof_id( const mesh::Edge * aEdge, const uint aDofType ) const
        {
            return mEdgeIdOffset +  mNumberOfDOFsPerEdge * ( aEdge->id() -1 ) + aDofType + 1;
        }

//------------------------------------------------------------------------------

        inline uint
        Field::sideset_integration_order() const
        {
            return mSideSetIntegrationOrder;
        }

//------------------------------------------------------------------------------

        inline IntegrationScheme
        Field::integration_scheme() const
        {
            return mIntegrationScheme ;
        }

//------------------------------------------------------------------------------

        inline uint
        Field::block_integration_order() const
        {
            return mBlockIntegrationOrder;
        }

//------------------------------------------------------------------------------

        inline void
        Field::assemble_jacobian( Element * aElement, Matrix< real > & aJacobian )
        {
            SpMatrix & J       =  *mJacobian;
            SpMatrix & D       =  *mDirichletMatrix;

            // get dimension of element Jacobian
            uint tN = aJacobian.n_rows();

            // add element jacobian to master
            for ( uint i = 0; i < tN; ++i )
            {
                Dof * tRow = aElement->dof( i );

                if ( !tRow->is_fixed() )
                {
                    for ( uint j = 0; j < tN; ++j )
                    {
                        Dof * tCol = aElement->dof( j );

                        if ( !tCol->is_fixed() )
                        {
                            J( tRow->index(), tCol->index()) += aJacobian( i, j );
                        }
                        else
                        {
                            D( tRow->index(), tCol->index()) -= aJacobian( i, j );
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_FIELD_HPP
