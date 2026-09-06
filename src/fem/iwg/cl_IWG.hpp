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

#ifndef BELFEM_CL_IWG_HPP
#define BELFEM_CL_IWG_HPP

#include "typedefs.hpp"
#include "cl_Bitset.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_Material.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"
#include "en_IWG_SideSetDofLinkMode.hpp"
#include "en_DomainType.hpp"
#include "en_FEM_GroupActivationMode.hpp"
#include "cl_TimestepMatrices.hpp"

#define BELFEM_MAX_DOFTYPES  32
#define BELFEM_MAX_NUMPROCS  128
//#define BELFEM_FERROAIR_ENRICHED

namespace belfem
{
    class Mesh;

    enum class DofMode
    {
        AllBlocksEqual = 0,
        BlockSpecific  = 1,
        UNDEFINED      = 2
    };



    namespace fem
    {
        class Dof ;
        class DofManagerBase;
        class Group;
        class Block ;
        class SideSet ;
        class Element;
        class BoundaryCondition ;
        class Calculator ;

//------------------------------------------------------------------------------

        /**
          * the dof table contains the dofs per block and sideset
          */
        struct DofTable
        {
            Vector< index_t > Node ;
            Vector< index_t > Edge ;
            Vector< index_t > Face ;
            Vector< index_t > Cell ;
            Vector< index_t > Lambda ;
        };

//------------------------------------------------------------------------------

        /**
         * Prototype for
         * Integrator of Weak Form Governing Equation
         *
         * @ingroup grp_fem_iwg
         * @see @ref fem_iwg_iwg_usage_guide
         */
        class IWG
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // rank of this proc
            const proc_t mCommRank ;

            // relaxation parameter for Newton-Raphson
            real mOmega     = 0.9;

            // penalty factor for weak BC
            Vector< real > mPenalty ;

            //! flag telling if we have been initialized
            bool mIsInitialized = false ;

            //! flag telling if we want to use bubbles on sidesets
            //! this works only if air elements are master
            bool mEnrichSideSets = false ;

            //! value for timeloop, needed if we don't want to write every timestep
            uint mTimeLoop = 0 ;

            bool mComputeJacobianOnBlock   = true ;

            //! no writer anywhere in the tree, by design: the sideset branch
            //! of DofManager::compute_jacobian is dormant capability for a
            //! future problem where surface matrices matter — not needed for
            //! the maxwell-thermal problem. Any future writer must also give
            //! that loop the is_active() guard it alone lacks
            bool mComputeJacobianOnSideset = false ;

            //! this enum tells which equation object is used
            const IwgType    mType;

            //! this enum tells which dimensionality we have
            const ModelDimensionality mDimensionality;

            //! this enum tells if we have to perform a Newton-Raphson
            const IwgMode    mMode;


            //! needed for the solver
            SymmetryMode mSymmetryMode ;

            //! except maxwell, most IWGs are AllBlocksEqual
            const DofMode    mDofMode ;

            //! mode how sidesets are linked
            const SideSetDofLinkMode mSideSetDofLinkMode ;

            //! tells which solver algorithm is to be used
            //! per default, the IWG sets Newton Raphson
            SolverAlgorithm     mSolverAlgorithm = SolverAlgorithm::UNDEFINED ;

            Mesh              * mMesh     = nullptr;
            DofManagerBase    * mField    = nullptr;
            Group             * mGroup    = nullptr;
            Calculator        * mCalc     = nullptr ;

            const Material    * mMaterial = nullptr;


            uint mNumberOfDofsPerNode = 0 ;
            uint mNumberOfDofsPerEdge = 0 ;  // counts dofs on LHS
            uint mNumberOfDofsPerFace = 0 ;  // counts dofs on LHS
            uint mNumberOfRhsDofsPerEdge = 0 ;  // counts dofs on RHS, needed for L2
            uint mNumberOfRhsDofsPerFace = 0 ;  // counts dofs on RHS, needed for L2
            uint mNumberOfThinShellLayers = 0 ;

            // dimension for N. Sentinel-initialized like mNumberOfNodesPerElement
            // below: derived constructors set this, but not all of them do it
            // before something reads it, and an indeterminate uint here reads
            // as a plausible dimension rather than as garbage
            uint mNumberOfSpatialDimensions = BELFEM_UINT_MAX ;

            // dimension for B
            uint mNumberOfDerivativeDimensions = 0 ;

            // matrix size
            uint   mNumberOfNodesPerElement = BELFEM_UINT_MAX ;
            uint   mNumberOfNodesPerMaster  = BELFEM_UINT_MAX ;
            uint   mNumberOfNodesPerSlave   = BELFEM_UINT_MAX ;

            uint   mNumberOfEdgesPerElement = 0 ;
            uint   mNumberOfFacesPerElement = 0 ;

            //! set by link_to_group(), not by the constructor -- the sentinel
            //! makes a premature read deterministic instead of undefined
            uint   mNumberOfDofsPerElement = BELFEM_UINT_MAX ;

            uint   mNumberOfNodeDofsPerElement = 0 ;

            // in this context, face dofs also count as edge dofs
            uint   mNumberOfEdgeDofsPerElement = 0 ;

            // needed for L2 projection
            uint   mNumberOfRhsEdgeDofsPerElement = 0 ;

            // label of DOF fields, has dimension of mNumberOfDofsPerNode
            Cell< string > mDofFields ;

            // label of Flux fields, has dimension of mNumberOfDofsPerNode
            Cell< string > mFluxFields ;

            // special purpose if RHS is a matrix
            Cell< string > mTensorFields ;

            // All the other fields
            Cell< string > mOtherFields ;

            Cell< string > mAllFields ;

            // list of fields that are not to be saved
            Cell< string > mHiddenFields ;

            // number of columns for rhs side
            uint   mNumberOfRhsCols = 1;

            // sideset is needed for convective terms, these are the FEM ones
            Vector< id_t > mWettedSidesets;

            // timestep
            real mDeltaTime = 1.0 ;

            // vector with indices for nodes on wetted surfaces
            Vector< index_t > mNodesOnWettedSidesets ;

            // alpha is a special boundary condition for convective flow
            bool mHasConvection = false ;

            // list of selected block ids
            Vector< id_t > mBlockIDs;

            // list of selected sideset IDs
            Vector< id_t > mSideSetIDs;

            // lookup table for block indices
            Map< id_t, index_t > mBlockIndices ;

            // lookup table for sideset indices
            Map< id_t, index_t > mSideSetIndices ;

            // map for dof types
            Map< string, uint > mDofTypeMap ;
            Map< uint, uint >   mDofFieldMap ;

            // default dof types, assumes all dofs are on all selected blocks
            Vector< index_t > mDefaultDofTypes ;

            // this contains the entity types for the dofs
            Vector< index_t > mDofEntityTypes ;

            // table containing elements per block and entity
            Cell< Vector < index_t > > mDofsPerBlock ;

            // table containing elements per sideset and entity
            Cell< Vector < index_t > > mDofsPerSideSet ;

            //! contains the block dofs
            Cell< DofTable * > mBlockDofs ;

            //! contains the sideset dofs
            Cell< DofTable * > mSideSetDofs ;

            //! contains the sideset dofs that sit only on the sideset,
            //! but not on master or slave
            Cell< DofTable * > mSideSetOnlyDofs ;

            Vector< index_t > mEdgeFieldIndices ;
            Vector< index_t > mFaceFieldIndices ;

            // non-node DOFs need to be handled differently
            // the multiplicity depends on the order of the edge element
            // we assume linear elements for now
            index_t mEdgeDofMultiplicity   = 0 ;
            index_t mFaceDofMultiplicity   = 0 ;
            index_t mCellDofMultiplicity   = 0 ;
            index_t mLambdaDofMultiplicity = 0 ;

            Map< string, uint > mDofMap ;

            //! needed for debug output
            Cell< string > mDofLabels ;

            //! this list contains the sidesets that are shells.
            //! needed by dof manager. Must be populated by
            //! set_sidesets of child class

            //! links block IDs with the designated types
            //! must be populated by child class
            Map< id_t, DomainType > mBlockTypes ;

            //! links sideset IDs with the designated types
            Map< id_t, DomainType > mSideSetTypes ;

            //! links sideset IDs with the designated types
            Map< id_t, DomainType > mSideSetSubTypes ;

            //! maps domain types to activation modes for blocks
            Map< DomainType, GroupActivationMode > mBlockActivationModes ;

            //! maps domain types to activation modes for sidesets
            Map< DomainType, GroupActivationMode > mSideSetActivationModes ;

            //! normal, if this is a 2d problem, todo: delete
            Vector< real > mNormal2D = { 0., 0. };

            //! normal, if this is a 3d problem, todo: delete
            Vector< real > mNormal3D = { 0., 0., 0. };

            InterpolationType mInterpolationType = InterpolationType::LAGRANGE ;

            //! abstract nodes needed for cut BCs in MAXWELL
            Cell< mesh::Node * > mAbstractNodes ;

            //! nodes that belong to no active element but still carry a phi dof ( see Mesh::orphaned_nodes() )
            Cell< mesh::Node * > mOrphanedNodes ;

            //! list with special dofs that sit on abstract nodes
            Cell< Dof * > mAbstractNodeDofs ;

            DofTable mAbstractNodeTable ;
            index_t mAbstractDofType = gNoIndex ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // The symmetry default is UNSYMMETRIC on purpose, and must stay
            // that way. BELFEM assembles and stores the FULL matrix, while
            // MUMPS with SYM != 0 ( PositiveDefiniteSymmetric = SYM 1,
            // GeneralSymmetric = SYM 2 ) expects only ONE half of it. Either
            // half will do -- the user guide 5.5.1 §5.2.2.1 accepts the lower
            // or the upper triangle including the diagonal -- but only one,
            // because "duplicate entries are summed" applies to the symmetric
            // pair: "if both aij and aji are provided, they will be summed".
            // Handing MUMPS the full matrix under SYM != 0 therefore
            // factorizes A with every OFF-DIAGONAL DOUBLED and the diagonal
            // counted once, which is a different operator. Until
            // 2026-08-29 it did so SILENTLY: measured on the eigen path,
            // ||Ax-b||/||b|| = 7.4e16 and a "converged" lambda_min of -2.5e-19
            // against a true 2.46e-6, with the eigensolver reporting success.
            // cl_SolverMUMPS.cpp now refuses SYM != 0 with an always-active
            // error, so that particular failure is loud rather than silent --
            // but the default belongs here regardless, because this is where
            // an IWG that never thinks about symmetry gets its answer.
            //
            // Scope, deliberately narrow: the mode set here travels to the
            // solver on the ordinary DofManager solve path
            // ( cl_FEM_DofMgr_SolverData.cpp, set_symmetry_mode from the
            // active IWG ). It is NOT universal -- the eigenvalue
            // shift-invert path overrides it with Unsymmetric of its own
            // accord, and the triangle requirement is MUMPS's, not a property
            // of every backend that reads SymmetryMode.
            //
            // Symmetric must therefore be an explicit, deliberate choice by a
            // caller that has arranged to supply one half -- never something an
            // IWG inherits by forgetting to pass an argument. That is not
            // hypothetical: deriving straight from this class is the DOCUMENTED
            // extension point for a user writing their own IWG
            // ( dof_manager_usage_guide.md, "class IWG_Custom : public IWG" ),
            // and every in-tree IWG happens to reach Unsymmetric only through
            // IWG_Timestep's own default. A user following the guide never
            // passes this argument, so before 2026-08-29 the documented path
            // handed them SYM = 1 and a silently wrong factorization.
            // This class OWNS the DofTable* held in mBlockDofs, mSideSetDofs
            // and mSideSetOnlyDofs -- allocated with new, deleted in the
            // destructor. The implicit copy would shallow-copy those pointers
            // and double-free them, so copying and moving are deleted rather
            // than left to the compiler. Deep-copy semantics are not wanted:
            // an IWG is created once, linked, and used in place.
            IWG( const IWG & ) = delete ;
            IWG & operator=( const IWG & ) = delete ;
            IWG( IWG && ) = delete ;
            IWG & operator=( IWG && ) = delete ;

            IWG( const IwgType aType,
                 const ModelDimensionality aDimensionality,
                 const IwgMode aMode=IwgMode::Iterative,
                 const SymmetryMode aSymmetryMode=SymmetryMode::Unsymmetric,
                 const DofMode      aDofMode=DofMode::AllBlocksEqual,
                 const SideSetDofLinkMode aSideSetDofLinkMode=SideSetDofLinkMode::FacetOnly );

//------------------------------------------------------------------------------

            virtual ~IWG() ;

//------------------------------------------------------------------------------

            /**
             * the blocks are selected over the IWG object.
             * the equation object knows how many DOFs sit on each entity
             *  per block
             */
            void
            select_blocks( const Vector< id_t > & aBlockIDs );

            void
            select_block( const id_t aBlockID );

//------------------------------------------------------------------------------

            /**
             * the sidesets are selected over the IWG object.
             * the equation object knows how many DOFs sit on each entity
             *  per sideset
             */
            void
            select_sidesets( const Vector< id_t > & aSidesetIDs );

//------------------------------------------------------------------------------

            /**
             * returns the list of selected blocks
             */
             const Vector< id_t > &
             selected_blocks() const ;

//------------------------------------------------------------------------------

            /**
             * returns the list of selected sidesets
             */
            const Vector< id_t > &
            selected_sidesets() const ;

//------------------------------------------------------------------------------

            virtual void
            compute_jacobian(
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            virtual void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            // convection term, eg for heat load
            virtual void
            compute_convection(
                    Element        * aElement,
                    Vector< real > & aConvection );

//------------------------------------------------------------------------------

            // convection term with alpha boundary condition
            virtual void
            compute_alpha_boundary_condition(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            virtual void
            compute_rhs(
                    Element        * aElement,
                    Vector< real > & aRHS );


//------------------------------------------------------------------------------

            virtual void
            compute_rhs(
                    Element        * aElement,
                    Matrix< real > & aRHS );

//------------------------------------------------------------------------------

            virtual void
            link_to_group( Group * aGroup );

//------------------------------------------------------------------------------

            /**
             * return the names of the potential fields
             */
            const Cell< string > &
            dof_fields() const ;

//------------------------------------------------------------------------------

            /**
             * return the names of the flux fields
             */
            const Cell< string > &
            flux_fields() const;

//------------------------------------------------------------------------------

            /**
             * special purpose if RHS is a matrix
             */
            const Cell< string > &
            tensor_fields() const;

//------------------------------------------------------------------------------

            /**
             * return the names of the other fields
             */
            const Cell< string > &
            other_fields() const;

//------------------------------------------------------------------------------

            /**
             * return the names of all fields
             */
            const Cell< string > &
            all_fields() const ;

//------------------------------------------------------------------------------

            /**
             * return the names of a specific field
             */
            const string &
            field( const index_t aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * return the number of connected fiends
             */
            index_t
            number_of_fields() const ;

//------------------------------------------------------------------------------

            /**
             * DEPRECATED!
             * @return
             */
            uint
            number_of_dofs_per_node() const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_node( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_node_on_sideset( const id_t aSideSetID ) const ;

//------------------------------------------------------------------------------


            uint
            number_of_dofs_per_edge( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_face( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_cell( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_dofs_per_element( Block * aBlock ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_nodes_per_element( SideSet * aSideSet ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_edges_per_element( SideSet * aSideSet ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_faces_per_element( SideSet * aSideSet ) const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_dofs_per_element( SideSet * aSideSet ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_lambda_dofs( const id_t aSideSetID  ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_node( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_edge( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_face( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_cell( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            // Vector< index_t > &
            // lambda_dofs( const id_t aSideSetID );

//------------------------------------------------------------------------------

            const Vector< index_t > &
            lambda_dofs( const id_t aSideSetID ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_node_on_sideset( const id_t aSideSetID, const bool aSideSetOnly = false ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_edge_on_sideset( const id_t aSideSetID, const bool aSideSetOnly = false ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_face_on_sideset( const id_t aSideSetID, const bool aSideSetOnly = false ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            dofs_per_cell_on_sideset( const id_t , const bool aSideSetOnly = false ) const ;

//------------------------------------------------------------------------------

            const Vector< index_t  > &
            dof_entity_types() const ;

//------------------------------------------------------------------------------

            const Vector< index_t  > &
            default_dof_types() const ;

//------------------------------------------------------------------------------

            /**
             * tells how many dofs sit on one edge
             */
             index_t
             edge_multiplicity() const ;

//------------------------------------------------------------------------------

            /**
             * tells how many dofs sit on one face
             */
            index_t
            face_multiplicity() const ;

//------------------------------------------------------------------------------

            /**
             * tells how many dofs sit on one cell
             */
            index_t
            cell_multiplicity() const ;

//------------------------------------------------------------------------------

            index_t
            edge_field_index( const index_t aDofType ) const;

//------------------------------------------------------------------------------

            index_t
            face_field_index( const index_t aDofType ) const;

//------------------------------------------------------------------------------

            /**
             * intended to be used for interface dof creation
             * @param aSidesetID
             * @return
             */
            inline const Vector< index_t > &
            dofs_per_sideset( const id_t aSidesetID ) const ;

//------------------------------------------------------------------------------

            uint
            num_rhs_cols() const;

//------------------------------------------------------------------------------

            void
            virtual set_field( DofManagerBase * aField );

//------------------------------------------------------------------------------

            /**
             * return the type of this IWG
             */
             IwgType
             type() const ;

//------------------------------------------------------------------------------

            /**
             * return the dimensionality of this IWG
             */
            ModelDimensionality
            model_dimensionality() const ;

//------------------------------------------------------------------------------

            /**
             * set the interpolation type of node elements
             */
            void
            set_interpolation_type( const InterpolationType aType );

//------------------------------------------------------------------------------

            /**
             * return the interpolation type of node elements
             */
            InterpolationType
            interpolation_type() const ;

//------------------------------------------------------------------------------

            /**
             * return the calculation mode of this IWG
             */
            IwgMode
            mode() const ;

//------------------------------------------------------------------------------

            /**
             * get relaxation parameter
             */
             real
             omega() const;

//------------------------------------------------------------------------------

            /**
              * set the relaxation parameter
              */
            void
            set_omega( const real & aOmega );

//------------------------------------------------------------------------------

            /**
              * set the penalty parameter
              */
            void
            set_penalty( const real aPsi, const uint aIndex=0 );

//------------------------------------------------------------------------------

            /**
             * get penalty parameter
             */
            real
            penalty( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * flag telling if matrices are symmetric
             */
             SymmetryMode
             symmetry_mode() const;

//------------------------------------------------------------------------------

            /**
             * special function, must be called before init->jacobian() is called
             */
            void
            set_num_rhs_cols( const uint & aNumRhsCols );

//------------------------------------------------------------------------------

            virtual void
            set_wetted_sidesets( const Vector< id_t > & aSideSets );

//------------------------------------------------------------------------------

            const Vector< id_t >  &
            wetted_sidesets() const ;

//------------------------------------------------------------------------------

            // timestep, if this is a transient problem
            real &
            delta_time() ;

//---------------------------------------------------------------------------------

            uint &
            time_loop() ;

//------------------------------------------------------------------------------

            /**
             * add additional non-dof fields to mOtherFields and mAllFields; a label "alpha" also sets mHasConvection
             * @param aFieldLabels
             */
            void
            add_fields( const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

            /**
             * tells if this field has an alpha boundary condition
             */
             bool
             has_convection() const ;

//------------------------------------------------------------------------------

            /**
             * special function to compute the boundary flux in magnetics
             * impose zero as weak BC
             */
            virtual void
            compute_boundary_flux_matrix(
                    Element        * aElement,
                    const uint       aDirection,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            const Matrix< real > &
            N( const uint & aIntegrationPoint );

//------------------------------------------------------------------------------

            void
            set_blocks(
                    const Vector< id_t >      & aBlockIDs,
                    const Cell< DomainType >  & aBlockTypes ) ;

//------------------------------------------------------------------------------

            void
            set_sidesets(
                const Vector< id_t > & aSideSetIDs,
                const Cell< DomainType > & aSideSetTypes ) ;

//------------------------------------------------------------------------------

            DomainType
            block_type( const id_t aID ) const ;

//------------------------------------------------------------------------------

            DomainType
            sideset_type( const id_t aID ) const ;

//------------------------------------------------------------------------------

            /**
             * returns the activation mode for a block based on its domain type
             */
            GroupActivationMode
            block_activation_mode( const DomainType aType ) const ;

//------------------------------------------------------------------------------

            /**
             * returns the activation mode for a sideset based on its domain type
             */
            GroupActivationMode
            sideset_activation_mode( const DomainType aType ) const ;

//------------------------------------------------------------------------------

            virtual void
            initialize();

//------------------------------------------------------------------------------

            virtual void
            initialize( const IwgType aType );

//------------------------------------------------------------------------------

            bool
            is_initialized() const ;

            inline bool
            enrich_sidesets() const ;

//------------------------------------------------------------------------------

            /**
             * This routine makes sure that the mesh fulfills the requirements.
             * Unless those are not explicitly specified by the IWG,
             * this routine only checks if the mesh is not a nullpointer.
             */
            virtual int
            check_mesh( Mesh * aMesh, const proc_t aMasterRank=0 );

//------------------------------------------------------------------------------

            /**
             * hides fields that are not to be saved to exodus
             */
             void
             hide_fields_from_exodus( Mesh * aMesh );

//-----------------------------------------------------------------------------

            /**
             * flag telling if the IWG has edge dofs
             */
             virtual bool
             has_edge_dofs() const ;

//-----------------------------------------------------------------------------

             /**
              * if this number is > 0, we have langrande dofs such as contact
              */
             uint
             lambda_multiplicity() const ;

//-----------------------------------------------------------------------------

            /**
             * return the dof linking mode for sideset
             */
             SideSetDofLinkMode
             sideset_dof_link_mode() const ;

//-----------------------------------------------------------------------------

            /**
             * sets the mode for the solver algorithm
             */
             void
             set_algorithm( const SolverAlgorithm aAlgorithm  );

//-----------------------------------------------------------------------------

             /**
              * returns the solver algorithm
              */
             SolverAlgorithm
             algorithm() const ;

//------------------------------------------------------------------------------

            /**
             * called by main file to copy fields into last timestep
             * eg. T0 = T
             */
            virtual void
            shift_fields();

//------------------------------------------------------------------------------

            /**
             * called by main file to copy fields from last timestep
             * eg. T = T0
             */
            virtual void
            reset_fields();

//------------------------------------------------------------------------------

            void
            collect_node_data(
                    Element        * aElement,
                    Cell< string > & aFieldLabels,
                    Matrix< real > & aData );

//------------------------------------------------------------------------------

            void
            collect_node_data(
                    Element        * aElement,
                    const string   & aFieldLabel );

//------------------------------------------------------------------------------

            void
            collect_node_data(
                    Element        * aElement,
                    const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

            // todo: old function, should be obsolete soon
            void
            collect_node_data(
                    Element        * aElement,
                    const string   & aFieldLabel,
                    Vector< real > & aData );

//------------------------------------------------------------------------------

            // todo: old function, should be obsolete soon
            void
            collect_node_data(
                    Element        * aElement,
                    const string   & aFieldLabel,
                    Vector< real > & aData,
                              uint & aOffset );

//------------------------------------------------------------------------------

            void
            collect_edge_data(
                    Element        * aElement,
                    const string   & aEdgeFieldLabel,
                    Vector< real > & aData );

//------------------------------------------------------------------------------

            void
            collect_edge_data(
                    Element        * aElement,
                    const string   & aEdgeFieldLabel,
                    const string   & aFaceFieldLabel,
                    Vector< real > & aData );

//------------------------------------------------------------------------------

            void
            collect_lambda_data(
                    Element        * aElement,
                    const string   & aFieldLabel,
                    real & aData ) ;

//------------------------------------------------------------------------------

            /**
             * return the type id of a dof
             */
            uint
            doftype( const string & aDofLabel ) const ;

//---------------------------------------------------------------------------------

            /**
             * makes sure that dof list is unique and also creares dofmap
             * must be accessible by Mawell_FieldList as well
             */
            void
            unique_and_rearrange( Cell< string > & aDofs, const bool aMakeMap=false );

//---------------------------------------------------------------------------------

            /**
             * called by dof manager
             */
             bool
             compute_jacobian_on_sideset() const ;

//---------------------------------------------------------------------------------

             /**
              * called by dof manager
              */
             bool
             compute_jacobian_on_block() const ;

//---------------------------------------------------------------------------------

            virtual void
            compute_mkf( Element * aElement );

//------------------------------------------------------------------------------

            virtual TimestepMatrices *
            matrices();

//------------------------------------------------------------------------------

            /**
             * for debugging
             */
            void
            print_dofs( Element * aElement, const bool aLocal=false );

//---------------------------------------------------------------------------------

            /**
             * for debugging
             */
             const string &
             dof_label( const index_t aDofIndex );

//---------------------------------------------------------------------------------
            /**
             * @return returns the calculator object
             */
            Calculator *
            calc();

//------------------------------------------------------------------------------
            // these nodes contain for example currents in maxwell cuts

            void
            set_abstract_nodes( Cell< mesh::Node * > & aNodes );

            void
            set_orphaned_nodes( Cell< mesh::Node * > & aNodes );

            Cell< mesh::Node * > &
            abstract_nodes();

            Cell< mesh::Node * > &
            orphaned_nodes();

            index_t
            abstract_dof_type() const ;

            virtual void
            collect_abstract_node_dofs() ;

//------------------------------------------------------------------------------

            virtual void
            create_custom_vectors_and_matrices( Calculator * aCalc );

            uint
            get_field_index( const uint aDofType ) const ;

//------------------------------------------------------------------------------

            virtual void
            custom_postprocess();

//------------------------------------------------------------------------------

            void
            set_abstract_dof_type( const uint aDofType );



//------------------------------------------------------------------------------

            // set the timestepping method
            virtual void
            set_timestepping_method(
                const EulerMethod aMethod,
                const bool aHaveStiffness=true );

            // return the timestepping method
            virtual EulerMethod
            method() const ;

//------------------------------------------------------------------------------

            virtual uint
            timestepping_order() const ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * Reserved; not defined and not called by the destructor.
             */
            void
            delete_boundary_conditions();

//------------------------------------------------------------------------------

            /**
             * Tidy up memory.Called by destructor.
             */
            void
            delete_block_dof_tables();

//------------------------------------------------------------------------------

            /**
             * Tidy up memory.Called by destructor.
             */
            void
            delete_sideset_dof_tables();

//------------------------------------------------------------------------------
            /**
             * this function assumes that all blocks are treated equally.
             * special purpose IWGs might overload this function
             */
            virtual void
            assign_dofs_per_block( const Vector< id_t > & aBlockIDs );

//------------------------------------------------------------------------------

            /**
             * this function assumes that all sidesets are treated equally.
             * special purpose IWGs might overload this function
             */
            virtual void
            assign_dofs_per_sideset( const Vector< id_t > & aSideSetIDs );

//------------------------------------------------------------------------------

            /**
             * populates the activation mode maps for blocks and sidesets
             * based on domain types. Called during initialization.
             * Default implementation maps DomainType::Inactive to
             * GroupActivationMode::Inactive and DomainType::GeometryOnly to
             * GroupActivationMode::GeometryOnly; every other domain type falls
             * back to GeometryAndDofs in block_activation_mode() /
             * sideset_activation_mode(). Override in derived classes to add entries.
             */
            virtual void
            init_activation_maps();

//------------------------------------------------------------------------------

            virtual void
            allocate_work_matrices( Group * aGroup );

//------------------------------------------------------------------------------

            /**
             * populate the node indices for the nodes that sit on wetted sidesets
             */
            void
            collect_nodes_on_wetted_sitdesets( Mesh * aMesh, const Vector< id_t > & aSideSets );

//------------------------------------------------------------------------------

            void
            count_dofs_per_block();

//------------------------------------------------------------------------------

            void
            count_dofs_per_sideset();

//------------------------------------------------------------------------------

            void
            collect_node_coords( Element * aElement, Matrix< real > & aX );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_doftype_map() ;

//---------------------------------------------------------------------------------

            void
            concatenate_field_lists();

//---------------------------------------------------------------------------------

            void
            create_block_dof_tables( const uint aNumBlocks );

//---------------------------------------------------------------------------------

            void
            create_sideset_dof_tables( const uint aNumSideSets );

//---------------------------------------------------------------------------------

            void
            count_sideset_dofs_per_sideset(
                    Vector< index_t >             & aDofsPerSideSet,
                    DofTable                      * aDofTable,
                    Vector< uint >                & aCount,
                    Bitset< BELFEM_MAX_DOFTYPES > & aBitset,
                    const bool                      aUseBitset );

//---------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        IWG::selected_blocks() const
        {
            return mBlockIDs ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        IWG::selected_sidesets() const
        {
            return mSideSetIDs ;
        }

//------------------------------------------------------------------------------

         inline bool
         IWG::has_convection() const
         {
            return mHasConvection ;
         }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        IWG::dof_fields() const
        {
         return mDofFields ;
        }


//------------------------------------------------------------------------------

        inline const Cell< string > &
        IWG::flux_fields() const
        {
            return mFluxFields ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        IWG::tensor_fields() const
        {
            return mTensorFields ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        IWG::other_fields() const
        {
            return mOtherFields ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        IWG::all_fields() const
        {
            return mAllFields ;
        }

//------------------------------------------------------------------------------

        inline const string &
        IWG::field( const index_t aIndex ) const
        {
            return mAllFields( aIndex );
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::number_of_fields() const
        {
            return mAllFields.size() ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG::omega() const
        {
            return mOmega ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG::penalty( const uint aIndex ) const
        {
            return mPenalty( aIndex ) ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::number_of_dofs_per_node( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Node.length() ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::number_of_dofs_per_node_on_sideset( const id_t aSideSetID ) const
        {
            return mSideSetDofs( mSideSetIndices( aSideSetID ) )->Node.length() ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::number_of_dofs_per_edge( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Edge.length() * mEdgeDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::number_of_dofs_per_face( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Face.length() * mFaceDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::number_of_dofs_per_cell( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Cell.length() * mCellDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_node( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID  ) )->Node;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_edge( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID  ) )->Edge;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_face( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID  ) )->Face;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::edge_multiplicity() const
        {
            return mEdgeDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::face_multiplicity() const
        {
            return mFaceDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::cell_multiplicity() const
        {
            return mCellDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_cell( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID  ) )->Cell;
        }

//------------------------------------------------------------------------------

       // inline Vector< index_t > &
        //IWG::lambda_dofs( const id_t aSideSetID )
        //{
         //   return mSideSetLambdaDofs( mSideSetIndices( aSideSetID ) );
        //}

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_node_on_sideset( const id_t aSideSetID, const bool aSideSetOnly ) const
        {
            return aSideSetOnly ?
                   mSideSetOnlyDofs( mSideSetIndices( aSideSetID ) )->Node :
                   mSideSetDofs( mSideSetIndices( aSideSetID ) )->Node ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_edge_on_sideset( const id_t aSideSetID, const bool aSideSetOnly ) const
        {
            return aSideSetOnly ?
                   mSideSetOnlyDofs( mSideSetIndices( aSideSetID ) )->Edge :
                   mSideSetDofs( mSideSetIndices( aSideSetID ) )->Edge ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_face_on_sideset( const id_t aSideSetID, const bool aSideSetOnly  ) const
        {
            return aSideSetOnly ?
                   mSideSetOnlyDofs( mSideSetIndices( aSideSetID ) )->Face :
                   mSideSetDofs( mSideSetIndices( aSideSetID ) )->Face ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_cell_on_sideset( const id_t aSideSetID, const bool aSideSetOnly  ) const
        {
            return aSideSetOnly ?
                   mSideSetOnlyDofs( mSideSetIndices( aSideSetID ) )->Cell :
                   mSideSetDofs( mSideSetIndices( aSideSetID ) )->Cell ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::lambda_dofs( const id_t aSideSetID ) const
        {
            return mSideSetDofs( mSideSetIndices( aSideSetID ) )->Lambda ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t  > &
        IWG::dof_entity_types() const
        {
            return mDofEntityTypes ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t  > &
        IWG::default_dof_types() const
        {
            return mDefaultDofTypes ;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::edge_field_index( const index_t aDofType ) const
        {
            return mEdgeFieldIndices( aDofType );
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::face_field_index( const index_t aDofType ) const
        {
            return mFaceFieldIndices( aDofType );
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        IWG::dofs_per_sideset( const id_t aSidesetID ) const
        {
            return mDofsPerSideSet( mSideSetIndices( aSidesetID ) );
        }

//------------------------------------------------------------------------------

        inline bool
        IWG::is_initialized() const
        {
            return mIsInitialized ;
        }

//------------------------------------------------------------------------------

        inline bool
        IWG::enrich_sidesets() const
        {
            return mEnrichSideSets ;
        }

//------------------------------------------------------------------------------

        inline int
        IWG::check_mesh( Mesh * aMesh, const proc_t aMasterRank )
        {
            return aMesh == nullptr ? 1 : 0 ;
        }


//------------------------------------------------------------------------------

        inline bool
        IWG::has_edge_dofs() const
        {
            return false ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::lambda_multiplicity() const
        {
            return mLambdaDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline SideSetDofLinkMode
        IWG::sideset_dof_link_mode() const
        {
            return mSideSetDofLinkMode ;
        }

//------------------------------------------------------------------------------

        inline SolverAlgorithm
        IWG::algorithm() const
        {
            return mSolverAlgorithm ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG::doftype( const string & aDofLabel ) const
        {
            if( mDofTypeMap.key_exists( aDofLabel ) )
            {
                return mDofTypeMap( aDofLabel );
            }
            else
            {
                return BELFEM_UINT_MAX ;
            }
        }

//------------------------------------------------------------------------------

        inline bool
        IWG::compute_jacobian_on_sideset() const
        {
            return mComputeJacobianOnSideset ;
        }


//---------------------------------------------------------------------------------

        inline bool
        IWG::compute_jacobian_on_block() const
        {
            return mComputeJacobianOnBlock ;
        }

//---------------------------------------------------------------------------------

        inline const string &
        IWG::dof_label( const index_t aDofIndex )
        {
            return mDofLabels( aDofIndex );
        }

//---------------------------------------------------------------------------------

        inline uint &
        IWG::time_loop()
        {
            return mTimeLoop ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG::set_interpolation_type( const InterpolationType aType )
        {
            mInterpolationType = aType ;
        }

//------------------------------------------------------------------------------

        inline InterpolationType
        IWG::interpolation_type() const
        {
            return mInterpolationType ;
        }

//------------------------------------------------------------------------------

        inline Calculator *
        IWG::calc()
        {
            return mCalc ;
        }

//------------------------------------------------------------------------------

        inline Cell< mesh::Node * > &
        IWG::abstract_nodes()
        {
            return mAbstractNodes ;
        }


        inline index_t
        IWG::abstract_dof_type() const
        {
            return mAbstractDofType ;
        }

        inline Cell< mesh::Node * > &
        IWG::orphaned_nodes()
        {
            return mOrphanedNodes ;
        }

        inline uint
        IWG::get_field_index( const uint aDofType ) const
        {
            return mDofFieldMap( aDofType );
        }

//---------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IWG_HPP
