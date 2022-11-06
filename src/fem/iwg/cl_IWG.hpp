//
// Created by Christian Messe on 04.11.19.
//

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

#define BELFEM_MAX_DOFTYPES  32
#define BELFEM_MAX_NUMPROCS  128
#define BELFEM_FERROAIR_ENRICHED

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
        class DofManagerBase;
        enum class DomainType ;
        class Group;
        class Block ;
        class SideSet ;
        class Element;
        class BoundaryCondition ;

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
         */
        class IWG
        {
            // rank of this proc
            const proc_t mRank ;

            // relaxation parameter for Newton-Raphson
            real mOmega     = 0.9;

            // penalty factor for weak BC
            real mPsi       = 1.0 ;

            // penalty factor for tape smoothing
            real mTau       = 0.0 ;

            // flag telling if we have been initialized
            bool mIsInitialized = false ;

            // value for timeloop, needed if we don't want to write every timestep
            uint mTimeLoop = 0 ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            bool mComputeJacobianOnBlock   = true ;
            bool mComputeJacobianOnSideset = false ;

            //!  boundary conditions
            Cell< BoundaryCondition * > mBoundaryConditions ;

            //! this enum tells which equation object is used
            const IwgType    mType;

            //! this enum tells if we have to perform a Newton-Raphson
            const IwgMode    mMode;


            //! needed for the solver
                SymmetryMode mSymmetryMode ;

            //! except maxwell, most IWGs are AllBlocksEqual
            const DofMode    mDofMode ;

            //! mode how sidesets are linked
            const SideSetDofLinkMode mSideSetDofLinkMode ;


            //! mode how thin shells are linked
                SideSetDofLinkMode mThinShellDofLinkMode = SideSetDofLinkMode::UNDEFINED ;

            //! tells which solver algorithm is to be used
            //! per default, the IWG sets Newton Raphson
            SolverAlgorithm     mSolverAlgorithm = SolverAlgorithm::UNDEFINED ;

            Mesh              * mMesh     = nullptr;
            DofManagerBase    * mField    = nullptr;
            Group             * mGroup    = nullptr;

            const Material * mMaterial = nullptr;


            uint mNumberOfDofsPerNode = 0 ;
            uint mNumberOfDofsPerEdge = 0 ;  // counts dofs on LHS
            uint mNumberOfDofsPerFace = 0 ;  // counts dofs on LHS
            uint mNumberOfRhsDofsPerEdge = 0 ;  // counts dofs on RHS, needed for L2
            uint mNumberOfRhsDofsPerFace = 0 ;  // counts dofs on RHS, needed for L2
            uint mNumberOfThinShellLayers = 0 ;

            // dimension for N
            uint mNumberOfSpatialDimensions;

            // dimension for B
            uint mNumberOfDerivativeDimensions = 0 ;

            // matrix size
            uint   mNumberOfNodesPerElement = BELFEM_UINT_MAX ;
            uint   mNumberOfNodesPerMaster = BELFEM_UINT_MAX ;
            uint   mNumberOfNodesPerSlave = BELFEM_UINT_MAX ;

            uint   mNumberOfEdgesPerElement = 0 ;
            uint   mNumberOfFacesPerElement = 0 ;
            uint   mNumberOfIntegrationPoints = 0 ;

            uint   mNumberOfDofsPerElement;

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
            bool mHasAlpha = false ;

            // list of selected block ids
            Vector< id_t > mBlockIDs;

            // list of selected sideset IDs
            Vector< id_t > mSideSetIDs;

            // master sideset id for ghosts
            id_t mGhostSideSetMaster = 0 ;

            // list of selected ghost sideset IDs
            Vector< id_t > mGhostSideSetIDs;

            // list of selected ghost block IDs
            Vector< id_t > mGhostBlockIDs;

            // lookup table for block indices
            Map< id_t, index_t > mBlockIndices ;

            // lookup table for sideset indices
            Map< id_t, index_t > mSideSetIndices ;

            // map for dof types
            Map< string, uint > mDofTypeMap ;

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

            // normal, if this is a 2d problem
            Vector< real > mNormal2D = { 0., 0., };

            // normal, if this is a 3d problem
            Vector< real > mNormal3D = { 0., 0., 0. };

            InterpolationType mInterpolationType = InterpolationType::LAGRANGE ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG( const IwgType aType,
                 const IwgMode aMode=IwgMode::Iterative,
                 const SymmetryMode aSymmetryMode=SymmetryMode::PositiveDefiniteSymmetric,
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
            select_block( const id_t aBlockIDs );

//------------------------------------------------------------------------------

            /**
             * the blocks are selected over the IWG object.
             * the equation object knows how many DOFs sit on each entity
             *  per block
             */
            void
            select_sidesets( const Vector< id_t > & aSidesetIDs );

//------------------------------------------------------------------------------

            void
            set_thin_shell_link_mode( const SideSetDofLinkMode aLinkMode );

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

            const Vector< id_t > &
            ghost_sidesets() const ;

            const Vector< id_t > &
            ghost_blocks() const ;

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
            number_of_dofs_per_edge( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_face( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            uint
            number_of_dofs_per_cell( const id_t aBlockID ) const ;

//------------------------------------------------------------------------------

            uint
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

            uint
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
            set_field( DofManagerBase * aField );

//------------------------------------------------------------------------------

            /**
             * return the type of this IWG
             */
             IwgType
             type() const ;

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
            set_psi( const real aPsi );

//------------------------------------------------------------------------------

            /**
             * get penalty parameter
             */
            real
            psi() const;

//------------------------------------------------------------------------------

            /**
              * set the penalty parameter
              */
            void
            set_tau( const real aPsi );

//------------------------------------------------------------------------------

            /**
             * get penalty parameter
             */
            real
            tau() const;

//------------------------------------------------------------------------------

            /**
             * flag telling if matrices are symmetric
             */
             SymmetryMode
             symmetry_mode() const;

//------------------------------------------------------------------------------

            /**
             * special function, must be called before init->jacobian() is called
             * @return
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
            const real &
            timestep() const;

//---------------------------------------------------------------------------------

            uint &
            time_loop() ;

//------------------------------------------------------------------------------

            /**
             * add additional fields to mFieldLabels and mSurfaceFieldLabels
             * @param aFieldLabels
             */
            void
            add_fields( const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

            /**
             * tells if this field has an alpha boundary condition
             */
             bool
             has_alpha() const ;

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
            set_sidesets( const Vector< id_t >    & aSideSetIDs,
                          const Cell< DomainType > & aSideSetTypes );

//------------------------------------------------------------------------------

            /**
             *
             * @param aMasterSideSet   : master sideset to copy dofs from
             * @param aGhostSideSetIDs : ids as created by tape roller
             */
            void
            set_ghost_sidesets(
                    const id_t                aMasterSideSet,
                    const Vector< id_t >    & aGhostSideSetIDs );

//------------------------------------------------------------------------------

            /**
             * @param aGhostBlockIDs : ids as created by tape roller
             */
            void
            set_ghost_blocks( const Vector< id_t >    & aGhostBlockIDs );

//------------------------------------------------------------------------------

            /**
             * return the ghost sideset IDs
             */
            const Vector< id_t > &
            ghost_sideset_ids() const ;

//------------------------------------------------------------------------------

            /**
             * return the ghost sideset IDs
             */
            uint
            number_of_ghost_sidesets() const ;

//------------------------------------------------------------------------------

            DomainType
            block_type( const id_t aID ) const ;

//------------------------------------------------------------------------------

            DomainType
            sideset_type( const id_t aID ) const ;

//------------------------------------------------------------------------------

            /**
             * called by dof manager
             */
            virtual void
            initialize();

//------------------------------------------------------------------------------

            bool
            is_initialized() const ;

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
             * return the dof linking mode for thin shells
             */
            SideSetDofLinkMode
            thin_shell_dof_link_mode() const ;

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
             * Adds a boundary condition to the equation.
             * Called by constructor
             */
            virtual void
            add_boundary_condition( BoundaryCondition * aBoundaryCondition );

//------------------------------------------------------------------------------

            /*
             * updates the boundary conditions based on timestep
             */
            void
            compute_boundary_conditions( const real aTimeStamp=0.0 );

//------------------------------------------------------------------------------

            /**
             * called by main file to copy fields into last timestep
             * eg. T0 = T
             */
            virtual void
            shift_fields();

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
                    const string   & aFieldLabel,
                    Vector< real > & aData );

//------------------------------------------------------------------------------

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

            void
            collect_node_data_from_layer(
                    Element        * aElement,
                    const string   & aLabel,
                    const uint       aLayer,
                    Vector< real > & aData );

//------------------------------------------------------------------------------

            void
            collect_edge_data_from_layer(
                    Element        * aElement,
                    const string   & aEdgeFieldLabel,
                    const uint       aLayer,
                    Vector< real > & aData );

//------------------------------------------------------------------------------

            void
            collect_edge_data_from_layer(
                    Element        * aElement,
                    const string   & aEdgeFieldLabel,
                    const string   & aFaceFieldLabel,
                    const uint       aLayer,
                    Vector< real > & aData );

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

//------------------------------------------------------------------------------

            /**
             * for debugging
             */
            void
            print_dofs( Element * aElement );

//---------------------------------------------------------------------------------

            /**
             * for debugging
             */
             const string &
             dof_label( const index_t aDofIndex );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * Tidy up memory.Called by destructor.
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

            virtual void
            allocate_work_matrices( Group    * aGroup );

//------------------------------------------------------------------------------

            /**
             * populate the node indices for the nodes that sit on wetted sidesteds
             * @return
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

            /**
             * computes the normal of the facet.
             * @param aElement
             * @param aIndex
             * @return
             */
            const Vector< real > &
            normal_tri3( Element * aElement, const uint aIndex=0  ) ;

//------------------------------------------------------------------------------
            /**
             * computes the normal of the facet.
             * We assume that the coordinates of the master element
             * have been collected in mGroup->work_Xm();
             * @param aElement
             * @param aIndex
             * @return
             */
            const Vector< real > &
            normal_tri6( Element * aElement, const uint aIndex  ) ;

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
            create_sideset_dof_tables( const uint aNumSideSets, const uint aNumGhostSidesets );

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

        inline const Vector< id_t > &
        IWG::ghost_sidesets() const
        {
            return mGhostSideSetIDs ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        IWG::ghost_blocks() const
        {
            return mGhostBlockIDs ;
        }

//------------------------------------------------------------------------------

         inline bool
         IWG::has_alpha() const
         {
            return mHasAlpha ;
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
        IWG::psi() const
        {
            return mPsi ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG::tau() const
        {
            return mTau ;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::number_of_dofs_per_node( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Node.length() ;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::number_of_dofs_per_edge( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Edge.length() * mEdgeDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline index_t
        IWG::number_of_dofs_per_face( const id_t aBlockID ) const
        {
            return mBlockDofs( mBlockIndices( aBlockID ) )->Face.length() * mFaceDofMultiplicity ;
        }

//------------------------------------------------------------------------------

        inline index_t
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

        inline SideSetDofLinkMode
        IWG::thin_shell_dof_link_mode() const
        {
            return mThinShellDofLinkMode ;
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

        inline const Vector< id_t > &
        IWG::ghost_sideset_ids() const
        {
            return mGhostSideSetIDs ;
        }

//---------------------------------------------------------------------------------

        inline uint
        IWG::number_of_ghost_sidesets() const
        {
            return mGhostSideSetIDs.length();
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

//---------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IWG_HPP
