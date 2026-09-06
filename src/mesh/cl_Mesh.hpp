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

#ifndef BELFEM_CL_MESH_HPP
#define BELFEM_CL_MESH_HPP

#include "typedefs.hpp"
#include "cl_Hash.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Element.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Face.hpp"
#include "cl_ControlPoint.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Mesh_GlobalVariable.hpp"
#include "cl_Mesh_Field.hpp"
#include "cl_ThinShell.hpp"
#include "cl_Curve.hpp"
#include "cl_Bitset.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_TensorMeshConfig.hpp"
#include "cl_Mesh_Periodicity.hpp"
#include "hdf5_types.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace mesh
    {
        class ProtoMesh ;

        class PeriodicityFactory ;

        class GmshReader;

        uint
        compute_facet_index( Facet * aFacet, Element * aElement, Cell< Node * > & aNodes );
    }


    /**
     * @brief Top-level container for all mesh entities.
     *
     * @ingroup grp_mesh
     * @see @ref mesh_mesh_usage_guide
     */
    class Mesh
    {
//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------

        // master proc that owns this mesh, default: 0
        const proc_t mMasterProc;
        const proc_t mCommRank ;
        const proc_t mCommSize ;


        // hash value for unique fingerprint, must be computed first
        // refers to original unaltered mesh without cuts and thin shells
        Hash mHash ;

        //! fingerprint of the SETTINGS this mesh was enriched with ( cuts, thin
        //! shell layers, coating walls ). Opaque here on purpose: the mesh only
        //! carries it into and out of its file, the meaning belongs to whoever
        //! built the mesh. Zero means "not set", which is what every mesh that
        //! was never enriched from an input deck reports
        uint64_t mConfigTag = 0 ;
        string   mConfigText ;

        string mPath ;

        Cell< mesh::Element * > mElements;
        Cell< mesh::Node * >    mNodes;
        Cell< mesh::Facet * >   mFacets;   // facets refer to sidesets
        Cell< mesh::Element * > mVertices;  // not to be confused with mesh vertex
        Cell< mesh::Edge * >    mEdges ;
        Cell< mesh::Face * >    mFaces ;
        Cell< mesh::ControlPoint * > mControlPoints;

        Cell< mesh::Element * > mBoundaryEdges ;

        Cell< mesh::Node * >         mHangingNodes;
        Cell< mesh::Edge * >         mHangingEdges ;
        Cell< mesh::Face * >         mHangingFaces ;
        Cell< mesh::Facet * >        mHangingFacets;
        Cell< mesh::Element * >      mHangingElements ;
        Cell< mesh::ControlPoint * > mHangingControlPoints;

        Cell< mesh::Block * >   mBlocks;
        Cell< mesh::SideSet * > mSideSets;
        Cell< mesh::Curve * >   mCurves ;
        Cell< mesh::ThinShell * > mThinShells ;

        Cell< mesh::GlobalVariable * > mGlobalVariables;
        Cell< mesh::Field * > mFields;

        Map< string, mesh::GlobalVariable * > mGlobalVariableMap;
        Map< id_t, mesh::GlobalVariable * >   mGlobalVariableIDMap;
        Map< string, mesh::Field * > mFieldMap;

        // how many partitions is this mesh split into?
        proc_t mNumberOfPartitions = 1;

        uint mNumberOfDimensions = 0;
        uint mNumberOfGlobalVariables = 0;
        uint mNumberOfFields = 0;

        real mTimeStamp = 0.0;
        uint mTimeStep = 1; // << -- timestep is 1-based for exodus compatibility


        friend class mesh::GmshReader;
        friend class mesh::PeriodicityFactory;
        friend class mesh::ProtoMesh;

        Map< id_t, mesh::Node * >    mNodeMap;
        Map< id_t, mesh::Element * > mElementMap;
        Map< id_t, mesh::Facet * >   mFacetMap;
        Map< id_t, mesh::Block * >   mBlockMap;
        Map< id_t, mesh::SideSet * > mSideSetMap;
        Map< id_t, mesh::Element * > mVertexMap;
        Map< id_t, mesh::Edge * >    mEdgeMap;
        Map< id_t, mesh::Face * >    mFaceMap;
        Map< id_t, mesh::ControlPoint * > mControlPointMap;
        Map< id_t, mesh::Curve * >   mCurveMap ;

        // maximum element order
        uint mMaxElementOrder = 0 ;

        //! flag telling if this is a mesh on the kernel
        bool mIsKernelMesh = false ;

        //! flag telling if mesh is finalized (master only)
        bool mIsFinalized = false ;
        bool mEdgesAreFinalized = false ;
        bool mFacesAreFinalized = false ;

        // we only call the symrcm once per mesh
        bool mNodesAreSorted = false ;

        // compute orientations on finalize
        bool mComputeFacetOrientationsWhenFinalizing = true ;

        //! for special purpose
        Cell< mesh::Node * > mAutoPins ;

        //! for special purpose
        Cell< mesh::Node * > mAbstractNodes ;

        //! for special purpose
        Cell< mesh::Node * > mOrphanedNodes ;

        Bitset< static_cast< size_t>( Connectivity::UNDEFINED )> mConnectivities ;

        const TensorMeshConfig * mTensorConfig = nullptr ;

        mesh::Periodicity * mPeriodicity = nullptr;

        bool mMeshCheckerFlag = false ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * creates an empty mesh container but specifies number of dimension
         */
        Mesh( const uint aNumberOfDimensions, const proc_t aMasterProc=0, const bool aComputeConnectivities = true );

//------------------------------------------------------------------------------

        /**
         * reads a mesh from a file
         */
        Mesh( const string & aPath,
            const proc_t aMasterProc = 0,
            const bool aComputeConnectivities = true,
            const bool aParallelMode = true );

//------------------------------------------------------------------------------

        /**
         * creates a tensor mesh of first, second or third order
         * ( B-splines only for order > 1 )
         */
        Mesh( const uint aOrder,
              const Vector< index_t > aNumNodes,
              const Vector< real > aStep,
              const Vector< real > aOrigin = {},
              const proc_t aMasterProc = 0 );

//------------------------------------------------------------------------------

        /**
         * Multiply all node coordinates with a factor, e.g. to convert from mm to m.
         */
         void
         scale_mesh( const real aFactor );

//------------------------------------------------------------------------------

        ~Mesh();

//------------------------------------------------------------------------------

        const string &
        path() const ;

//------------------------------------------------------------------------------

        bool
        is_tensormesh() const ;

//------------------------------------------------------------------------------

        const TensorMeshConfig *
        tensorconf() const;

//------------------------------------------------------------------------------

        void
        save( const string & aFilePath  );

//------------------------------------------------------------------------------

        uint
        number_of_dimensions() const;

//------------------------------------------------------------------------------

        void
        set_number_of_dimensions( const uint & aNumberOfDimensions );

//------------------------------------------------------------------------------

        index_t
        number_of_nodes() const;

//------------------------------------------------------------------------------

        index_t
        number_of_elements() const;

//------------------------------------------------------------------------------

        index_t
        number_of_facets() const;

//------------------------------------------------------------------------------

        index_t
        number_of_edges() const;

//------------------------------------------------------------------------------

        index_t
        number_of_faces() const;

//------------------------------------------------------------------------------

        index_t
        number_of_control_points() const;

//------------------------------------------------------------------------------

        uint
        number_of_blocks() const;

//------------------------------------------------------------------------------

        uint
        number_of_fields() const;

//------------------------------------------------------------------------------

        uint
        number_of_global_variables() const;

//------------------------------------------------------------------------------

        mesh::Node *
        node( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Node *
        node( const index_t i, const index_t j );

//------------------------------------------------------------------------------

        mesh::Node *
        node( const index_t i, const index_t j, const index_t k );

//------------------------------------------------------------------------------

        mesh::Edge *
        edge( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Face *
        face( const id_t aID );

//------------------------------------------------------------------------------

        bool
        element_exists( const id_t aID ) const;

        mesh::Element *
        element( const id_t aID );

        mesh::Element *
        element( const index_t i, const index_t j );

        mesh::Element *
        element( const index_t i, const index_t j, const index_t k );

//------------------------------------------------------------------------------

        mesh::ControlPoint *
        control_point( const id_t aID );

        mesh::ControlPoint *
        control_point( const index_t i, const index_t j );

        mesh::ControlPoint *
        control_point( const index_t i, const index_t j, const index_t k );

//------------------------------------------------------------------------------

        bool
        facet_exists( const id_t aID ) const;

        mesh::Facet *
        facet( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Block *
        block( const id_t aID );

//------------------------------------------------------------------------------

        mesh::SideSet *
        sideset( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Element *
        vertex( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Basis *
        basis_by_index( const EntityType aType, const index_t aIndex );

        mesh::Basis *
        basis( const EntityType aType, const id_t aID );

//------------------------------------------------------------------------------

        /**
         * get field by label
         */
        mesh::Field *
        field( const string & aLabel );

//------------------------------------------------------------------------------

        mesh::Field *
        field( const index_t aIndex );

//------------------------------------------------------------------------------

        // test if a field has already been created
        bool
        field_exists( const string & aLabel );

//------------------------------------------------------------------------------

        Vector< real > &
        field_data( const string & aLabel );

//------------------------------------------------------------------------------

        Cell< mesh::GlobalVariable * > &
        global_variables() ;

//------------------------------------------------------------------------------

        mesh::GlobalVariable *
        global_variable( const id_t aID );

//------------------------------------------------------------------------------

        mesh::GlobalVariable *
        global_variable( const string & aLabel );

//------------------------------------------------------------------------------

        real &
        global_variable_data( const string & aLabel );

//------------------------------------------------------------------------------

        // test if a global variable has already been created
        bool
        global_variable_exists( const string & aLabel );

//------------------------------------------------------------------------------

        uint
        number_of_sidesets() const;

//------------------------------------------------------------------------------

        Cell< mesh::Node * > & hanging_nodes();
        Cell< mesh::Edge * > & hanging_edges();
        Cell< mesh::Face * > & hanging_faces();
        Cell< mesh::Facet * > & hanging_facets();
        Cell< mesh::Element * > & hanging_elements();
        Cell< mesh::ControlPoint * > & hanging_control_points();

        index_t
        number_of_hanging_basis() const ;

//------------------------------------------------------------------------------

        void
        add_block( mesh::Block * aBlock );

//------------------------------------------------------------------------------

        void
        add_sideset( mesh::SideSet * aSideSet );

//------------------------------------------------------------------------------

        Vector< real > &
        create_field(
                const string & aLabel,
                const EntityType aEntity=EntityType::NODE,
                const id_t aID = 0 );

//------------------------------------------------------------------------------

        real &
        create_global_variable( const string & aLabel,
                                const real aValue = 0.0,
                                const id_t aID = 0 );

//------------------------------------------------------------------------------

        void
        collect_elements_from_blocks();

//------------------------------------------------------------------------------

        void
        collect_facets_from_sidesets();

//------------------------------------------------------------------------------

        void
        update_facet_nodes();

//------------------------------------------------------------------------------

        void
        collect_elements_from_group(
                Cell< mesh::Element * > & aElements,
                const id_t aGroupId,
                const ElementType aType = ElementType::EMPTY );

//------------------------------------------------------------------------------

        void
        compute_facet_orientations();

//------------------------------------------------------------------------------

        /**
         * finalizes all aspects of the mesh; edges and faces are finalized too
         * when connectivities are computed
         */
        void
        finalize();

//------------------------------------------------------------------------------

        bool
        is_finalized() const ;

        bool
        edges_are_finalized() const ;

        bool
        faces_are_finalized() const ;

//------------------------------------------------------------------------------

        void
        collect_hanging_basis();

//------------------------------------------------------------------------------

        void
        expand_hanging_basis_sources();

//------------------------------------------------------------------------------

        /**
         * unfinalizes all aspects of mesh but edges
         */
        void
        unfinalize();

//------------------------------------------------------------------------------

        void
        finalize_edges();

//------------------------------------------------------------------------------

        void
        finalize_edges( Cell< mesh::Element * > & aElements );

//------------------------------------------------------------------------------

        void
        finalize_faces();

//------------------------------------------------------------------------------

        void
        unflag_everything( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_nodes( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_edges( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_faces( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_facets( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_elements( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_vertices( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        void
        unflag_all_control_points( const uint aFlagIndex=0 );

//------------------------------------------------------------------------------

        /**
         * expose Node container
         */
        Cell< mesh::Node * > &
        nodes();

//------------------------------------------------------------------------------

        /**
         * expose Block container
         */
        Cell< mesh::Block * > &
        blocks();

//------------------------------------------------------------------------------

        /**
         * check if block exists
         */
         bool
         block_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------

        /**
         * check if sideset exists
         */
        bool
        sideset_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------

        /**
         * check if node exists
         */
        bool
        node_exists( const id_t aID ) const ;

//------------------------------------------------------------------------------

        /**
         * expose the thinshell container
         */
        Cell< mesh::ThinShell * > &
        thin_shells();

//------------------------------------------------------------------------------

        /**
         * expose Element container
         */
        Cell< mesh::Element * > &
        elements();

//------------------------------------------------------------------------------

        /**
         * expose Facet container
         */
        Cell< mesh::Facet * > &
        facets();

//------------------------------------------------------------------------------

        /**
         * expose Control Point container
         */
        Cell< mesh::ControlPoint * > &
        control_points();

//------------------------------------------------------------------------------

        /**
         * expose Sideset container
         */
        Cell< mesh::SideSet * > &
        sidesets();

//------------------------------------------------------------------------------

        /**
         * expose Curve container
         */
        Cell< mesh::Curve * > &
        curves();

//------------------------------------------------------------------------------

        /**
         * return a curve by its ID
         */
        mesh::Curve *
        curve( const id_t aID );

//------------------------------------------------------------------------------

        void
        create_curve_map();

//------------------------------------------------------------------------------

        /**
         * expose edge container
         */
        Cell< mesh::Edge * > &
        edges();

//------------------------------------------------------------------------------

        /**
         * expose face container
         */
        Cell< mesh::Face * > &
        faces();

//------------------------------------------------------------------------------

        /**
         * expose Vertex container
         */
        Cell< mesh::Element * > &
        vertices();

 //------------------------------------------------------------------------------

         /**
          * expose container for boundary edges
          */
         Cell< mesh::Element * > &
         boundary_edges();

//------------------------------------------------------------------------------

        /**
         * partition the mesh and set element and node ownerships
         */
        void
        partition(
            const uint & aNumberOfPartitions,
            const bool aSetProcOwners = true,
            const bool aForceContinuousPartitions = true,
            const bool aResetVertexContainers = true );

//------------------------------------------------------------------------------

        /**
         * partition the mesh and set element and node ownerships, but also pass selected blocks
         */
        void
        partition( const uint & aNumberOfPartitions,
                   const Vector< id_t > & aSelectedBlocks,
                   const bool aSetProcOwners = true,
                   const bool aForceContinuousPartitions = false,
                   const bool aResetVertexContainers = true );

//------------------------------------------------------------------------------

        /**
         * partition the mesh and set element and node ownerships, but also pass selected blocks
         */
        void
        partition( const uint & aNumberOfPartitions,
                   const Vector< id_t > & aSelectedBlocks,
                   const Vector< id_t > & aSelectedSideSets,
                   const bool aSetProcOwners = true,
                   const bool aForceContinuousPartitions = false,
                   const bool aResetVertexContainers = true );


//------------------------------------------------------------------------------

        void
        update_node_indices();

//------------------------------------------------------------------------------

        void
        update_edge_indices();

//------------------------------------------------------------------------------

        void
        update_face_indices();

//------------------------------------------------------------------------------

        void
        update_facet_indices();

//------------------------------------------------------------------------------

        /**
         * creates the edges on the mesh
         * @param aPrint  : prints the result if flag is set
         * @param aNedelecBlocks   : blocks that carry edge dofs
         * @param aNedelecSideSets : sidesets that carry edge dofs
         * @param aCreateEdgesOnAllSideSets : also create edges on sidesets
         *        that are not listed in aNedelecSideSets
         */
        void
        create_edges( const bool aPrint=false,
                      const Vector< id_t > aNedelecBlocks = Vector< id_t >(),
                      const Vector< id_t > aNedelecSideSets =  Vector< id_t >(),
                      const bool aCreateEdgesOnAllSideSets = true );

//------------------------------------------------------------------------------

        void
        reset_edges();

//------------------------------------------------------------------------------

        /**
         * creates the faces on the mesh
         * @param aPrint  : prints the result if flag is set
         * @param aNedelecBlocks   : blocks that carry face dofs
         * @param aNedelecSideSets : sidesets that carry face dofs
         */
        void
        create_faces(  const bool aPrint=false,
                       const Vector< id_t > aNedelecBlocks = Vector< id_t >(),
                       const Vector< id_t > aNedelecSideSets =  Vector< id_t >() );

        void
        reset_faces();

//------------------------------------------------------------------------------

        void
        update_element_indices();

//------------------------------------------------------------------------------

        void
        update_vertex_indices();

//------------------------------------------------------------------------------

        void
        update_control_point_indices();

//------------------------------------------------------------------------------

        /**
         * return the current timestamp
         */
        real &
        time_stamp();

//------------------------------------------------------------------------------

        /**
         * return the index of the time
         */
        const uint &
        time_step() const;

        uint &
        time_step() ;

//------------------------------------------------------------------------------

        /**
         * set the index of the time step
         */
        void
        set_time_step( const uint aTimeStep );


//------------------------------------------------------------------------------

        void
        set_connectivity( const Connectivity aConnectivity ) ;

//------------------------------------------------------------------------------

        void
        reset_connectivity( const Connectivity aConnectivity ) ;

//------------------------------------------------------------------------------

        bool
        test_connectivity( const Connectivity aConnectivity );

//------------------------------------------------------------------------------

        /**
         * return the ID of the proc that owns this mesh
         */
        const proc_t &
        master() const;

//------------------------------------------------------------------------------

        void
        update_ownerships();

//------------------------------------------------------------------------------

        void
        set_vertex_owners();

//------------------------------------------------------------------------------

        /**
         * to be called by partitioner
         */
        void
        set_number_of_partitions( const proc_t & aNumberOfPartitions );

//------------------------------------------------------------------------------

        /**
         * tells if Nedelec edges exist on this mesh
         * @return
         */
        bool
        edges_exist() const ;

//------------------------------------------------------------------------------

        /**
         * tells if Nedelec faces exist on this mesh
         */
        bool
        faces_exist() const ;

//------------------------------------------------------------------------------

        /**
         * special function called by Kernel
         */
        void
        create_edge_map();

//------------------------------------------------------------------------------

        /**
         * special function called by Kernel
         */
        void
        create_face_map();

//------------------------------------------------------------------------------

        /**
         * return the max interpolation order on the mesh
         */
        const uint &
        max_element_order() const ;

//------------------------------------------------------------------------------

        /**
         * called by kernel
         */
         void
         set_kernel_flag() ;

//------------------------------------------------------------------------------

         bool
         is_kernel_mesh() const ;

//------------------------------------------------------------------------------

        /**
         * flags elements that are curved
         */
         void
         flag_curved_elements();

//------------------------------------------------------------------------------

        void
        distribute_edge_directions();

//------------------------------------------------------------------------------

        void
        save_faces( const string & aPath );

//------------------------------------------------------------------------------

        id_t
        max_node_id() ;

//------------------------------------------------------------------------------

        id_t
        max_element_id() ;

//------------------------------------------------------------------------------

        id_t
        max_block_and_sideset_id() ;

//------------------------------------------------------------------------------

        void
        set_abstract_nodes( Cell< mesh::Node * > & aNodes );

//------------------------------------------------------------------------------

        Cell< mesh::Node * > &
        autopins();

//------------------------------------------------------------------------------

        Cell< mesh::Node * > &
        abstract_nodes();

//------------------------------------------------------------------------------

        Cell< mesh::Node * > &
        orphaned_nodes();

//------------------------------------------------------------------------------

        // compute and return the checksum
        std::size_t
        checksum() ;

//------------------------------------------------------------------------------

        // enforces a checksum instead of computing it; used by BfmFile::load()
        // and by the controller to stamp a derived mesh
        void
        force_checksum( const std::size_t aChecksum ) ;

//------------------------------------------------------------------------------

        // the settings fingerprint that travels into the mesh file, see mConfigTag
        void
        set_config_tag( const uint64_t aTag, const string & aText ) ;

        uint64_t
        config_tag() const ;

        const string &
        config_text() const ;

//------------------------------------------------------------------------------

        void
        update_element_map();

        void
        update_block_map();

        void
        update_sideset_map();

        id_t
        max_node_id() const ;

        id_t
        max_element_id() const ;


        Mesh *
        extract_thin_shell_mesh();

        void
        populate_element_neighbors();

        void
        set_compute_facet_orientation_flag( const bool aFlag );

//------------------------------------------------------------------------------

        /**
         * Compute the total memory used by the mesh in bytes
         */
        size_t
        memory() const ;

//------------------------------------------------------------------------------

        bool
        has_periodicity() const;

        mesh::Periodicity *
        periodicity();

        void
        set_periodicity( mesh::Periodicity * aPeriodicity );

        //void
        //save_fields( const string & aFilename, const uint aRunningTimestep=0 );

        //uint
        //load_fields( const string & aFilename );

        void
        save_meta( hid_t aFile, const uint aRunningTimestep=0 );

        uint
        load_meta( hid_t aFile );

        void
        save_fields( hid_t aFile );

        void
        load_fields( hid_t aFile );

        void
        save_globals( hid_t aFile );

        void
        load_globals( hid_t aFile );

        void
        set_mesh_checker_flag();

        void
        reset_mesh_checker_flag();

        bool
        mesh_checker_flag() const;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        // Resize every existing field whose entity type matches aEntityType
        // to aSize and reset its values to 0.0. Called from finalize() for
        // NODE and ELEMENT, finalize_edges() for EDGE, and finalize_faces()
        // for FACE so that fields registered before mesh-topology updates
        // (cut duplication, partitioning, etc.) stay in sync with the
        // current entity counts.
        void
        resize_and_reset_fields( const EntityType aEntityType, const index_t aSize );

//------------------------------------------------------------------------------

        void
        create_maps();

//------------------------------------------------------------------------------

        void
        reset_maps();

//------------------------------------------------------------------------------

        // set the block ids of each element
        void
        set_block_ids();

//------------------------------------------------------------------------------

        // set the sideset ids of each facet
        void
        set_sideset_ids();



//------------------------------------------------------------------------------

        // set the indices for the sidesets (needed to create curves)
        void
        set_sideset_indices();

//------------------------------------------------------------------------------

        // set the indices for the blocks
        void
        set_block_indices();

//------------------------------------------------------------------------------

        /**
         * compute the max interpolation order on the mesh
         */
        void
        compute_max_element_order();

//------------------------------------------------------------------------------

         void
         compute_edge_directions();

//------------------------------------------------------------------------------

         void
         compute_checksum();

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------

    inline const string &
    Mesh::path() const
    {
        return mPath ;
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::is_tensormesh() const
    {
        return mTensorConfig != nullptr;
    }

//------------------------------------------------------------------------------

    inline const TensorMeshConfig *
    Mesh::tensorconf() const
    {
        return mTensorConfig;
    }

//------------------------------------------------------------------------------

    inline uint
    Mesh::number_of_dimensions() const
    {
        return mNumberOfDimensions;
    }

//------------------------------------------------------------------------------

    inline void
    Mesh::set_number_of_dimensions( const uint & aNumberOfDimensions )
    {
        mNumberOfDimensions = aNumberOfDimensions;
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_nodes() const
    {
        return mNodes.size();
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_elements() const
    {
        return mElements.size();
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_facets() const
    {
        return mFacets.size();
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_edges() const
    {
        return mEdges.size();
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_faces() const
    {
        return mFaces.size();
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_control_points() const
    {
        return mControlPoints.size();
    }

//------------------------------------------------------------------------------

    inline uint
    Mesh::number_of_blocks() const
    {
        return mBlocks.size();
    }

//------------------------------------------------------------------------------

    inline uint
    Mesh::number_of_sidesets() const
    {
        return mSideSets.size();
    }

//------------------------------------------------------------------------------

    inline mesh::Node *
    Mesh::node( const id_t aID )
    {
        BELFEM_ASSERT( mNodeMap.key_exists( aID ),
                      "Tried to access invalid node id: %lu.",
                      ( long unsigned int ) aID );

        return mNodeMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Node *
    Mesh::node( const index_t i, const index_t j )
    {
        BELFEM_ASSERT( this->is_tensormesh() && this->number_of_dimensions() == 2,
            "can only access node over 2d index on 2d-tensormeshes" );

        return mNodes( mTensorConfig->node_index( i, j ) );
    }

//------------------------------------------------------------------------------

    inline mesh::Node *
    Mesh::node( const index_t i, const index_t j, const index_t k )
    {
        BELFEM_ASSERT( this->is_tensormesh() && this->number_of_dimensions() == 3,
            "can only access node over 3d index on 3d-tensormeshes" );

        return mNodes( mTensorConfig->node_index( i, j, k ) );
    }
//------------------------------------------------------------------------------

    inline mesh::Edge *
    Mesh::edge( const id_t aID )
    {
        BELFEM_ASSERT( mEdgeMap.key_exists( aID ),
                      "Tried to access invalid edge id: %lu.",
                      ( long unsigned int ) aID );

        return mEdgeMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Face *
    Mesh::face( const id_t aID )
    {
        BELFEM_ASSERT( mFaceMap.key_exists( aID ),
                      "Tried to access invalid face id: %lu.",
                      ( long unsigned int ) aID );

        return mFaceMap( aID );
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::element_exists( const id_t aID ) const
    {
        return mElementMap.key_exists( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Element *
    Mesh::element( const id_t aID )
    {
        BELFEM_ASSERT( mElementMap.key_exists( aID ),
                      "Tried to access invalid element id: %lu.",
                      ( long unsigned int ) aID );

        return mElementMap( aID ) ;
    }

//------------------------------------------------------------------------------

    inline mesh::Element *
    Mesh::element( const index_t i, const index_t j )
    {
        BELFEM_ASSERT( this->is_tensormesh() && this->number_of_dimensions() == 2,
            "can only access element over 2d index on 2d-tensormeshes" );

        return mElements( mTensorConfig->element_index( i, j ) );
    }

//------------------------------------------------------------------------------

    inline mesh::Element *
    Mesh::element( const index_t i, const index_t j, const index_t k )
    {
        BELFEM_ASSERT( this->is_tensormesh() && this->number_of_dimensions() == 3,
            "can only access element over 3d index on 3d-tensormeshes" );

        return mElements( mTensorConfig->element_index( i, j, k ) );
    }

//------------------------------------------------------------------------------

    inline mesh::ControlPoint *
    Mesh::control_point( const id_t aID )
    {
        BELFEM_ASSERT( mControlPointMap.key_exists( aID ),
                      "Tried to access invalid control point id: %lu.",
                      ( long unsigned int ) aID );

        return mControlPointMap( aID );
    }

    inline mesh::ControlPoint *
    Mesh::control_point( const index_t i, const index_t j )
    {
        BELFEM_ASSERT( this->is_tensormesh() && this->number_of_dimensions() == 2,
            "can only access control point over 2d index on 2d-tensormeshes" );

        return mControlPoints(
            i * mTensorConfig->num_control_points( 1 ) + j  );
    }

    inline mesh::ControlPoint *
    Mesh::control_point( const index_t i, const index_t j, const index_t k )
    {
        BELFEM_ASSERT( this->is_tensormesh() && this->number_of_dimensions() == 3,
            "can only access control point over 3d index on 3d-tensormeshes" );

        return mControlPoints(
            mTensorConfig->num_control_points( 2 ) *
            ( i *  mTensorConfig->num_control_points( 1 ) + j ) + k );
    }


//------------------------------------------------------------------------------

    inline bool
    Mesh::facet_exists( const id_t aID ) const
    {
        return mFacetMap.key_exists( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Facet *
    Mesh::facet( const id_t aID )
    {
        BELFEM_ASSERT( mFacetMap.key_exists( aID ),
                      "Tried to access invalid facet id: %lu.",
                      ( long unsigned int ) aID );

        return mFacetMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Block *
    Mesh::block( const id_t aID )
    {
        BELFEM_ASSERT( mBlockMap.key_exists( aID ),
                      "Tried to access invalid block id: %lu.",
                      ( long unsigned int ) aID );

        return mBlockMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::SideSet *
    Mesh::sideset( const id_t aID )
    {
        BELFEM_ASSERT( mSideSetMap.key_exists( aID ),
                      "Tried to access invalid sideset id: %lu.",
                      ( long unsigned int ) aID );

        return mSideSetMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Element *
    Mesh::vertex( const id_t aID )
    {
        BELFEM_ASSERT( mVertexMap.key_exists( aID ),
                      "Tried to access invalid vertex id: %lu.",
                      ( long unsigned int ) aID );

        return mVertexMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Field *
    Mesh::field( const index_t aIndex )
    {
        return mFields( aIndex );
    }

//------------------------------------------------------------------------------

    inline mesh::Field *
    Mesh::field( const string & aLabel )
    {
        BELFEM_ASSERT( mFieldMap.key_exists( aLabel ),
                      "Field '%s' does not exist on mesh", aLabel.c_str());

        return mFieldMap( aLabel );
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::field_exists( const string & aLabel )
    {
        return mFieldMap.key_exists( aLabel ) ;
    }


//------------------------------------------------------------------------------

    inline Vector< real > &
    Mesh::field_data( const string & aLabel )
    {
        BELFEM_ASSERT( mFieldMap.key_exists( aLabel ),
                      "Field '%s' does not exist on mesh", aLabel.c_str());

        return mFieldMap( aLabel )->data();
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::GlobalVariable * > &
    Mesh::global_variables()
    {
        return mGlobalVariables;
    }

//------------------------------------------------------------------------------

    inline mesh::GlobalVariable *
    Mesh::global_variable( const id_t aID )
    {
        return mGlobalVariableIDMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::GlobalVariable *
    Mesh::global_variable( const string & aLabel )
    {
        return mGlobalVariableMap( aLabel );
    }

//------------------------------------------------------------------------------

    inline real &
    Mesh::global_variable_data( const string & aLabel )
    {
        return mGlobalVariableMap( aLabel )->value();
    }

//------------------------------------------------------------------------------

    // test if a global variable has already been created
    inline bool
    Mesh::global_variable_exists( const string & aLabel )
    {
        return mGlobalVariableMap.key_exists( aLabel );
    }

//------------------------------------------------------------------------------

    inline uint
    Mesh::number_of_fields() const
    {
        return mFields.size() ;
    }

//------------------------------------------------------------------------------

    inline uint
    Mesh::number_of_global_variables() const
    {
        return mGlobalVariables.size();
    }

//------------------------------------------------------------------------------

    /**
     * return the current timestamp
     */
    inline real &
    Mesh::time_stamp()
    {
        return mTimeStamp;
    }

//------------------------------------------------------------------------------

    /**
     * return the current timestep
     */
    inline const uint &
    Mesh::time_step() const
    {
        return mTimeStep;
    }

    inline uint &
    Mesh::time_step()
    {
        return mTimeStep;
    }


//------------------------------------------------------------------------------

    inline void
    Mesh::set_time_step( const uint aTimeStep )
    {
        mTimeStep = aTimeStep;
    }

//------------------------------------------------------------------------------

    inline void
    Mesh::set_connectivity( const Connectivity aConnectivity )
    {
        mConnectivities.set( static_cast<index_t>( aConnectivity ));
    }

//------------------------------------------------------------------------------

    inline void
    Mesh::reset_connectivity( const Connectivity aConnectivity )
    {
        mConnectivities.reset( static_cast<index_t>( aConnectivity ));
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::test_connectivity( const Connectivity aConnectivity )
    {
        return mConnectivities.test( static_cast<index_t>( aConnectivity ));
    }

//------------------------------------------------------------------------------

    inline const proc_t &
    Mesh::master() const
    {
        return mMasterProc;
    }

//------------------------------------------------------------------------------

    inline void
    Mesh::set_number_of_partitions( const proc_t & aNumberOfPartitions )
    {
        mNumberOfPartitions = aNumberOfPartitions;
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::block_exists( const id_t aID ) const
    {
        return mBlockMap.key_exists( aID );
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::sideset_exists( const id_t aID ) const
    {
        return mSideSetMap.key_exists( aID );
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::node_exists( const id_t aID ) const
    {
        return mNodeMap.key_exists( aID );
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::edges_exist() const
    {
        return mEdges.size() > 0 ;
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::faces_exist() const
    {
        return mFaces.size() > 0 ;
    }

//------------------------------------------------------------------------------

    inline const uint &
    Mesh::max_element_order() const
    {
        return mMaxElementOrder;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Node * > &
    Mesh::nodes()
    {
        return mNodes;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Block * > &
    Mesh::blocks()
    {
        return mBlocks;
    }


//------------------------------------------------------------------------------

    inline Cell< mesh::SideSet * > &
    Mesh::sidesets()
    {
        return mSideSets;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Curve * > &
    Mesh::curves()
    {
        return mCurves;
    }

//------------------------------------------------------------------------------

    inline mesh::Curve *
    Mesh::curve( const id_t aID )
    {
        return mCurveMap( aID );
    }


//------------------------------------------------------------------------------

    inline Cell< mesh::ThinShell * > &
    Mesh::thin_shells()
    {
        return mThinShells;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Element * > &
    Mesh::elements()
    {
        return mElements;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Facet * > &
    Mesh::facets()
    {
        return mFacets ;
    }

    inline Cell< mesh::ControlPoint * > &
    Mesh::control_points()
    {
        return mControlPoints;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Edge * > &
    Mesh::edges()
    {
        return mEdges ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Face * > &
    Mesh::faces()
    {
        return mFaces ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Element * > &
    Mesh::boundary_edges()
    {
        return mBoundaryEdges;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Element * > &
    Mesh::vertices()
    {
        return mVertices;
    }

//------------------------------------------------------------------------------

    inline void
    Mesh::set_kernel_flag()
    {
        mIsKernelMesh = true ;
    }

//------------------------------------------------------------------------------

    inline bool
    Mesh::is_kernel_mesh() const
    {
        return mIsKernelMesh ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Node * > &
    Mesh::autopins()
    {
        return mAutoPins ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Node * > &
    Mesh::abstract_nodes()
    {
        return mAbstractNodes ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Node * > &
    Mesh::orphaned_nodes()
    {
        return mOrphanedNodes ;
    }

//------------------------------------------------------------------------------

    inline mesh::Basis *
    Mesh::basis_by_index( const EntityType aType, const index_t aIndex )
    {
        switch ( aType )
        {
            case EntityType::NODE :
            {
                return mNodes( aIndex );
            }
            case EntityType::EDGE :
            {
                return mEdges( aIndex );
            }
            case EntityType::FACE :
            {
                return mFaces( aIndex );
            }
            case EntityType::ELEMENT :
            {
                return mElements( aIndex );
            }
            case EntityType::FACET :
            {
                return mFacets( aIndex );
            }
            case EntityType::CONTROLPOINT :
            {
                return mControlPoints( aIndex );
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid entity type" );
                return nullptr;
            }
        }

    }

    inline mesh::Basis* Mesh::basis( const EntityType aType, const id_t aID )
    {
        switch ( aType )
        {
            case EntityType::NODE :
            {
                return this->node( aID );
            }
            case EntityType::EDGE :
            {
                return this->edge( aID );
            }
            case EntityType::FACE :
            {
                return this->face( aID );
            }
            case EntityType::ELEMENT :
            {
                return this->element( aID );
            }
            case EntityType::FACET :
            {
                return this->facet( aID );
            }
            case EntityType::CONTROLPOINT :
            {
                return this->control_point( aID );
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid entity type" );
                return nullptr;
            }
        }
    }

    inline void
    Mesh::set_compute_facet_orientation_flag( const bool aFlag )
    {
        mComputeFacetOrientationsWhenFinalizing = aFlag ;
    }

    inline bool
    Mesh::is_finalized() const
    {
        return mIsFinalized ;
    }

    inline bool
    Mesh::edges_are_finalized() const
    {
        return mEdgesAreFinalized ;
    }


    inline bool
    Mesh::faces_are_finalized() const
    {
        return mFacesAreFinalized ;
    }


//------------------------------------------------------------------------------

    inline Cell< mesh::Node * > & Mesh::hanging_nodes() { return mHangingNodes; }
    inline Cell< mesh::Edge * > & Mesh::hanging_edges() { return mHangingEdges; }
    inline Cell< mesh::Face * > & Mesh::hanging_faces() { return mHangingFaces; }
    inline  Cell< mesh::Facet * > & Mesh::hanging_facets() { return mHangingFacets; }
    inline Cell< mesh::Element * > & Mesh::hanging_elements() { return mHangingElements; }
    inline Cell< mesh::ControlPoint * > & Mesh::hanging_control_points() { return mHangingControlPoints; }

    inline index_t
    Mesh::number_of_hanging_basis() const
    {
        return   mHangingNodes.size()
               + mHangingEdges.size()
               + mHangingFaces.size()
               + mHangingFacets.size()
               + mHangingElements.size()
               + mHangingControlPoints.size();
    }


    inline void
    Mesh::set_periodicity( mesh::Periodicity * aPeriodicity )
    {
        BELFEM_ASSERT( mPeriodicity == nullptr, "Periodicity object has already been set for this mesh" );
        mPeriodicity = aPeriodicity;
    }

    inline bool
    Mesh::has_periodicity() const
    {
        return mPeriodicity != nullptr;
    }

    inline mesh::Periodicity *
    Mesh::periodicity()
    {
        BELFEM_ASSERT( mPeriodicity != nullptr, "Periodicity object has not been set for this mesh" );

        return mPeriodicity;
    }


    inline void
    Mesh::set_mesh_checker_flag()
    {
        mMeshCheckerFlag = true;
    }

    inline void
    Mesh::reset_mesh_checker_flag()
    {
        mMeshCheckerFlag = false;
    }

    inline bool
    Mesh::mesh_checker_flag() const
    {
        return mMeshCheckerFlag;
    }


   namespace mesh
   {
        template< typename T >
        void
        collect_hanging_basis( Cell< T * > & aBasis, Cell< T * > & aHangingBasis )
        {
            index_t tCount = 0 ;
            for( T * tBasis : aBasis )
            {
                if( tBasis->is_hanging() )
                {
                    ++tCount ;
                }
            }

            if( tCount == 0 )
            {
                aHangingBasis.clear();
            }
            else
            {
                aHangingBasis.set_size( tCount, nullptr );
                tCount = 0;
                for ( T * tBasis: aBasis )
                {
                    if ( tBasis->is_hanging())
                    {
                        aHangingBasis( tCount++ ) = tBasis;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
   }
}

#endif //BELFEM_CL_MESH_HPP
