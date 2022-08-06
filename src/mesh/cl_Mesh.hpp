//
// Created by Christian Messe on 2019-07-25.
//


#ifndef BELFEM_CL_MESH_HPP
#define BELFEM_CL_MESH_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Element.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Face.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Mesh_GlobalVariable.hpp"
#include "cl_Mesh_Field.hpp"
#include "st_Mesh_Layer.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    namespace mesh
    {
        class GmshReader;
    }

    class Mesh
    {
//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------

        // master proc that owns this mesh, default: 0
        const proc_t mMasterProc;


        Cell< mesh::Element * > mElements;

        Cell< mesh::Facet * >   mConnectors;

        Cell< mesh::Node * >    mNodes;
        Cell< mesh::Facet * >   mFacets;   // facets refer to sidesets
        Cell< mesh::Element * > mVertices;  // not to be confused with mesh vertex
        Cell< mesh::Element * > mBoundaryEdges ;
        Cell< mesh::Edge * >    mEdges ;
        Cell< mesh::Face * >    mFaces ;

        Cell< mesh::Block * >   mBlocks;
        Cell< mesh::SideSet * > mCuts; // special for Maxwell cuts
        Cell< mesh::SideSet * > mSideSets;

        Cell< mesh::GlobalVariable * > mGlobalVariables;
        Cell< mesh::Field * > mFields;

        Map< string, mesh::GlobalVariable * > mGlobalVariableMap;
        Map< string, mesh::Field * > mFieldMap;

        // how many partitions is this mesh split into?
        proc_t mNumberOfPartitions = 1;

        uint mNumberOfDimensions = 0;
        uint mNumberOfGlobalVariables = 0;
        uint mNumberOfFields = 0;

        real mTimeStamp = 0.0;
        uint mTimeStep = 1; // << -- timestep is 1-based for exodus compatibility

        // flag that tells if connectivities other than element to node are to be computed
        // this flag must be on if a FEM model is supposed to be run on the mesh
        const bool mComputeConnectivities;

        // flag telling if facets have already been linked
        bool mFacetsAreLinked = false;

        friend mesh::GmshReader;

        Map< id_t, mesh::Node * >    mNodeMap;
        Map< id_t, mesh::Element * > mElementMap;
        Map< id_t, mesh::Facet * >   mFacetMap;
        Map< id_t, mesh::Block * >   mBlockMap;
        Map< id_t, mesh::SideSet * > mSideSetMap;
        Map< id_t, mesh::SideSet * > mCutMap;
        Map< id_t, mesh::Element * > mVertexMap;
        Map< id_t, mesh::Edge * >    mEdgeMap;
        Map< id_t, mesh::Face * >    mFaceMap;

        // contains list of original ids ( row 0 ) vs duplicates (row 1)
        // written by scissors if cut is performed
        Matrix< id_t > mNodeCutTable ;
        Map< id_t, mesh::Node * > mNodeCutMap ;

        // contains a list of cloned facets vs originals (for tapes)
        Matrix< id_t > mTapeFacetTable ;
        Map< id_t, mesh::Facet * > mTapeFacetMap ;

        // maximum element order
        uint mMaxElementOrder = 0 ;

        //! flag telling if this is a mesh on the kernel
        bool mIsKernelMesh = false ;

        //! flag telling if mesh is finalized (master only)
        bool mIsFinalized = false ;

        //! list of blocks with Nedelec Elements
        Vector< id_t > mNedelecBlocks ;

        //! list of sidesets with Nedelec Elements
        Vector< id_t > mNedelecSideSets ;

        //! id of ghost sidesets
        Vector< id_t > mGhostSideSetIDs ;

        //! ids of original facets that have ghosts
        Vector< id_t > mGhostFacetIDs ;

        //! the map links the element ids to the index in ghost elements
        Map< id_t, index_t > mGhostFacetMap ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * creates an empty mesh container
         */
        Mesh( const proc_t aMasterProc = 0, const bool aComputeConnectivities = true );

//------------------------------------------------------------------------------

        /**
         * creates an empty mesh container but specifies number of dimension
         */
        Mesh( const uint & aNumberOfDimensions, const proc_t aMasterProc, const bool aComputeConnectivities = true );

//------------------------------------------------------------------------------

        /**
         * reads a mesh from a file
         */
        Mesh( const string & aPath, const proc_t aMasterProc = 0, const bool aComputeConnectivities = true );

//------------------------------------------------------------------------------

        /**
         * Multiply all node coordinates with a factor, e.g. to convert from mm to m.
         */
         void
         scale_mesh( const real aFactor );

//------------------------------------------------------------------------------

        ~Mesh();

//------------------------------------------------------------------------------

        void
        save( const string & aFilePath );

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
        number_of_connectors() const;

//------------------------------------------------------------------------------

        index_t
        number_of_edges() const;

//------------------------------------------------------------------------------

        index_t
        number_of_faces() const;

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

        mesh::Edge *
        edge( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Face *
        face( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Element *
        element( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Facet *
        facet( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Block *
        block( const id_t aID );

//------------------------------------------------------------------------------

        mesh::SideSet *
        sideset( const id_t aID );

//------------------------------------------------------------------------------

        mesh::SideSet *
        cut( const id_t aID );

//------------------------------------------------------------------------------

        mesh::Element *
        vertex( const id_t aID );

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

        mesh::GlobalVariable *
        global_variable( const index_t aIndex );

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

        /*mesh::SideSet *
        sideset( const index_t aIndex );*/

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
        collect_facets_from_cuts();

//------------------------------------------------------------------------------

        void
        collect_elements_from_group(
                Cell< mesh::Element * > & aElements,
                const id_t aGroupId,
                const ElementType aType = ElementType::EMPTY );

//------------------------------------------------------------------------------

        void
        collect_edges_from_group(
                Cell< mesh::Element * > & aElements,
                const id_t aGroupId );

//------------------------------------------------------------------------------

        void
        compute_facet_orientations();

//------------------------------------------------------------------------------

        /**
         * finalizes all aspects of mesh but edges
         */
        void
        finalize();

//------------------------------------------------------------------------------

        /**
         * unfinalizes all aspects of mesh but edges
         */
        void
        unfinalize();

//------------------------------------------------------------------------------

        void
        finalize_edges( Cell< mesh::Element * > & aElements );

//------------------------------------------------------------------------------

        void
        finalize_faces();

//------------------------------------------------------------------------------

        void
        unflag_everything();

//------------------------------------------------------------------------------

        void
        unflag_all_nodes();

//------------------------------------------------------------------------------

        void
        unflag_all_edges();

//------------------------------------------------------------------------------

        void
        unflag_all_faces();

//------------------------------------------------------------------------------

        void
        unflag_all_facets();

//------------------------------------------------------------------------------

        void
        unflag_all_connectors();

//------------------------------------------------------------------------------

        void
        unflag_all_elements();

//------------------------------------------------------------------------------

        void
        unflag_all_vertices();

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
         * expose Cut container
         */
        Cell< mesh::SideSet * > &
        cuts();

//------------------------------------------------------------------------------

        /**
         * returns the cut duplicate of a node
         */
        mesh::Node *
        duplicate_node( const id_t aID );
//------------------------------------------------------------------------------

        /**
         * returns the original facet the ghost was cloned from
         */

        mesh::Facet *
        original_facet( const id_t aID ) ;

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
         * check if cut exists
         */
        bool
        cut_exists( const id_t aID ) const ;

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
         * expose Connector container
         */
        Cell< mesh::Facet * > &
        connectors();

//------------------------------------------------------------------------------

        /**
         * expose Sideset container
         */
        Cell< mesh::SideSet * > &
        sidesets();

//------------------------------------------------------------------------------

        /**
         * expose edge container
         */
        Cell< mesh::Edge * > &
        edges();

//------------------------------------------------------------------------------

        /**
         * expose boundary edge container
         */
        Cell< mesh::Element * > &
        boundary_edges();

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
         * returns the IDs of ghost sidesets needed for thin shells
         */
         const Vector< id_t > &
         ghost_sideset_ids() const ;

//------------------------------------------------------------------------------

        /**
         * partition the mesh and set element and node ownerships
         */
        void
        partition( const uint & aNumberOfPartitions, const bool aSetProcOwners = true );

//------------------------------------------------------------------------------

        /**
         * partition the mesh and set element and node ownerships, but also pass selected blocks
         */
        void
        partition( const uint & aNumberOfPartitions,
                   const Vector< id_t > & aSelectedBlocks,
                   const bool aSetProcOwners = true );

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
         */
        void
        create_edges( const bool aPrint=false,
                      const Vector< id_t > aNedelecBlocks = Vector< id_t >(),
                      const Vector< id_t > aNedelecSideSets =  Vector< id_t >() );

//------------------------------------------------------------------------------

        void
        reset_edges();

//------------------------------------------------------------------------------

        /**
         * creates the edges on the mesh
         * @param aPrint  : prints the result if flag is set
         */
        void
        create_faces(  const bool aPrint=false,
                       const Vector< id_t > aNedelecBlocks = Vector< id_t >(),
                       const Vector< id_t > aNedelecSideSets =  Vector< id_t >()  );

        void
        reset_faces();

//------------------------------------------------------------------------------

        void
        update_element_indices();

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
         * return the index of the time
         */
        void
        set_time_step( const uint aTimeStep );


//------------------------------------------------------------------------------

        /**
         * a reader may call this function to tell the mesh that the faces
         * have been linked already, no need to do this in finalize() again.
         */
        void
        set_facets_are_linked_flag();

//------------------------------------------------------------------------------

        /**
         * return the ID of the proc that owns this mesh
         */
        const proc_t &
        master() const;

//------------------------------------------------------------------------------

        void
        set_node_owners();

//------------------------------------------------------------------------------

        void
        set_vertex_owners();

//------------------------------------------------------------------------------

        void
        set_connector_owners();


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
         * @return
         */
        bool
        faces_exist() const ;


//------------------------------------------------------------------------------

        /**
         * special function called by Kernel
         * @return
         */
        void
        create_edge_map();

//------------------------------------------------------------------------------

        /**
         * special function called by Kernel
         * @return
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
         * expose the cut table for the nodes
         */
        Matrix< id_t > &
        node_cut_table() ;

//------------------------------------------------------------------------------

        /**
         * expose the facet table for the thin shells
         */
        Matrix< id_t > &
        tape_facet_table() ;

//------------------------------------------------------------------------------

        /**
         * return the list with blocks that contain edge elements
         */
         const Vector< id_t > &
         nedelec_blocks() const ;

//------------------------------------------------------------------------------

        /**
         * return the list with sidesets that contain edge elements
         */
        const Vector< id_t > &
        nedelec_sidesets() const ;

//------------------------------------------------------------------------------

        /**
         * flags elements that are curved
         */
         void
         flag_curved_elements();

//------------------------------------------------------------------------------

        /**
         * create the ghost layers for selected elements
         * called by tape roller on master proc.
         */
         void
         create_ghost_sidesets(
                 const Vector< id_t >    & aGhostSideSetIDs,
                 const Vector< id_t >    & aElementIDs,
                 Cell< mesh::Layer  * >  & aLayers );

//------------------------------------------------------------------------------

        /**
         * create the ghost layers for selected elements
         * called by non-master procs from kernel
         */
         void
         distribute_ghost_sidesets( const proc_t aTarget, const proc_t aMasterProc=0 );


//------------------------------------------------------------------------------

        /**
         * get the number of layers if this is thin shell
         */
         index_t
         number_of_thin_shell_layers() const ;

//------------------------------------------------------------------------------

        /**
         * needed to create thin shell element
         */
        mesh::Facet *
        ghost_facet( const id_t aID, const uint aLayer );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        connect_nodes_to_elements();

//------------------------------------------------------------------------------

        void
        connect_edges_to_elements( Cell< mesh::Element * > & aElements );

//------------------------------------------------------------------------------

        void
        connect_nodes_to_facets();

//------------------------------------------------------------------------------

        void
        connect_facets_to_elements();

//------------------------------------------------------------------------------

        void
        connect_elements_to_elements();

//------------------------------------------------------------------------------

        void
        connect_nodes_to_nodes();

//------------------------------------------------------------------------------

        void
        connect_edges_to_edges();

//------------------------------------------------------------------------------

        /**
         * connects both nodes and edges in vertex container
         */
        void
        connect_vertices_to_vertices();

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

        /**
         * compute the max interpolation order on the mesh
         */
        void
        compute_max_element_order();

//------------------------------------------------------------------------------
    };

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
    Mesh::number_of_connectors() const
    {
        return mConnectors.size();
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

    inline mesh::Element *
    Mesh::element( const id_t aID )
    {
        BELFEM_ASSERT( mElementMap.key_exists( aID ),
                      "Tried to access invalid element id: %lu.",
                      ( long unsigned int ) aID );

        return mElementMap( aID ) ;
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

    inline mesh::SideSet *
    Mesh::cut( const id_t aID )
    {
        BELFEM_ASSERT( mCutMap.key_exists( aID ),
                      "Tried to access invalid cut id: %lu.",
                      ( long unsigned int ) aID );

        return mCutMap( aID );
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

    inline mesh::GlobalVariable *
    Mesh::global_variable( const index_t aIndex )
    {
        return mGlobalVariables( aIndex );
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
    Mesh::set_facets_are_linked_flag()
    {
        mFacetsAreLinked = true;
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
    Mesh::cut_exists( const id_t aID ) const
    {
        return mCutMap.key_exists( aID );
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
    Mesh::cuts()
    {
        return mCuts;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::SideSet * > &
    Mesh::sidesets()
    {
        return mSideSets;
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

//------------------------------------------------------------------------------

    inline Cell< mesh::Facet * > &
    Mesh::connectors()
    {
        return mConnectors ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Edge * > &
    Mesh::edges()
    {
        return mEdges ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Element * > &
    Mesh::boundary_edges()
    {
        return mBoundaryEdges ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Face * > &
    Mesh::faces()
    {
        return mFaces ;
    }

//------------------------------------------------------------------------------

    inline Cell< mesh::Element * > &
    Mesh::vertices()
    {
        return mVertices;
    }

//------------------------------------------------------------------------------

    inline const Vector< id_t > &
    Mesh::ghost_sideset_ids() const
    {
        return mGhostSideSetIDs ;
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

    inline Matrix< id_t > &
    Mesh::node_cut_table()
    {
        return mNodeCutTable ;
    }

//------------------------------------------------------------------------------

    inline Matrix< id_t > &
    Mesh::tape_facet_table()
    {
        return mTapeFacetTable ;
    }

//------------------------------------------------------------------------------

    inline mesh::Node *
    Mesh::duplicate_node( const id_t aID )
    {
        return mNodeCutMap.key_exists( aID )
            ? mNodeCutMap( aID ) : mNodeMap( aID );
    }

//------------------------------------------------------------------------------

    inline mesh::Facet *
    Mesh::original_facet( const id_t aID )
    {
        return mTapeFacetMap.key_exists( aID )
            ? mTapeFacetMap( aID ) : mFacetMap( aID );
    }

//------------------------------------------------------------------------------

    inline const Vector< id_t > &
    Mesh::nedelec_blocks() const
    {
        return mNedelecBlocks ;
    }

//------------------------------------------------------------------------------

    inline const Vector< id_t > &
    Mesh::nedelec_sidesets() const
    {
        return mNedelecSideSets ;
    }

//------------------------------------------------------------------------------

    inline index_t
    Mesh::number_of_thin_shell_layers() const
    {
        return mGhostSideSetIDs.length() ;
    }

//------------------------------------------------------------------------------

    inline mesh::Facet *
    Mesh::ghost_facet( const id_t aID, const uint aLayer )
    {
        return mSideSetMap( mGhostSideSetIDs( aLayer ) )->facet_by_index( mGhostFacetMap( aID ) );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_CL_MESH_HPP
