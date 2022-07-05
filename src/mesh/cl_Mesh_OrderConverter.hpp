//
// Created by Christian Messe on 22.10.19.
//

#ifndef BELFEM_CL_MESH_ORDERCONVERTER_HPP
#define BELFEM_CL_MESH_ORDERCONVERTER_HPP

#include "cl_Mesh.hpp"
#include "cl_Matrix.hpp"
namespace belfem
{
    // can change the order of a mesh, eg create
    // second order form a linear mesh
    class OrderConverter
    {
        Mesh * mInMesh;
        Mesh * mMesh;

        // 0: Master index
        // 1: Slave index
        // 2: Facet index on master
        // 3: Facet index on slave
        Matrix< index_t > mSharedFacetTable;

        // 0: Master index
        // 1: Facet index on master
        Matrix< index_t > mUniqueFacetTable;

        // contains indices of first and second node
        Matrix< index_t > mEdges;

        // contains node ID per edge
        Map< index_t, index_t > mEdgeMap;

        Cell< mesh::Node * > mOriginalNodes;
        Cell< mesh::Node * > mEdgeNodes;
        Cell< mesh::Node * > mFacetNodes;
        Cell< mesh::Node * > mCenterNodes;

        Matrix< index_t > mFacesPerElement;

        id_t mNodeID;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        OrderConverter( Mesh * aInMesh, Mesh * aOutMesh = nullptr );

//------------------------------------------------------------------------------

        ~OrderConverter() = default;

//------------------------------------------------------------------------------

        Mesh *
        mesh();

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        bool
        check_input_mesh();

//------------------------------------------------------------------------------

        void
        create_shared_facets();

//------------------------------------------------------------------------------

        void
        create_unique_facets();

//------------------------------------------------------------------------------

        uint
        count_common_nodes( mesh::Element * aElementA, mesh::Element * aElementB );

//------------------------------------------------------------------------------

        void
        get_facet_nodes(
                mesh::Element * aElementA,
                mesh::Element * aElementB,
                Cell< mesh::Node * > & aNodes );

//------------------------------------------------------------------------------

        uint
        get_facet_index(  mesh::Element * aElement, Cell< mesh::Node * > & aNodes );

//------------------------------------------------------------------------------

        void
        create_edges();

//------------------------------------------------------------------------------

        void
        get_max_node_id();

//------------------------------------------------------------------------------

        void
        copy_nodes();

//------------------------------------------------------------------------------

        void
        create_edge_nodes();

//------------------------------------------------------------------------------

        void
        create_facet_nodes_2d();

//------------------------------------------------------------------------------

        void
        create_facet_nodes_3d();

//------------------------------------------------------------------------------

        void
        create_center_nodes_2d();

//------------------------------------------------------------------------------

        void
        create_center_nodes_3d();

//------------------------------------------------------------------------------

        void
        create_nodes();

//------------------------------------------------------------------------------

        void
        create_elements();

//------------------------------------------------------------------------------

        void
        link_edge_nodes();

//------------------------------------------------------------------------------

        void
        link_facet_nodes();

//------------------------------------------------------------------------------

        void
        link_center_nodes_2d();

//------------------------------------------------------------------------------

        void
        link_center_nodes_3d();

//------------------------------------------------------------------------------

        void
        facet_table( const ElementType & aType, Vector< uint > & aTable );

//------------------------------------------------------------------------------

        void
        create_blocks();

//------------------------------------------------------------------------------

        void
        create_sidesets();

//------------------------------------------------------------------------------

        ElementType
        upgrade_type( const ElementType & aType );

//------------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_MESH_ORDERCONVERTER_HPP
