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

#ifndef BELFEM_CL_MESH_DISTRIBUTOR_HPP
#define BELFEM_CL_MESH_DISTRIBUTOR_HPP

#include "typedefs.hpp"

#include "cl_DynamicBitset.hpp"
#include "cl_Cell.hpp"

#include "cl_Mesh.hpp"
#include "cl_CommTable.hpp"

#include "cl_ProtoMesh.hpp"

namespace belfem
{
    namespace mesh
    {

        class Distributor
        {
            const proc_t mCommRank ;
            const proc_t mCommSize ;

            Mesh * mMesh = nullptr ;

            bool mOwnMesh = false ;

            DynamicBitset * mNodeBitset     = nullptr ;
            DynamicBitset * mElementBitset  = nullptr ;
            DynamicBitset * mEdgeBitset     = nullptr ;
            DynamicBitset * mFaceBitset     = nullptr ;
            DynamicBitset * mFacetBitset    = nullptr ;
            DynamicBitset * mVertexBitset   = nullptr ;
            DynamicBitset * mControlPointBitset = nullptr ;

            // work Cell with indices
            index_t  mNumberOfAllEntities[ 14 ];

            index_t & mNumberOfDimensions       = mNumberOfAllEntities[ 0 ];
            index_t & mMaxElementOrder          = mNumberOfAllEntities[ 1 ];
            index_t & mNumberOfAllNodes         = mNumberOfAllEntities[ 2 ];
            index_t & mNumberOfAbstractNodes    = mNumberOfAllEntities[ 3 ];
            index_t & mNumberOfAllEdges         = mNumberOfAllEntities[ 4 ];
            index_t & mNumberOfAllFaces         = mNumberOfAllEntities[ 5 ];
            index_t & mNumberOfAllElements      = mNumberOfAllEntities[ 6 ];
            index_t & mNumberOfAllFacets        = mNumberOfAllEntities[ 7 ];
            index_t & mNumberOfAllVertices      = mNumberOfAllEntities[ 8 ];
            index_t & mNumberOfAllControlPoints = mNumberOfAllEntities[ 9 ];
            index_t & mNumberOfAllBlocks        = mNumberOfAllEntities[ 10 ];
            index_t & mNumberOfAllSideSets      = mNumberOfAllEntities[ 11 ];
            index_t & mNumberOfAllThinShells    = mNumberOfAllEntities[ 12 ];
            index_t & mNumberOfAllTMatrices     = mNumberOfAllEntities[ 13 ];

            Cell< index_t > mIndices ;

            Cell< proto::NodeData *  >        mNodeData ;
            Cell< proto::ElementData * >      mElementData ;
            Cell< proto::ElementExtra * >     mElementExtra ;
            Cell< proto::EdgeData *  >        mEdgeData ;
            Cell< proto::FaceData *  >        mFaceData ;
            Cell< proto::FacetData * >        mFacetData ;
            Cell< proto::FacetExtra * >       mFacetExtra ;
            Cell< proto::ElementData * >      mVertexData ;
            Cell< proto::ControlPointData * > mControlPointData ;
            Cell< proto::TMatrixData * >      mTMatrices ;

            Cell< CommTable * >        mTables ;

            ProtoMesh * mProtoMesh = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Distributor( Mesh * aMesh = nullptr ) ;

            ~Distributor() ;

            void
            run();

            Cell< CommTable * > &
            tables();

            Mesh *
            partial_mesh();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            count_entities();

            void
            select_entities( const proc_t aTarget );

            void
            select_sources( mesh::Basis * aBasis );

            void
            populate_t_matrices( const proc_t aTarget );

            void
            populate_node_data( const proc_t aTarget );

            void
            populate_element_data( const proc_t aTarget );

            void
            populate_element_extra( const proc_t aTarget );

            void
            populate_edge_data( const proc_t aTarget );

            void
            populate_face_data( const proc_t aTarget );

            void
            populate_facet_data( const proc_t aTarget );

            void
            populate_facet_extra( const proc_t aTarget );

            void
            populate_vertex_data( const proc_t aTarget );

            void
            populate_control_point_data( const proc_t aTarget );



//------------------------------------------------------------------------------

            void
            send_thinshell_data();

            void
            receive_thinshell_data();

//------------------------------------------------------------------------------

            void
            send_block_data();

            void
            receive_block_data();

            void
            send_sideset_data();

            void
            receive_sideset_data();

//------------------------------------------------------------------------------

            void
            send_node_data();

            void
            receive_node_data();

//------------------------------------------------------------------------------

            void
            send_element_data();

            void
            receive_element_data();

//------------------------------------------------------------------------------

            void
            send_edge_data();

            void
            receive_edge_data();

//------------------------------------------------------------------------------

            void
            send_face_data();

            void
            receive_face_data();

//------------------------------------------------------------------------------

            void
            send_facet_data();

            void
            receive_facet_data();

//------------------------------------------------------------------------------

            void
            send_element_extra( const proc_t aTarget = 0);

            void
            receive_element_extra();

//------------------------------------------------------------------------------

            void
            send_facet_extra();

            void
            receive_facet_extra();

//------------------------------------------------------------------------------

            void
            send_vertex_data();

            void
            receive_vertex_data();

//------------------------------------------------------------------------------

            void
            send_control_point_data();

            void
            receive_control_point_data();

//------------------------------------------------------------------------------

            void
            send_t_matrices();

            void
            receive_t_matrices();

//------------------------------------------------------------------------------

            void
            flag_curved_elements( const proc_t aTarget );

            void
            flag_curved_facets( const proc_t aTarget );

//------------------------------------------------------------------------------

            void
            create_bitsets();

            void
            delete_bitsets();

            void
            reset_bitsets();

//------------------------------------------------------------------------------

            void
            delete_node_data();

            void
            delete_element_data();

            void
            delete_element_extra();

            void
            delete_edge_data();

            void
            delete_face_data();

            void
            delete_facet_data();

            void
            delete_facet_extra();

            void
            delete_vertex_data();

            void
            delete_tables();

            void
            delete_mesh_data();

            void
            delete_t_matrices();

//------------------------------------------------------------------------------
        };

        inline
        Cell< CommTable * > &
        Distributor::tables()
        {
            return mTables;
        }


    }
}
#endif // BELFEM_CL_MESH_DISTRIBUTOR_HPP
