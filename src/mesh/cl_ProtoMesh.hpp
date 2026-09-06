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

#ifndef BELFEM_CL_MESHDATA_HPP
#define BELFEM_CL_MESHDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Map.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Face.hpp"
#include "cl_Element.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_ThinShell.hpp"

#include "st_ProtoMesh.hpp"

namespace belfem
{
    class Mesh ;

    namespace mesh
    {
//------------------------------------------------------------------------------

        class ProtoMesh
        {
            Mesh * mMesh        = nullptr ;

            proto::MetaData     * mMetaData     = nullptr ;
            proto::NodeData     * mNodeData     = nullptr ;
            proto::ElementData  * mElementData  = nullptr ;

            proto::EdgeData     * mEdgeData    = nullptr ;
            proto::FaceData     * mFaceData    = nullptr ;
            proto::FacetData    * mFacetData   = nullptr ;

            proto::ElementExtra * mElementExtra = nullptr ;
            proto::FacetExtra   * mFacetExtra = nullptr ;

            proto::ElementData      * mVertexData  = nullptr ;
            proto::ControlPointData * mControlPointData = nullptr ;
            proto::PeriodicityData  * mPeriodicityData = nullptr ;
            proto::TMatrixData      * mTMatrixData  = nullptr ;

            Cell< proto::GroupData >     mBlockData ;
            Cell< proto::GroupData >     mSideSetData ;
            Cell< proto::ThinShellData > mThinShellData ;

            Map< id_t, Node * >         mNodeMap ;
            Map< id_t, Element * >      mElementMap ;
            Map< id_t, Edge * >         mEdgeMap ;

            Map< id_t, Face * >         mFaceMap ;
            Map< id_t, Facet * >        mFacetMap ;
            Map< id_t, Element * >      mVertexMap ;
            Map< id_t, Block * >        mBlockMap ;
            Map< id_t, SideSet * >      mSideSetMap ;
            Map< id_t, ThinShell * >    mThinShellMap ;
            Map< id_t, ControlPoint * > mControlPointMap ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ProtoMesh( Mesh * aMesh );

            ~ProtoMesh();

            Cell< proto::GroupData > &
            block_data();

            proto::GroupData &
            block_data( const index_t aIndex );

            Cell< proto::GroupData > &
            sideset_data();

            proto::GroupData &
            sideset_data( const index_t aIndex );

            Cell< proto::ThinShellData > &
            thin_shell_data();

            proto::ThinShellData &
            thin_shell_data( const index_t aIndex );

            proto::MetaData *
            meta_data();

            proto::NodeData *
            node_data();

            proto::ElementData *
            element_data();

            proto::EdgeData *
            edge_data();

            proto::FaceData *
            face_data();

            proto::FacetData *
            facet_data();

            proto::ElementData *
            vertex_data();

            proto::ElementExtra *
            element_extra();

            proto::FacetExtra *
            facet_extra();

            proto::ControlPointData *
            control_point_data();

            proto::PeriodicityData *
            periodicity_data();

            proto::TMatrixData *
            t_matrix_data();

//------------------------------------------------------------------------------

            void
            create_nodes();

            void
            create_elements();

            void
            create_edges();

            void
            create_faces();

            void
            create_facets();

            void
            create_vertices();

            void
            create_control_points();

            void
            create_element_extra();

            void
            create_facet_extra();

            void
            edge();

            // aKeepEmpty: the distributor drops sidesets with no facets on this proc,
            // but the BFM loader must keep them ( input domain types are re-applied
            // by id on reload, and thin-shell source sidesets are empty by design )
            void
            create_sidesets( const bool aKeepEmpty = false );

            void
            create_thinshells();

            void
            create_periodicitiy( const bool aCrosslinkEntities = false );

            void
            create_t_matrices();

//------------------------------------------------------------------------------

            void
            populate_meta_data();

            void
            populate_node_data();

            void
            populate_block_data();

            void
            populate_sideset_data();

            void
            populate_element_data( const bool aPopulateTopology, const bool aPopulateGeo );

            void
            populate_facet_data( const bool aPopulateElements, const bool aPopulateTopology, const bool aPopulateGeo );

            void
            populate_edge_data( const bool aPopulateTopology );

            void
            populate_face_data( const bool aPopulateElements );

            void
            populate_control_point_data( const bool aPopulateTopology );

//------------------------------------------------------------------------------

            void
            reconstruct_edge_connectivity();

            //! re-point 1:1 edge-on-edge hangs to the twin the slot
            //! reconstruction favors; must run AFTER load_hanging_entities()
            void
            normalize_edge_hangs();

            void
            reconstruct_face_connectivity();

//------------------------------------------------------------------------------

            void
            reset_node_data();

            void
            reset_element_data();

            void
            reset_facet_data();

            void
            reset_edge_data();

            void
            reset_face_data();

            void
            reset_control_point_data();

            void
            reset_periodicity_data();

//------------------------------------------------------------------------------

            Block *
            block( const id_t aID );

            SideSet *
            sideset( const id_t aID );

            Node *
            node( const id_t aID );

            Edge *
            edge( const id_t aID );

            Face *
            face( const id_t aID ) ;

            Facet *
            facet( const id_t aID ) ;

            Element *
            element( const id_t aID ) ;

            ControlPoint *
            control_point( const id_t aID );

            Basis *
            basis( const EntityType aType, const id_t aID );

        };

//------------------------------------------------------------------------------

        inline proto::MetaData *
        ProtoMesh::meta_data()
        {
            return mMetaData ;
        }

        inline Cell< proto::GroupData >&
        ProtoMesh::block_data()
        {
            return mBlockData ;
        }

        inline proto::GroupData &
        ProtoMesh::block_data( const index_t aIndex )
        {
            return mBlockData( aIndex );
        }

        inline Cell<  proto::GroupData > &
        ProtoMesh::sideset_data()
        {
            return mSideSetData ;
        }

        inline proto::GroupData &
        ProtoMesh::sideset_data( const index_t aIndex )
        {
            return mSideSetData( aIndex );
        }

        inline Block *
        ProtoMesh::block( const id_t aID )
        {
            return mBlockMap( aID );
        }

        inline SideSet *
        ProtoMesh::sideset( const id_t aID )
        {
            return mSideSetMap( aID );
        }

        inline Cell< proto::ThinShellData > &
        ProtoMesh::thin_shell_data()
        {
            return mThinShellData ;
        }

        inline proto::ThinShellData &
        ProtoMesh::thin_shell_data( const index_t aIndex )
        {
            return mThinShellData( aIndex );
        }


        inline proto::NodeData *
        ProtoMesh::node_data()
        {
            return mNodeData ;
        }

        inline proto::ElementData *
        ProtoMesh::element_data()
        {
            return mElementData ;
        }

        inline proto::EdgeData *
        ProtoMesh::edge_data()
        {
            return mEdgeData ;
        }

        inline proto::FaceData *
        ProtoMesh::face_data()
        {
            return mFaceData ;
        }

        inline proto::FacetData *
        ProtoMesh::facet_data()
        {
            return mFacetData ;
        }

        inline proto::ElementExtra *
        ProtoMesh::element_extra()
        {
            return mElementExtra ;
        }

        inline proto::FacetExtra *
        ProtoMesh::facet_extra()
        {
            return mFacetExtra ;
        }

        inline proto::ElementData *
        ProtoMesh::vertex_data()
        {
            return mVertexData ;
        }

        inline proto::ControlPointData *
        ProtoMesh::control_point_data()
        {
            return mControlPointData ;
        }

        inline proto::PeriodicityData *
        ProtoMesh::periodicity_data()
        {
            return mPeriodicityData ;
        }

        inline proto::TMatrixData *
        ProtoMesh::t_matrix_data()
        {
            return mTMatrixData ;
        }

        inline Node *
        ProtoMesh::node( const id_t aID )
        {
            return mNodeMap( aID );
        }

        inline Edge *
        ProtoMesh::edge( const id_t aID )
        {
            return mEdgeMap( aID );
        }

        inline Face *
        ProtoMesh::face( const id_t aID )
        {
            return mFaceMap( aID );
        }

        inline Facet *
        ProtoMesh::facet( const id_t aID )
        {
            return mFacetMap( aID );
        }

        inline Element *
        ProtoMesh::element( const id_t aID )
        {
            return mElementMap( aID );
        }

        inline ControlPoint *
        ProtoMesh::control_point( const id_t aID )
        {
            return mControlPointMap( aID );
        }

        inline Basis *
        ProtoMesh::basis( const EntityType aType, const id_t aID )
        {
            switch ( aType )
            {
                case EntityType::NODE :
                {
                    return mNodeMap( aID );
                }
                case EntityType::EDGE :
                {
                    return mEdgeMap( aID );
                }
                case EntityType::FACE :
                {
                    return mFaceMap( aID );
                }
                case EntityType::FACET :
                {
                    return mFacetMap( aID );
                }
                case EntityType::ELEMENT :
                {
                    return mElementMap( aID );
                }
                case EntityType::CONTROLPOINT :
                {
                    return mControlPointMap( aID );
                }
                default:
                {
                    BELFEM_ERROR( false, "unsupported entity type" );
                    return nullptr ;
                }
            }
        }

    }
}
#endif //BELFEM_CL_MESHDATA_HPP